using Dates
using Printf
using Statistics

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))

include(joinpath(@__DIR__, "write_ecrad_candidates.jl"))

const ALL_SKY_REFERENCE = "validation/reference/ecrad/ecckd_all_sky_tropical_column.nc"
const ECRAD_ALL_SKY_PROPERTIES =
    "validation/external/ecrad/test/ifs/radiative_properties_ecckd_tc.nc"

const ALL_SKY_OPTICS_ENV = Dict(
    "RH_CANDIDATE_GAS_OPTICS" => "official_ecckd",
    "RH_CLOUD_EFFECTIVE_RADIUS_OPTICS" => "true",
    "RH_LIQUID_CLOUD_LW_MASS_ABSORPTION" => "100",
    "RH_ICE_CLOUD_LW_MASS_ABSORPTION" => "50",
    "RH_LIQUID_CLOUD_SW_SINGLE_SCATTERING_ALBEDO" => "1.0",
    "RH_ICE_CLOUD_SW_SINGLE_SCATTERING_ALBEDO" => "1.0",
    "RH_LIQUID_CLOUD_SW_SCATTERING_ASYMMETRY" => "0.85",
    "RH_ICE_CLOUD_SW_SCATTERING_ASYMMETRY" => "0.75",
    "RH_CLOUD_FRACTION_EXPONENT" => "0.5",
    "RH_LIQUID_CLOUD_EFFECTIVE_RADIUS_SW_SCALE" => "1.0",
    "RH_ICE_CLOUD_EFFECTIVE_RADIUS_SW_SCALE" => "1.0",
    "RH_CLOUD_SCATTERING_TABLE_OPTICS" => "true",
    "RH_CLOUD_SCATTERING_MAPPING_METHOD" => "ecrad",
    "RH_CLOUD_SCATTERING_DELTA_EDDINGTON_AVERAGE" => "true",
    "RH_CLOUD_SCATTERING_THICK_AVERAGING" => "true",
    "RH_CLOUD_SCATTERING_DELTA_EDDINGTON_SCALE" => "true",
    "RH_CLOUD_SCATTERING_SW_SSA_SCALE" => "1.0",
    "RH_CLOUD_LW_TABLE_EXTINCTION_AS_ABSORPTION" => "false",
    "RH_CLOUD_OVERLAP_SHORTWAVE" => "false",
    "RH_CLOUD_OVERLAP_RULE" => "maximum",
    "RH_AEROSOL_OPTICS" => "false",
)

function with_all_sky_optics_env(f)
    previous = Dict(key => get(ENV, key, nothing) for key in keys(ALL_SKY_OPTICS_ENV))
    try
        for (key, value) in ALL_SKY_OPTICS_ENV
            ENV[key] = value
        end
        return f()
    finally
        for (key, value) in previous
            if value === nothing
                delete!(ENV, key)
            else
                ENV[key] = value
            end
        end
    end
end

function compare_stats(candidate, reference)
    error = vec(candidate .- reference)
    return (
        rmse = sqrt(mean(abs2, error)),
        max_abs = maximum(abs, error),
        mean_bias = mean(error),
        mean_abs = mean(abs.(error)),
        reference_mean = mean(vec(reference)),
        candidate_mean = mean(vec(candidate)),
    )
end

function total_sw_optics(shortwave)
    total = shortwave.optical_depth .+ shortwave.rayleigh_optical_depth
    ssa = similar(total)
    asymmetry = similar(total)
    for index in eachindex(total)
        ssa[index] = total[index] == 0 ? 0 : shortwave.rayleigh_optical_depth[index] / total[index]
        asymmetry[index] = shortwave.scattering_asymmetry[index]
    end
    return total, ssa, asymmetry
end

function cloudy_region_cloud_arrays(dataset, column, ng_lw, ng_sw)
    variables = String.(collect(keys(dataset)))
    nlayers = size(dataset["pressure_layer"], 1)
    lw_cloud = zeros(Float64, ng_lw, nlayers)
    lw_scat_cloud = zeros(Float64, ng_lw, nlayers)
    lw_asymmetry_cloud = zeros(Float64, ng_lw, nlayers)
    sw_abs_cloud = zeros(Float64, ng_sw, nlayers)
    sw_scat_cloud = zeros(Float64, ng_sw, nlayers)
    sw_asymmetry_cloud = zeros(Float64, ng_sw, nlayers)
    all(name in variables for name in ("cloud_fraction", "liquid_water_path", "ice_water_path")) ||
        return lw_cloud, lw_scat_cloud, lw_asymmetry_cloud,
               sw_abs_cloud, sw_scat_cloud, sw_asymmetry_cloud

    liquid_water_path = Array(dataset["liquid_water_path"][:, column])
    ice_water_path = Array(dataset["ice_water_path"][:, column])
    cloud_fraction = Array(dataset["cloud_fraction"][:, column])
    liquid_lw_absorption =
        env_float("RH_LIQUID_CLOUD_LW_MASS_ABSORPTION",
                  env_float("RH_CLOUD_LW_MASS_ABSORPTION", 0.0))
    ice_lw_absorption =
        env_float("RH_ICE_CLOUD_LW_MASS_ABSORPTION",
                  env_float("RH_CLOUD_LW_MASS_ABSORPTION", 0.0))
    liquid_ssa = clamp(env_float("RH_LIQUID_CLOUD_SW_SINGLE_SCATTERING_ALBEDO",
                                 env_float("RH_CLOUD_SW_SINGLE_SCATTERING_ALBEDO", 0.0)),
                       0.0, 1.0)
    ice_ssa = clamp(env_float("RH_ICE_CLOUD_SW_SINGLE_SCATTERING_ALBEDO",
                              env_float("RH_CLOUD_SW_SINGLE_SCATTERING_ALBEDO", 0.0)),
                    0.0, 1.0)
    liquid_asymmetry =
        clamp(env_float("RH_LIQUID_CLOUD_SW_SCATTERING_ASYMMETRY",
                        env_float("RH_CLOUD_SW_SCATTERING_ASYMMETRY", 0.0)),
              -1.0, 1.0)
    ice_asymmetry =
        clamp(env_float("RH_ICE_CLOUD_SW_SCATTERING_ASYMMETRY",
                        env_float("RH_CLOUD_SW_SCATTERING_ASYMMETRY", 0.0)),
              -1.0, 1.0)

    if env_bool("RH_CLOUD_SCATTERING_TABLE_OPTICS", false)
        mapped = mapped_cloud_scattering_properties()
        re_liquid = "re_liquid" in variables ?
            Array(dataset["re_liquid"][:, column]) :
            fill(mapped.liquid_effective_radius, nlayers)
        re_ice = "re_ice" in variables ?
            Array(dataset["re_ice"][:, column]) :
            fill(mapped.ice_effective_radius, nlayers)
        liquid_scale = env_float("RH_LIQUID_CLOUD_SCATTERING_SW_SCALE", 1.0)
        ice_scale = env_float("RH_ICE_CLOUD_SCATTERING_SW_SCALE", 1.0)
        sw_ssa_scale = max(env_float("RH_CLOUD_SCATTERING_SW_SSA_SCALE", 1.0), 0.0)
        liquid_lw_scale = env_float("RH_LIQUID_CLOUD_SCATTERING_LW_SCALE", 1.0)
        ice_lw_scale = env_float("RH_ICE_CLOUD_SCATTERING_LW_SCALE", 1.0)
        for k in 1:nlayers
            fraction = clamp(cloud_fraction[k], 0.0, 1.0)
            path_scale = fraction > 0 ? inv(fraction) : 0.0
            cloudy_liquid_path = max(liquid_water_path[k], 0.0) * path_scale
            cloudy_ice_path = max(ice_water_path[k], 0.0) * path_scale
            liquid_lw = cloud_scattering_gpoint_properties(
                mapped.liquid_table, mapped.lw_mapping, re_liquid[k];
                mapping_method = mapped.mapping_method,
                delta_eddington_average = mapped.delta_eddington_average,
                thick_averaging = mapped.thick_averaging)
            ice_lw = cloud_scattering_gpoint_properties(
                mapped.ice_table, mapped.lw_mapping, re_ice[k];
                mapping_method = mapped.mapping_method,
                delta_eddington_average = mapped.delta_eddington_average,
                thick_averaging = mapped.thick_averaging)
            liquid_sw = cloud_scattering_gpoint_properties(
                mapped.liquid_table, mapped.sw_mapping, re_liquid[k];
                mapping_method = mapped.mapping_method,
                delta_eddington_average = mapped.delta_eddington_average,
                thick_averaging = mapped.thick_averaging)
            ice_sw = cloud_scattering_gpoint_properties(
                mapped.ice_table, mapped.sw_mapping, re_ice[k];
                mapping_method = mapped.mapping_method,
                delta_eddington_average = mapped.delta_eddington_average,
                thick_averaging = mapped.thick_averaging)
            for ig in 1:ng_lw
                liquid_extinction =
                    liquid_lw_scale * liquid_lw.mass_extinction_coefficient[ig] *
                    cloudy_liquid_path
                ice_extinction =
                    ice_lw_scale * ice_lw.mass_extinction_coefficient[ig] *
                    cloudy_ice_path
                liquid_scat =
                    liquid_lw.single_scattering_albedo[ig] * liquid_extinction
                ice_scat =
                    ice_lw.single_scattering_albedo[ig] * ice_extinction
                total_scat = liquid_scat + ice_scat
                total_extinction = liquid_extinction + ice_extinction
                lw_cloud[ig, k] = total_extinction
                lw_scat_cloud[ig, k] = total_scat
                lw_asymmetry_cloud[ig, k] = total_scat == 0 ? 0.0 :
                    (liquid_lw.asymmetry_factor[ig] * liquid_scat +
                     ice_lw.asymmetry_factor[ig] * ice_scat) / total_scat
            end
            for ig in 1:ng_sw
                liquid_extinction =
                    liquid_scale * liquid_sw.mass_extinction_coefficient[ig] *
                    cloudy_liquid_path
                ice_extinction =
                    ice_scale * ice_sw.mass_extinction_coefficient[ig] *
                    cloudy_ice_path
                liquid_scat =
                    liquid_sw.single_scattering_albedo[ig] * liquid_extinction
                ice_scat =
                    ice_sw.single_scattering_albedo[ig] * ice_extinction
                total_scat = liquid_scat + ice_scat
                total_extinction = liquid_extinction + ice_extinction
                asymmetry = total_scat == 0 ? 0.0 :
                    (liquid_sw.asymmetry_factor[ig] * liquid_scat +
                     ice_sw.asymmetry_factor[ig] * ice_scat) / total_scat
                if env_bool("RH_CLOUD_SCATTERING_DELTA_EDDINGTON_SCALE", true) &&
                   total_scat > 0
                    forward_fraction = asymmetry^2
                    total_extinction -= total_scat * forward_fraction
                    total_scat *= 1.0 - forward_fraction
                    asymmetry /= 1.0 + asymmetry
                end
                if sw_ssa_scale != 1.0
                    total_scat = min(total_scat * sw_ssa_scale,
                                     max(total_extinction, 0.0))
                end
                sw_abs_cloud[ig, k] = max(total_extinction - total_scat, 0.0)
                sw_scat_cloud[ig, k] = total_scat
                sw_asymmetry_cloud[ig, k] = asymmetry
            end
        end
        return lw_cloud, lw_scat_cloud, lw_asymmetry_cloud,
               sw_abs_cloud, sw_scat_cloud, sw_asymmetry_cloud
    elseif env_bool("RH_CLOUD_EFFECTIVE_RADIUS_OPTICS", false) &&
       all(name in variables for name in ("re_liquid", "re_ice"))
        re_liquid = Array(dataset["re_liquid"][:, column])
        re_ice = Array(dataset["re_ice"][:, column])
        liquid_sw_scale = env_float("RH_LIQUID_CLOUD_EFFECTIVE_RADIUS_SW_SCALE", 1.0)
        ice_sw_scale = env_float("RH_ICE_CLOUD_EFFECTIVE_RADIUS_SW_SCALE", 1.0)
        for k in 1:nlayers
            fraction = clamp(cloud_fraction[k], 0.0, 1.0)
            path_scale = fraction > 0 ? inv(fraction) : 0.0
            cloudy_liquid_path = max(liquid_water_path[k], 0.0) * path_scale
            cloudy_ice_path = max(ice_water_path[k], 0.0) * path_scale
            liquid_extinction = liquid_sw_scale * 1.5 * cloudy_liquid_path /
                (LIQUID_WATER_DENSITY * max(re_liquid[k], eps(Float64)))
            ice_extinction = ice_sw_scale * 1.5 * cloudy_ice_path /
                (ICE_DENSITY * max(re_ice[k], eps(Float64)))
            liquid_scat = liquid_ssa * liquid_extinction
            ice_scat = ice_ssa * ice_extinction
            total_scat = liquid_scat + ice_scat
            lw = liquid_lw_absorption * cloudy_liquid_path +
                 ice_lw_absorption * cloudy_ice_path
            sw_abs = (1.0 - liquid_ssa) * liquid_extinction +
                     (1.0 - ice_ssa) * ice_extinction
            sw_scat = total_scat
            asymmetry = total_scat == 0 ? 0.0 :
                (liquid_asymmetry * liquid_scat + ice_asymmetry * ice_scat) /
                total_scat
            lw_cloud[:, k] .= lw
            sw_abs_cloud[:, k] .= sw_abs
            sw_scat_cloud[:, k] .= sw_scat
            sw_asymmetry_cloud[:, k] .= asymmetry
        end
        return lw_cloud, lw_scat_cloud, lw_asymmetry_cloud,
               sw_abs_cloud, sw_scat_cloud, sw_asymmetry_cloud
    end

    liquid_sw_extinction =
        env_float("RH_LIQUID_CLOUD_SW_MASS_EXTINCTION",
                  env_float("RH_CLOUD_SW_MASS_EXTINCTION", 0.0))
    ice_sw_extinction =
        env_float("RH_ICE_CLOUD_SW_MASS_EXTINCTION",
                  env_float("RH_CLOUD_SW_MASS_EXTINCTION", 0.0))
    for k in 1:nlayers
        fraction = clamp(cloud_fraction[k], 0.0, 1.0)
        path_scale = fraction > 0 ? inv(fraction) : 0.0
        cloudy_liquid_path = max(liquid_water_path[k], 0.0) * path_scale
        cloudy_ice_path = max(ice_water_path[k], 0.0) * path_scale
        liquid_extinction = liquid_sw_extinction * cloudy_liquid_path
        ice_extinction = ice_sw_extinction * cloudy_ice_path
        liquid_scat = liquid_ssa * liquid_extinction
        ice_scat = ice_ssa * ice_extinction
        total_scat = liquid_scat + ice_scat
        lw_cloud[:, k] .= liquid_lw_absorption * cloudy_liquid_path +
                          ice_lw_absorption * cloudy_ice_path
        sw_abs_cloud[:, k] .= (1.0 - liquid_ssa) * liquid_extinction +
                              (1.0 - ice_ssa) * ice_extinction
        sw_scat_cloud[:, k] .= total_scat
        sw_asymmetry_cloud[:, k] .= total_scat == 0 ? 0.0 :
            (liquid_asymmetry * liquid_scat + ice_asymmetry * ice_scat) /
            total_scat
    end
    return lw_cloud, lw_scat_cloud, lw_asymmetry_cloud,
           sw_abs_cloud, sw_scat_cloud, sw_asymmetry_cloud
end

function candidate_all_sky_optics(reference)
    gas_optics = candidate_gas_optics(Float64)
    pressure_layers = Array(reference["pressure_layer"])
    pressure_interfaces = Array(reference["pressure_interface"])
    temperature_layers = Array(reference["temperature_layer"])
    temperature_interfaces = Array(reference["temperature_interface"])
    gas_amounts = gas_column_amounts(reference, pressure_interfaces)
    variables = String.(collect(keys(reference)))
    solar_irradiance = "solar_irradiance" in variables ?
        Array(reference["solar_irradiance"]) :
        fill(DEFAULT_SOLAR_IRRADIANCE, size(pressure_layers, 2))
    toa_shortwave_down = Array(reference["sw_down"])[1, :]
    cos_zenith = "cos_solar_zenith_angle" in variables ?
        Array(reference["cos_solar_zenith_angle"]) :
        clamp.(toa_shortwave_down ./ solar_irradiance, 0.0, 1.0)

    nlayers, ncolumns = size(pressure_layers)
    ng_lw = size(gas_optics.longwave_absorption, 1)
    ng_sw = size(gas_optics.shortwave_absorption, 1)

    lw_total = zeros(Float64, ng_lw, nlayers, ncolumns)
    lw_cloud = zeros(Float64, ng_lw, nlayers, ncolumns)
    sw_total = zeros(Float64, ng_sw, nlayers, ncolumns)
    sw_ssa = zeros(Float64, ng_sw, nlayers, ncolumns)
    sw_asymmetry = zeros(Float64, ng_sw, nlayers, ncolumns)
    sw_cloud = zeros(Float64, ng_sw, nlayers, ncolumns)
    sw_cloud_ssa = zeros(Float64, ng_sw, nlayers, ncolumns)
    sw_cloud_asymmetry = zeros(Float64, ng_sw, nlayers, ncolumns)
    cloudy_region_lw_cloud = zeros(Float64, ng_lw, nlayers, ncolumns)
    cloudy_region_lw_cloud_ssa = zeros(Float64, ng_lw, nlayers, ncolumns)
    cloudy_region_lw_cloud_asymmetry = zeros(Float64, ng_lw, nlayers, ncolumns)
    cloudy_region_sw_cloud = zeros(Float64, ng_sw, nlayers, ncolumns)
    cloudy_region_sw_cloud_ssa = zeros(Float64, ng_sw, nlayers, ncolumns)
    cloudy_region_sw_cloud_asymmetry = zeros(Float64, ng_sw, nlayers, ncolumns)
    cloudy_region_od_lw = zeros(Float64, ng_lw, nlayers, ncolumns)
    cloudy_region_od_sw = zeros(Float64, ng_sw, nlayers, ncolumns)
    cloudy_region_ssa_sw = zeros(Float64, ng_sw, nlayers, ncolumns)
    cloudy_region_asymmetry_sw = zeros(Float64, ng_sw, nlayers, ncolumns)

    for j in 1:ncolumns
        atmosphere = ColumnAtmosphere(
            pressure_layers = pressure_layers[:, j],
            pressure_interfaces = pressure_interfaces[:, j],
            temperature_layers = temperature_layers[:, j],
            temperature_interfaces = temperature_interfaces[:, j],
            gases = Dict(name => values[:, j] for (name, values) in gas_amounts.amounts),
            surface = (;),
            geometry = (; cos_zenith = cos_zenith[j]),
        )
        longwave = LongwaveOpticalProperties(zeros(Float64, ng_lw, nlayers),
                                             zeros(Float64, ng_lw, nlayers))
        shortwave = ShortwaveOpticalProperties(zeros(Float64, ng_sw, nlayers))
        optical_properties!(longwave, shortwave, gas_optics, atmosphere)
        gas_lw = copy(longwave.optical_depth)
        gas_sw_abs = copy(shortwave.optical_depth)
        gas_sw_scat = copy(shortwave.rayleigh_optical_depth)
        gas_sw_total = gas_sw_abs .+ gas_sw_scat

        region_lw, region_lw_scat, region_lw_asymmetry,
            region_sw_abs, region_sw_scat, region_sw_asymmetry =
            cloudy_region_cloud_arrays(reference, j, ng_lw, ng_sw)
        region_total = gas_sw_total .+ region_sw_abs .+ region_sw_scat
        region_ssa = similar(region_total)
        region_asymmetry = similar(region_total)
        region_cloud_total = region_sw_abs .+ region_sw_scat
        region_cloud_ssa = similar(region_total)
        region_cloud_asymmetry = similar(region_total)
        for index in eachindex(region_total)
            region_ssa[index] =
                region_total[index] == 0 ? 0 :
                (gas_sw_scat[index] + region_sw_scat[index]) / region_total[index]
            region_asymmetry[index] =
                region_sw_scat[index] == 0 ? 0 : region_sw_asymmetry[index]
            region_cloud_ssa[index] =
                region_cloud_total[index] == 0 ? 0 :
                region_sw_scat[index] / region_cloud_total[index]
            region_cloud_asymmetry[index] =
                region_sw_scat[index] == 0 ? 0 : region_sw_asymmetry[index]
        end
        cloudy_region_od_lw[:, :, j] = gas_lw .+ region_lw
        cloudy_region_lw_cloud[:, :, j] = region_lw
        region_lw_ssa = similar(region_lw)
        for index in eachindex(region_lw)
            region_lw_ssa[index] =
                region_lw[index] == 0 ? 0.0 : region_lw_scat[index] / region_lw[index]
        end
        cloudy_region_lw_cloud_ssa[:, :, j] = region_lw_ssa
        cloudy_region_lw_cloud_asymmetry[:, :, j] = region_lw_asymmetry
        cloudy_region_sw_cloud[:, :, j] = region_cloud_total
        cloudy_region_sw_cloud_ssa[:, :, j] = region_cloud_ssa
        cloudy_region_sw_cloud_asymmetry[:, :, j] = region_cloud_asymmetry
        cloudy_region_od_sw[:, :, j] = region_total
        cloudy_region_ssa_sw[:, :, j] = region_ssa
        cloudy_region_asymmetry_sw[:, :, j] = region_asymmetry

        maybe_add_cloud_optics!(longwave, shortwave, reference, j)
        maybe_add_aerosol_optics!(longwave, shortwave, reference, pressure_interfaces, j)

        total, ssa, asymmetry = total_sw_optics(shortwave)
        cloud_total = total .- gas_sw_abs .- gas_sw_scat
        cloud_scat = shortwave.rayleigh_optical_depth .- gas_sw_scat
        cloud_ssa = similar(cloud_total)
        cloud_asymmetry = similar(cloud_total)
        for index in eachindex(cloud_total)
            cloud_ssa[index] = cloud_total[index] == 0 ? 0 : cloud_scat[index] / cloud_total[index]
            cloud_asymmetry[index] = shortwave.scattering_asymmetry[index]
        end

        lw_total[:, :, j] = longwave.optical_depth
        lw_cloud[:, :, j] = longwave.optical_depth .- gas_lw
        sw_total[:, :, j] = total
        sw_ssa[:, :, j] = ssa
        sw_asymmetry[:, :, j] = asymmetry
        sw_cloud[:, :, j] = cloud_total
        sw_cloud_ssa[:, :, j] = cloud_ssa
        sw_cloud_asymmetry[:, :, j] = cloud_asymmetry
    end

    return (
        od_lw = lw_total,
        od_lw_cloud = lw_cloud,
        od_sw = sw_total,
        ssa_sw = sw_ssa,
        asymmetry_sw = sw_asymmetry,
        od_sw_cloud = sw_cloud,
        ssa_sw_cloud = sw_cloud_ssa,
        asymmetry_sw_cloud = sw_cloud_asymmetry,
        cloudy_region_od_lw = cloudy_region_od_lw,
        cloudy_region_od_lw_cloud = cloudy_region_lw_cloud,
        cloudy_region_ssa_lw_cloud = cloudy_region_lw_cloud_ssa,
        cloudy_region_asymmetry_lw_cloud = cloudy_region_lw_cloud_asymmetry,
        cloudy_region_od_sw = cloudy_region_od_sw,
        cloudy_region_ssa_sw = cloudy_region_ssa_sw,
        cloudy_region_asymmetry_sw = cloudy_region_asymmetry_sw,
        cloudy_region_od_sw_cloud = cloudy_region_sw_cloud,
        cloudy_region_ssa_sw_cloud = cloudy_region_sw_cloud_ssa,
        cloudy_region_asymmetry_sw_cloud = cloudy_region_sw_cloud_asymmetry,
    )
end

function push_optics_gap_rows!(rows, properties, candidate, original_columns, candidate_kind, mapping)
    for variable in mapping
        haskey(properties, variable.reference) || continue
        reference_values = Array(properties[variable.reference][:, :, original_columns])
        candidate_values = getproperty(candidate, Symbol(variable.candidate))
        stats = compare_stats(candidate_values, reference_values)
        push!(rows, (
            candidate_kind = candidate_kind,
            variable = variable.reference,
            candidate_variable = variable.candidate,
            units = variable.units,
            rmse = stats.rmse,
            max_abs = stats.max_abs,
            mean_bias = stats.mean_bias,
            mean_abs = stats.mean_abs,
            reference_mean = stats.reference_mean,
            candidate_mean = stats.candidate_mean,
        ))
    end
    return rows
end

function optics_gap_rows(reference, properties, candidate)
    original_columns = round.(Int, Array(reference["column"]))
    rows = NamedTuple[]
    grid_mean_mapping = (
        (reference = "od_lw", candidate = "od_lw", units = "1"),
        (reference = "od_lw_cloud", candidate = "od_lw_cloud", units = "1"),
        (reference = "od_sw", candidate = "od_sw", units = "1"),
        (reference = "ssa_sw", candidate = "ssa_sw", units = "1"),
        (reference = "asymmetry_sw", candidate = "asymmetry_sw", units = "1"),
        (reference = "od_sw_cloud", candidate = "od_sw_cloud", units = "1"),
        (reference = "ssa_sw_cloud", candidate = "ssa_sw_cloud", units = "1"),
        (reference = "asymmetry_sw_cloud", candidate = "asymmetry_sw_cloud", units = "1"),
    )
    cloudy_region_mapping = (
        (reference = "od_lw", candidate = "cloudy_region_od_lw", units = "1"),
        (reference = "od_lw_cloud", candidate = "cloudy_region_od_lw_cloud", units = "1"),
        (reference = "ssa_lw_cloud", candidate = "cloudy_region_ssa_lw_cloud", units = "1"),
        (reference = "asymmetry_lw_cloud", candidate = "cloudy_region_asymmetry_lw_cloud", units = "1"),
        (reference = "od_sw", candidate = "cloudy_region_od_sw", units = "1"),
        (reference = "ssa_sw", candidate = "cloudy_region_ssa_sw", units = "1"),
        (reference = "asymmetry_sw", candidate = "cloudy_region_asymmetry_sw", units = "1"),
        (reference = "od_sw_cloud", candidate = "cloudy_region_od_sw_cloud", units = "1"),
        (reference = "ssa_sw_cloud", candidate = "cloudy_region_ssa_sw_cloud", units = "1"),
        (reference = "asymmetry_sw_cloud", candidate = "cloudy_region_asymmetry_sw_cloud", units = "1"),
    )
    push_optics_gap_rows!(rows, properties, candidate, original_columns,
                          "grid_mean_current", grid_mean_mapping)
    cloudy_region_kind = env_bool("RH_CLOUD_SCATTERING_DELTA_EDDINGTON_SCALE", true) ?
        "cloudy_region_delta_scaled" : "cloudy_region_unscaled"
    push_optics_gap_rows!(rows, properties, candidate, original_columns,
                          cloudy_region_kind, cloudy_region_mapping)
    return rows
end

function run_all_sky_optics_gap()
    return with_all_sky_optics_env() do
        NCDATASETS.NCDataset(reference_path(ALL_SKY_REFERENCE)) do reference
            NCDATASETS.NCDataset(reference_path(ECRAD_ALL_SKY_PROPERTIES)) do properties
                candidate = candidate_all_sky_optics(reference)
                rows = optics_gap_rows(reference, properties, candidate)
                return (
                    case = "ecrad_all_sky_optics_gap",
                    date = string(Dates.now()),
                    reference_case = ALL_SKY_REFERENCE,
                    ecrad_properties = ECRAD_ALL_SKY_PROPERTIES,
                    candidate_configuration = [
                        (variable = key, value = ALL_SKY_OPTICS_ENV[key])
                        for key in sort(collect(keys(ALL_SKY_OPTICS_ENV)))
                    ],
                    rows = rows,
                )
            end
        end
    end
end

function markdown_all_sky_optics_gap(result)
    lines = String[
        "# ecRad All-Sky Optics Gap",
        "",
        "This diagnostic compares candidate optical properties against ecRad's saved all-sky `radiative_properties_ecckd_tc.nc` before radiative-transfer solving.",
        "",
        "| Candidate kind | Variable | RMSE | Max abs | Mean bias | Mean abs | Reference mean | Candidate mean |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in result.rows
        push!(lines, "| `$(row.candidate_kind)` | `$(row.variable)` | $(@sprintf("%.12g", row.rmse)) | $(@sprintf("%.12g", row.max_abs)) | $(@sprintf("%.12g", row.mean_bias)) | $(@sprintf("%.12g", row.mean_abs)) | $(@sprintf("%.12g", row.reference_mean)) | $(@sprintf("%.12g", row.candidate_mean)) |")
    end
    append!(lines, [
        "",
        "## Candidate Configuration",
        "",
        "| Environment variable | Value |",
        "|---|---|",
    ])
    for row in result.candidate_configuration
        push!(lines, "| `$(row.variable)` | `$(row.value)` |")
    end
    return join(lines, "\n") * "\n"
end

function all_sky_optics_gap_main()
    result = run_all_sky_optics_gap()
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)
    json_path = joinpath(results_dir, "ecrad_all_sky_optics_gap.json")
    md_path = joinpath(results_dir, "ecrad_all_sky_optics_gap.md")
    write(json_path, json_object(result))
    write(md_path, markdown_all_sky_optics_gap(result))
    print(markdown_all_sky_optics_gap(result))
    println("Wrote $json_path")
    println("Wrote $md_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    all_sky_optics_gap_main()
end
