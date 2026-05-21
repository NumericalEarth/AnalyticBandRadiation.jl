using Dates
using Printf

using Lightflux

include(joinpath(@__DIR__, "ecrad_reference_manifest.jl"))

const CANDIDATE_PREFIX = "radiative_heating_"
const GRAVITY = 9.80665
const MOLAR_MASS_DRY_AIR = 0.0289647
const MOLAR_MASS_WATER = 0.01801528
const MOLAR_MASS_OZONE = 0.0479982
const DEFAULT_SOLAR_IRRADIANCE = 1361.0
const LIQUID_WATER_DENSITY = 1000.0
const ICE_DENSITY = 917.0
const WATER_VAPOR_EPSILON = 0.622
const OFFICIAL_ECCKD_GASES = (:composite, :h2o, :o3, :co2, :ch4, :n2o, :cfc11, :cfc12)
const CLOUD_SCATTERING_CACHE = Ref{Any}(nothing)
const LAYER_CLOUD_SCATTERING_CACHE = Ref(Dict{Any, Any}())
const IFS_AEROSOL_CACHE = Ref{Any}(nothing)
const IFS_AEROSOL_TYPE_MAP = (-1, -2, -3, 7, 8, 9, -4, 10, 11, 11, -5, 14)

env_float(name, default) = parse(Float64, get(ENV, name, string(default)))
env_bool(name, default) = lowercase(get(ENV, name, default ? "true" : "false")) in
    ("1", "true", "yes", "on")
candidate_mode() = get(ENV, "RH_CANDIDATE_GAS_OPTICS", "toy")

function require_ncdatasets()
    NCDATASETS === nothing &&
        error("NCDatasets is required. Run with `julia --project=test validation/write_ecrad_candidates.jl`.")
    return NCDATASETS
end

function candidate_gas_optics(FT)
    mode = candidate_mode()
    if mode == "official_ecckd"
        lw_path = reference_path(get(ENV, "RH_ECCKD_LW_PATH",
                                     "validation/external/ecrad/data/ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc"))
        sw_path = reference_path(get(ENV, "RH_ECCKD_SW_PATH",
                                     "validation/external/ecrad/data/ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc"))
        return read_ecckd_tabulated_gas_optics(lw_path, sw_path;
                                               gas_names = OFFICIAL_ECCKD_GASES,
                                               h2o_mole_fraction = env_float("RH_ECCKD_H2O_MOLE_FRACTION", 0.005))
    elseif mode != "toy"
        error("unsupported RH_CANDIDATE_GAS_OPTICS=$mode; expected `toy` or `official_ecckd`")
    end

    ng_lw = 32
    ng_sw = 16
    lw_h2o_scale = FT(env_float("RH_LW_H2O_SCALE", 20.0))
    lw_co2_scale = FT(env_float("RH_LW_CO2_SCALE", 20.0))
    sw_h2o_scale = FT(env_float("RH_SW_H2O_SCALE", 300.0))
    sw_co2_scale = FT(env_float("RH_SW_CO2_SCALE", 300.0))
    longwave_absorption = zeros(FT, ng_lw, 2)
    shortwave_absorption = zeros(FT, ng_sw, 2)

    for ig in 1:ng_lw
        longwave_absorption[ig, 1] = lw_h2o_scale * FT(3.0e-4 / sqrt(ig))
        longwave_absorption[ig, 2] = lw_co2_scale * FT(2.0e-2 + 8.0e-4ig)
    end
    for ig in 1:ng_sw
        shortwave_absorption[ig, 1] = sw_h2o_scale * FT(4.0e-5 / sqrt(ig))
        shortwave_absorption[ig, 2] = sw_co2_scale * FT(1.0e-5 + 2.0e-7ig)
    end

    return EcCKDGasOpticsModel(
        gas_names = (:h2o, :co2),
        longwave_absorption = longwave_absorption,
        shortwave_absorption = shortwave_absorption,
        longwave_source_scale = ones(FT, ng_lw),
        longwave_weights = fill(inv(FT(ng_lw)), ng_lw),
        shortwave_weights = fill(inv(FT(ng_sw)), ng_sw),
    )
end

function layer_amounts(mixing_ratio, pressure_interfaces)
    nlayers, ncolumns = size(mixing_ratio)
    amounts = zeros(Float64, nlayers, ncolumns)
    for j in 1:ncolumns, k in 1:nlayers
        Δp = pressure_interfaces[k + 1, j] - pressure_interfaces[k, j]
        amounts[k, j] = mixing_ratio[k, j] * Δp / GRAVITY
    end
    return amounts
end

function layer_moles_from_mass_mixing_ratio(mass_mixing_ratio, pressure_interfaces, molar_mass)
    nlayers, ncolumns = size(mass_mixing_ratio)
    amounts = zeros(Float64, nlayers, ncolumns)
    for j in 1:ncolumns, k in 1:nlayers
        Δp = pressure_interfaces[k + 1, j] - pressure_interfaces[k, j]
        amounts[k, j] = mass_mixing_ratio[k, j] * Δp / GRAVITY / molar_mass
    end
    return amounts
end

layer_moles_from_specific_humidity(specific_humidity, pressure_interfaces) =
    layer_moles_from_mass_mixing_ratio(specific_humidity, pressure_interfaces, MOLAR_MASS_WATER)

function layer_air_moles(pressure_interfaces)
    nlayers, ncolumns = size(pressure_interfaces, 1) - 1, size(pressure_interfaces, 2)
    amounts = zeros(Float64, nlayers, ncolumns)
    for j in 1:ncolumns, k in 1:nlayers
        Δp = pressure_interfaces[k + 1, j] - pressure_interfaces[k, j]
        amounts[k, j] = Δp / GRAVITY / MOLAR_MASS_DRY_AIR
    end
    return amounts
end

function layer_moles_from_vmr(vmr, dry_air_moles)
    nlayers, ncolumns = size(vmr)
    amounts = zeros(Float64, nlayers, ncolumns)
    for j in 1:ncolumns, k in 1:nlayers
        amounts[k, j] = vmr[k, j] * dry_air_moles[k, j]
    end
    return amounts
end

function gas_column_amounts(dataset, pressure_interfaces)
    if candidate_mode() == "official_ecckd"
        variables = String.(collect(keys(dataset)))
        h2o = Array(dataset["h2o"])
        air_moles = layer_air_moles(pressure_interfaces)
        amounts = Dict{Symbol, Matrix{Float64}}(
            :composite => air_moles,
            :h2o => layer_moles_from_specific_humidity(h2o, pressure_interfaces),
            :co2 => layer_moles_from_vmr(Array(dataset["co2"]), air_moles),
        )
        "o3" in variables &&
            (amounts[:o3] = layer_moles_from_mass_mixing_ratio(Array(dataset["o3"]),
                                                               pressure_interfaces,
                                                               MOLAR_MASS_OZONE))
        "ch4" in variables &&
            (amounts[:ch4] = layer_moles_from_vmr(Array(dataset["ch4"]), air_moles))
        "n2o" in variables &&
            (amounts[:n2o] = layer_moles_from_vmr(Array(dataset["n2o"]), air_moles))
        "cfc11" in variables &&
            (amounts[:cfc11] = layer_moles_from_vmr(Array(dataset["cfc11"]), air_moles))
        "cfc12" in variables &&
            (amounts[:cfc12] = layer_moles_from_vmr(Array(dataset["cfc12"]), air_moles))
        return (
            amounts = amounts,
            units = "mol m^-2",
        )
    end
    return (
        amounts = Dict{Symbol, Matrix{Float64}}(
            :h2o => layer_amounts(Array(dataset["h2o"]), pressure_interfaces),
            :co2 => layer_amounts(Array(dataset["co2"]), pressure_interfaces),
        ),
        units = "kg m^-2 proxy",
    )
end

function write_variable!(dataset, name, values, dims)
    variables = String.(collect(keys(dataset)))
    if name in variables
        dataset[name][:] = values
        return dataset[name]
    end
    var = NCDATASETS.defVar(dataset, name, Float64, dims)
    var[:] = values
    return var
end

function mapped_cloud_scattering_properties()
    cached = CLOUD_SCATTERING_CACHE[]
    cached === nothing || return cached

    sw_mapping_path = reference_path(get(ENV, "RH_CLOUD_SW_MAPPING_PATH",
        "validation/external/ecrad/data/ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc"))
    lw_mapping_path = reference_path(get(ENV, "RH_CLOUD_LW_MAPPING_PATH",
        "validation/external/ecrad/data/ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc"))
    liquid_path = reference_path(get(ENV, "RH_LIQUID_CLOUD_SCATTERING_PATH",
        "validation/external/ecrad/data/mie_droplet_scattering.nc"))
    ice_path = reference_path(get(ENV, "RH_ICE_CLOUD_SCATTERING_PATH",
        "validation/external/ecrad/data/baum-general-habit-mixture_ice_scattering.nc"))

    sw_mapping = read_ecckd_spectral_mapping(sw_mapping_path)
    lw_mapping = read_ecckd_spectral_mapping(lw_mapping_path)
    liquid = read_cloud_scattering_table(liquid_path)
    ice = read_cloud_scattering_table(ice_path)
    liquid_radius = env_float("RH_LIQUID_CLOUD_SCATTERING_EFFECTIVE_RADIUS", 10.0e-6)
    ice_radius = env_float("RH_ICE_CLOUD_SCATTERING_EFFECTIVE_RADIUS", 30.0e-6)
    mapping_method = Symbol(get(ENV, "RH_CLOUD_SCATTERING_MAPPING_METHOD", "ecrad"))
    delta_average = env_bool("RH_CLOUD_SCATTERING_DELTA_EDDINGTON_AVERAGE", true)
    thick_averaging = env_bool("RH_CLOUD_SCATTERING_THICK_AVERAGING", true)
    liquid_sw = cloud_scattering_gpoint_properties(
        liquid, sw_mapping, liquid_radius;
        mapping_method = mapping_method,
        delta_eddington_average = delta_average,
        thick_averaging = thick_averaging,
    )
    ice_sw = cloud_scattering_gpoint_properties(
        ice, sw_mapping, ice_radius;
        mapping_method = mapping_method,
        delta_eddington_average = delta_average,
        thick_averaging = thick_averaging,
    )
    liquid_lw = cloud_scattering_gpoint_properties(
        liquid, lw_mapping, liquid_radius;
        mapping_method = mapping_method,
        delta_eddington_average = delta_average,
        thick_averaging = thick_averaging,
    )
    ice_lw = cloud_scattering_gpoint_properties(
        ice, lw_mapping, ice_radius;
        mapping_method = mapping_method,
        delta_eddington_average = delta_average,
        thick_averaging = thick_averaging,
    )
    cached = (
        sw_mapping = sw_mapping,
        lw_mapping = lw_mapping,
        liquid_table = liquid,
        ice_table = ice,
        liquid = liquid_sw,
        ice = ice_sw,
        liquid_sw = liquid_sw,
        ice_sw = ice_sw,
        liquid_lw = liquid_lw,
        ice_lw = ice_lw,
        liquid_effective_radius = liquid_radius,
        ice_effective_radius = ice_radius,
        mapping_method = mapping_method,
        delta_eddington_average = delta_average,
        thick_averaging = thick_averaging,
        mapping_path = sw_mapping_path,
        sw_mapping_path = sw_mapping_path,
        lw_mapping_path = lw_mapping_path,
        liquid_path = liquid_path,
        ice_path = ice_path,
    )
    CLOUD_SCATTERING_CACHE[] = cached
    return cached
end

function add_mapped_cloud_longwave_absorption!(longwave,
                                               liquid_properties,
                                               ice_properties,
                                               liquid_water_path,
                                               ice_water_path,
                                               cloud_fraction;
                                               cloud_fraction_exponent = 1,
                                               liquid_extinction_scale = 1,
                                               ice_extinction_scale = 1,
                                               extinction_as_absorption = false)
    ng, nlayers = size(longwave.optical_depth)
    length(liquid_properties.mass_extinction_coefficient) == ng ||
        throw(DimensionMismatch("liquid cloud g-point properties must match longwave g-points"))
    length(ice_properties.mass_extinction_coefficient) == ng ||
        throw(DimensionMismatch("ice cloud g-point properties must match longwave g-points"))
    length(liquid_water_path) == nlayers ||
        throw(DimensionMismatch("liquid_water_path must match longwave layers"))
    length(ice_water_path) == nlayers ||
        throw(DimensionMismatch("ice_water_path must match longwave layers"))
    length(cloud_fraction) == nlayers ||
        throw(DimensionMismatch("cloud_fraction must match longwave layers"))

    exponent = max(Float64(cloud_fraction_exponent), 0.0)
    for k in 1:nlayers
        fraction_scale = clamp(Float64(cloud_fraction[k]), 0.0, 1.0)^exponent
        lwp = max(Float64(liquid_water_path[k]), 0.0)
        iwp = max(Float64(ice_water_path[k]), 0.0)
        for ig in 1:ng
            liquid_ext = liquid_extinction_scale *
                liquid_properties.mass_extinction_coefficient[ig] * lwp
            ice_ext = ice_extinction_scale *
                ice_properties.mass_extinction_coefficient[ig] * iwp
            if extinction_as_absorption
                tau = liquid_ext + ice_ext
            else
                liquid_ssa = clamp(liquid_properties.single_scattering_albedo[ig], 0.0, 1.0)
                ice_ssa = clamp(ice_properties.single_scattering_albedo[ig], 0.0, 1.0)
                tau = (1.0 - liquid_ssa) * liquid_ext + (1.0 - ice_ssa) * ice_ext
            end
            longwave.optical_depth[ig, k] += fraction_scale * tau
        end
    end
    return longwave
end

function add_mapped_cloud_longwave_scattering!(longwave,
                                               liquid_properties,
                                               ice_properties,
                                               liquid_water_path,
                                               ice_water_path,
                                               cloud_fraction;
                                               cloud_fraction_exponent = 1,
                                               liquid_extinction_scale = 1,
                                               ice_extinction_scale = 1)
    longwave.single_scattering_albedo === nothing &&
        throw(ArgumentError("longwave single_scattering_albedo storage is required"))
    longwave.scattering_asymmetry === nothing &&
        throw(ArgumentError("longwave scattering_asymmetry storage is required"))
    ng, nlayers = size(longwave.optical_depth)
    length(liquid_properties.mass_extinction_coefficient) == ng ||
        throw(DimensionMismatch("liquid cloud g-point properties must match longwave g-points"))
    length(ice_properties.mass_extinction_coefficient) == ng ||
        throw(DimensionMismatch("ice cloud g-point properties must match longwave g-points"))
    length(liquid_water_path) == nlayers ||
        throw(DimensionMismatch("liquid_water_path must match longwave layers"))
    length(ice_water_path) == nlayers ||
        throw(DimensionMismatch("ice_water_path must match longwave layers"))
    length(cloud_fraction) == nlayers ||
        throw(DimensionMismatch("cloud_fraction must match longwave layers"))

    exponent = max(Float64(cloud_fraction_exponent), 0.0)
    for k in 1:nlayers
        fraction_scale = clamp(Float64(cloud_fraction[k]), 0.0, 1.0)^exponent
        lwp = max(Float64(liquid_water_path[k]), 0.0)
        iwp = max(Float64(ice_water_path[k]), 0.0)
        for ig in 1:ng
            old_tau = max(Float64(longwave.optical_depth[ig, k]), 0.0)
            old_ssa = clamp(Float64(longwave.single_scattering_albedo[ig, k]), 0.0, 1.0)
            old_g = clamp(Float64(longwave.scattering_asymmetry[ig, k]), -1.0, 1.0)
            old_scattering = old_tau * old_ssa
            old_moment = old_scattering * old_g

            liquid_ext = fraction_scale * liquid_extinction_scale *
                liquid_properties.mass_extinction_coefficient[ig] * lwp
            ice_ext = fraction_scale * ice_extinction_scale *
                ice_properties.mass_extinction_coefficient[ig] * iwp
            liquid_ssa = clamp(liquid_properties.single_scattering_albedo[ig], 0.0, 1.0)
            ice_ssa = clamp(ice_properties.single_scattering_albedo[ig], 0.0, 1.0)
            liquid_g = clamp(liquid_properties.asymmetry_factor[ig], -1.0, 1.0)
            ice_g = clamp(ice_properties.asymmetry_factor[ig], -1.0, 1.0)
            added_scattering = liquid_ssa * liquid_ext + ice_ssa * ice_ext
            added_moment = liquid_ssa * liquid_ext * liquid_g +
                ice_ssa * ice_ext * ice_g
            total_tau = old_tau + liquid_ext + ice_ext
            total_scattering = old_scattering + added_scattering
            longwave.optical_depth[ig, k] = total_tau
            longwave.single_scattering_albedo[ig, k] =
                total_tau <= 0 ? 0.0 : clamp(total_scattering / total_tau, 0.0, 1.0)
            longwave.scattering_asymmetry[ig, k] =
                total_scattering <= 0 ? 0.0 :
                clamp((old_moment + added_moment) / total_scattering, -1.0, 1.0)
        end
    end
    return longwave
end

function layer_cloud_scattering_properties(mapped, re_liquid, re_ice, k; shortwave::Bool)
    key = (
        shortwave = shortwave,
        re_liquid = Float64(re_liquid[k]),
        re_ice = Float64(re_ice[k]),
        mapping_method = mapped.mapping_method,
        delta_eddington_average = mapped.delta_eddington_average,
        thick_averaging = mapped.thick_averaging,
    )
    cache = LAYER_CLOUD_SCATTERING_CACHE[]
    haskey(cache, key) && return cache[key]
    liquid_table = mapped.liquid_table
    ice_table = mapped.ice_table
    mapping = shortwave ? mapped.sw_mapping : mapped.lw_mapping
    properties = (
        liquid = cloud_scattering_gpoint_properties(
            liquid_table, mapping, re_liquid[k];
            mapping_method = mapped.mapping_method,
            delta_eddington_average = mapped.delta_eddington_average,
            thick_averaging = mapped.thick_averaging,
        ),
        ice = cloud_scattering_gpoint_properties(
            ice_table, mapping, re_ice[k];
            mapping_method = mapped.mapping_method,
            delta_eddington_average = mapped.delta_eddington_average,
            thick_averaging = mapped.thick_averaging,
        ),
    )
    cache[key] = properties
    return properties
end

function add_layer_mapped_cloud_scattering!(shortwave,
                                            mapped,
                                            liquid_water_path,
                                            ice_water_path,
                                            cloud_fraction,
                                            re_liquid,
                                            re_ice;
                                            cloud_fraction_exponent = 1,
                                            liquid_extinction_scale = 1,
                                            ice_extinction_scale = 1,
                                            shortwave_scattering_scale = 1,
                                            delta_eddington_scale = true)
    ng, nlayers = size(shortwave.optical_depth)
    exponent = max(Float64(cloud_fraction_exponent), 0.0)
    scattering_scale = max(Float64(shortwave_scattering_scale), 0.0)
    for k in 1:nlayers
        props = layer_cloud_scattering_properties(mapped, re_liquid, re_ice, k;
                                                  shortwave = true)
        fraction_scale = clamp(Float64(cloud_fraction[k]), 0.0, 1.0)^exponent
        lwp = max(Float64(liquid_water_path[k]), 0.0)
        iwp = max(Float64(ice_water_path[k]), 0.0)
        for ig in 1:ng
            liquid_ext = liquid_extinction_scale *
                props.liquid.mass_extinction_coefficient[ig] * lwp
            ice_ext = ice_extinction_scale *
                props.ice.mass_extinction_coefficient[ig] * iwp
            liquid_scat =
                clamp(props.liquid.single_scattering_albedo[ig], 0.0, 1.0) * liquid_ext
            ice_scat =
                clamp(props.ice.single_scattering_albedo[ig], 0.0, 1.0) * ice_ext
            scattering_sum = liquid_scat + ice_scat
            asymmetry = scattering_sum == 0.0 ? 0.0 :
                (clamp(props.liquid.asymmetry_factor[ig], -1.0, 1.0) * liquid_scat +
                 clamp(props.ice.asymmetry_factor[ig], -1.0, 1.0) * ice_scat) /
                scattering_sum
            total_extinction = liquid_ext + ice_ext
            if delta_eddington_scale && scattering_sum > 0.0
                forward_fraction = asymmetry^2
                total_extinction -= scattering_sum * forward_fraction
                scattering_sum *= 1.0 - forward_fraction
                asymmetry /= 1.0 + asymmetry
            end
            if scattering_scale != 1.0
                scattering_sum = min(scattering_sum * scattering_scale,
                                     max(total_extinction, 0.0))
            end
            absorption_tau = fraction_scale * max(total_extinction - scattering_sum, 0.0)
            scattering_tau = fraction_scale * scattering_sum
            shortwave.optical_depth[ig, k] += absorption_tau
            old_scattering = shortwave.rayleigh_optical_depth[ig, k]
            total_scattering = old_scattering + scattering_tau
            shortwave.scattering_asymmetry[ig, k] = total_scattering == 0.0 ?
                0.0 :
                (shortwave.scattering_asymmetry[ig, k] * old_scattering +
                 asymmetry * scattering_tau) / total_scattering
            shortwave.rayleigh_optical_depth[ig, k] = total_scattering
        end
    end
    return shortwave
end

function add_layer_mapped_cloud_longwave!(longwave,
                                          mapped,
                                          liquid_water_path,
                                          ice_water_path,
                                          cloud_fraction,
                                          re_liquid,
                                          re_ice;
                                          cloud_fraction_exponent = 1,
                                          liquid_extinction_scale = 1,
                                          ice_extinction_scale = 1,
                                          scattering = false,
                                          extinction_as_absorption = false)
    ng, nlayers = size(longwave.optical_depth)
    exponent = max(Float64(cloud_fraction_exponent), 0.0)
    for k in 1:nlayers
        props = layer_cloud_scattering_properties(mapped, re_liquid, re_ice, k;
                                                  shortwave = false)
        fraction_scale = clamp(Float64(cloud_fraction[k]), 0.0, 1.0)^exponent
        lwp = max(Float64(liquid_water_path[k]), 0.0)
        iwp = max(Float64(ice_water_path[k]), 0.0)
        for ig in 1:ng
            liquid_ext = fraction_scale * liquid_extinction_scale *
                props.liquid.mass_extinction_coefficient[ig] * lwp
            ice_ext = fraction_scale * ice_extinction_scale *
                props.ice.mass_extinction_coefficient[ig] * iwp
            liquid_ssa = clamp(props.liquid.single_scattering_albedo[ig], 0.0, 1.0)
            ice_ssa = clamp(props.ice.single_scattering_albedo[ig], 0.0, 1.0)
            if scattering
                longwave.single_scattering_albedo === nothing &&
                    throw(ArgumentError("longwave single_scattering_albedo storage is required"))
                longwave.scattering_asymmetry === nothing &&
                    throw(ArgumentError("longwave scattering_asymmetry storage is required"))
                old_tau = max(Float64(longwave.optical_depth[ig, k]), 0.0)
                old_ssa = clamp(Float64(longwave.single_scattering_albedo[ig, k]), 0.0, 1.0)
                old_g = clamp(Float64(longwave.scattering_asymmetry[ig, k]), -1.0, 1.0)
                old_scattering = old_tau * old_ssa
                added_scattering = liquid_ssa * liquid_ext + ice_ssa * ice_ext
                added_moment =
                    liquid_ssa * liquid_ext * clamp(props.liquid.asymmetry_factor[ig], -1.0, 1.0) +
                    ice_ssa * ice_ext * clamp(props.ice.asymmetry_factor[ig], -1.0, 1.0)
                total_tau = old_tau + liquid_ext + ice_ext
                total_scattering = old_scattering + added_scattering
                longwave.optical_depth[ig, k] = total_tau
                longwave.single_scattering_albedo[ig, k] =
                    total_tau == 0.0 ? 0.0 : total_scattering / total_tau
                longwave.scattering_asymmetry[ig, k] =
                    total_scattering == 0.0 ? 0.0 :
                    (old_scattering * old_g + added_moment) / total_scattering
            else
                tau = extinction_as_absorption ?
                    liquid_ext + ice_ext :
                    (1.0 - liquid_ssa) * liquid_ext + (1.0 - ice_ssa) * ice_ext
                longwave.optical_depth[ig, k] += tau
            end
        end
    end
    return longwave
end

function maybe_add_cloud_optics!(longwave, shortwave, dataset, column)
    variables = String.(collect(keys(dataset)))
    all(name in variables for name in ("cloud_fraction", "liquid_water_path", "ice_water_path")) ||
        return false

    cloud_fraction = Array(dataset["cloud_fraction"][:, column])
    liquid_water_path = Array(dataset["liquid_water_path"][:, column])
    ice_water_path = Array(dataset["ice_water_path"][:, column])
    cloud = CloudOpticalProperties(zeros(Float64, length(cloud_fraction)),
                                   zeros(Float64, length(cloud_fraction)))
    if env_bool("RH_CLOUD_SCATTERING_TABLE_OPTICS", false)
        mapped = mapped_cloud_scattering_properties()
        use_layer_radius = env_bool("RH_CLOUD_SCATTERING_LAYER_EFFECTIVE_RADIUS", true) &&
            all(name in variables for name in ("re_liquid", "re_ice"))
        if use_layer_radius
            add_layer_mapped_cloud_scattering!(
                shortwave, mapped, liquid_water_path, ice_water_path,
                cloud_fraction, Array(dataset["re_liquid"][:, column]),
                Array(dataset["re_ice"][:, column]);
                cloud_fraction_exponent = env_float("RH_CLOUD_FRACTION_EXPONENT", 1.0),
                liquid_extinction_scale = env_float("RH_LIQUID_CLOUD_SCATTERING_SW_SCALE", 1.0),
                ice_extinction_scale = env_float("RH_ICE_CLOUD_SCATTERING_SW_SCALE", 1.0),
                shortwave_scattering_scale = env_float("RH_CLOUD_SCATTERING_SW_SSA_SCALE", 1.0),
                delta_eddington_scale = env_bool("RH_CLOUD_SCATTERING_DELTA_EDDINGTON_SCALE", true),
            )
        else
            add_mapped_cloud_scattering!(shortwave, mapped.liquid_sw, mapped.ice_sw,
                                         liquid_water_path, ice_water_path,
                                         cloud_fraction;
                                         cloud_fraction_exponent =
                                             env_float("RH_CLOUD_FRACTION_EXPONENT", 1.0),
                                         liquid_extinction_scale =
                                             env_float("RH_LIQUID_CLOUD_SCATTERING_SW_SCALE", 1.0),
                                         ice_extinction_scale =
                                             env_float("RH_ICE_CLOUD_SCATTERING_SW_SCALE", 1.0),
                                         shortwave_scattering_scale =
                                             env_float("RH_CLOUD_SCATTERING_SW_SSA_SCALE", 1.0),
                                         delta_eddington_scale =
                                             env_bool("RH_CLOUD_SCATTERING_DELTA_EDDINGTON_SCALE", true))
        end
        if env_bool("RH_LW_CLOUD_SCATTERING", longwave.single_scattering_albedo !== nothing)
            if use_layer_radius
                add_layer_mapped_cloud_longwave!(
                    longwave, mapped, liquid_water_path, ice_water_path,
                    cloud_fraction, Array(dataset["re_liquid"][:, column]),
                    Array(dataset["re_ice"][:, column]);
                    cloud_fraction_exponent = env_float("RH_CLOUD_FRACTION_EXPONENT", 1.0),
                    liquid_extinction_scale = env_float("RH_LIQUID_CLOUD_SCATTERING_LW_SCALE", 1.0),
                    ice_extinction_scale = env_float("RH_ICE_CLOUD_SCATTERING_LW_SCALE", 1.0),
                    scattering = true,
                )
            else
                add_mapped_cloud_longwave_scattering!(
                    longwave, mapped.liquid_lw, mapped.ice_lw,
                    liquid_water_path, ice_water_path, cloud_fraction;
                    cloud_fraction_exponent = env_float("RH_CLOUD_FRACTION_EXPONENT", 1.0),
                    liquid_extinction_scale = env_float("RH_LIQUID_CLOUD_SCATTERING_LW_SCALE", 1.0),
                    ice_extinction_scale = env_float("RH_ICE_CLOUD_SCATTERING_LW_SCALE", 1.0),
                )
            end
        else
            if use_layer_radius
                add_layer_mapped_cloud_longwave!(
                    longwave, mapped, liquid_water_path, ice_water_path,
                    cloud_fraction, Array(dataset["re_liquid"][:, column]),
                    Array(dataset["re_ice"][:, column]);
                    cloud_fraction_exponent = env_float("RH_CLOUD_FRACTION_EXPONENT", 1.0),
                    liquid_extinction_scale = env_float("RH_LIQUID_CLOUD_SCATTERING_LW_SCALE", 1.0),
                    ice_extinction_scale = env_float("RH_ICE_CLOUD_SCATTERING_LW_SCALE", 1.0),
                    extinction_as_absorption =
                        env_bool("RH_CLOUD_LW_TABLE_EXTINCTION_AS_ABSORPTION", false),
                )
            else
                add_mapped_cloud_longwave_absorption!(
                    longwave, mapped.liquid_lw, mapped.ice_lw,
                    liquid_water_path, ice_water_path, cloud_fraction;
                    cloud_fraction_exponent = env_float("RH_CLOUD_FRACTION_EXPONENT", 1.0),
                    liquid_extinction_scale = env_float("RH_LIQUID_CLOUD_SCATTERING_LW_SCALE", 1.0),
                    ice_extinction_scale = env_float("RH_ICE_CLOUD_SCATTERING_LW_SCALE", 1.0),
                    extinction_as_absorption =
                        env_bool("RH_CLOUD_LW_TABLE_EXTINCTION_AS_ABSORPTION", false),
                )
            end
        end
        return true
    end
    if env_bool("RH_CLOUD_EFFECTIVE_RADIUS_OPTICS", false) &&
       all(name in variables for name in ("re_liquid", "re_ice"))
        re_liquid = Array(dataset["re_liquid"][:, column])
        re_ice = Array(dataset["re_ice"][:, column])
        fraction_exponent = env_float("RH_CLOUD_FRACTION_EXPONENT", 1.0)
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
        liquid_sw_scale = env_float("RH_LIQUID_CLOUD_EFFECTIVE_RADIUS_SW_SCALE", 1.0)
        ice_sw_scale = env_float("RH_ICE_CLOUD_EFFECTIVE_RADIUS_SW_SCALE", 1.0)
        liquid_lw_absorption =
            env_float("RH_LIQUID_CLOUD_LW_MASS_ABSORPTION",
                      env_float("RH_CLOUD_LW_MASS_ABSORPTION", 0.0))
        ice_lw_absorption =
            env_float("RH_ICE_CLOUD_LW_MASS_ABSORPTION",
                      env_float("RH_CLOUD_LW_MASS_ABSORPTION", 0.0))

        for k in eachindex(cloud_fraction)
            scale = clamp(cloud_fraction[k], 0.0, 1.0)^max(fraction_exponent, 0.0)
            liquid_extinction = liquid_sw_scale * 1.5 * liquid_water_path[k] /
                (LIQUID_WATER_DENSITY * max(re_liquid[k], eps(Float64)))
            ice_extinction = ice_sw_scale * 1.5 * ice_water_path[k] /
                (ICE_DENSITY * max(re_ice[k], eps(Float64)))
            liquid_scat = liquid_ssa * liquid_extinction
            ice_scat = ice_ssa * ice_extinction
            total_scat = liquid_scat + ice_scat
            cloud.longwave_optical_depth[k] = scale *
                (liquid_lw_absorption * liquid_water_path[k] +
                 ice_lw_absorption * ice_water_path[k])
            cloud.shortwave_optical_depth[k] = scale *
                ((1.0 - liquid_ssa) * liquid_extinction +
                 (1.0 - ice_ssa) * ice_extinction)
            cloud.shortwave_scattering_optical_depth[k] = scale * total_scat
            cloud.shortwave_scattering_asymmetry[k] = total_scat == 0 ?
                0.0 :
                (liquid_asymmetry * liquid_scat + ice_asymmetry * ice_scat) /
                total_scat
        end
        add_cloud_optical_depths!(longwave, shortwave, cloud)
        return true
    end

    cloud_model = LayerLiquidIceCloudOpticsModel(
        liquid_water_path = liquid_water_path,
        ice_water_path = ice_water_path,
        cloud_fraction = cloud_fraction,
        liquid_longwave_mass_absorption =
            env_float("RH_LIQUID_CLOUD_LW_MASS_ABSORPTION",
                      env_float("RH_CLOUD_LW_MASS_ABSORPTION", 0.0)),
        ice_longwave_mass_absorption =
            env_float("RH_ICE_CLOUD_LW_MASS_ABSORPTION",
                      env_float("RH_CLOUD_LW_MASS_ABSORPTION", 0.0)),
        liquid_shortwave_mass_extinction =
            env_float("RH_LIQUID_CLOUD_SW_MASS_EXTINCTION",
                      env_float("RH_CLOUD_SW_MASS_EXTINCTION", 0.0)),
        ice_shortwave_mass_extinction =
            env_float("RH_ICE_CLOUD_SW_MASS_EXTINCTION",
                      env_float("RH_CLOUD_SW_MASS_EXTINCTION", 0.0)),
        liquid_shortwave_single_scattering_albedo =
            env_float("RH_LIQUID_CLOUD_SW_SINGLE_SCATTERING_ALBEDO",
                      env_float("RH_CLOUD_SW_SINGLE_SCATTERING_ALBEDO", 0.0)),
        ice_shortwave_single_scattering_albedo =
            env_float("RH_ICE_CLOUD_SW_SINGLE_SCATTERING_ALBEDO",
                      env_float("RH_CLOUD_SW_SINGLE_SCATTERING_ALBEDO", 0.0)),
        liquid_shortwave_scattering_asymmetry =
            env_float("RH_LIQUID_CLOUD_SW_SCATTERING_ASYMMETRY",
                      env_float("RH_CLOUD_SW_SCATTERING_ASYMMETRY", 0.0)),
        ice_shortwave_scattering_asymmetry =
            env_float("RH_ICE_CLOUD_SW_SCATTERING_ASYMMETRY",
                      env_float("RH_CLOUD_SW_SCATTERING_ASYMMETRY", 0.0)),
        cloud_fraction_exponent =
            env_float("RH_CLOUD_FRACTION_EXPONENT", 1.0),
    )
    Lightflux.cloud_optical_properties!(cloud, cloud_model, (;))
    add_cloud_optical_depths!(longwave, shortwave, cloud)
    return true
end

function build_cloudy_region_shortwave(clear_shortwave, dataset, column)
    variables = String.(collect(keys(dataset)))
    all(name in variables for name in ("cloud_fraction", "liquid_water_path", "ice_water_path")) ||
        return clear_shortwave
    env_bool("RH_CLOUD_SCATTERING_TABLE_OPTICS", false) || return clear_shortwave

    cloudy = deepcopy(clear_shortwave)
    cloud_fraction = Array(dataset["cloud_fraction"][:, column])
    liquid_water_path = Array(dataset["liquid_water_path"][:, column])
    ice_water_path = Array(dataset["ice_water_path"][:, column])
    cloudy_liquid_water_path = similar(liquid_water_path)
    cloudy_ice_water_path = similar(ice_water_path)
    full_cloud_fraction = ones(Float64, length(cloud_fraction))
    for k in eachindex(cloud_fraction)
        fraction = clamp(cloud_fraction[k], 0.0, 1.0)
        scale = fraction > 0 ? inv(fraction) : 0.0
        cloudy_liquid_water_path[k] = max(liquid_water_path[k], 0.0) * scale
        cloudy_ice_water_path[k] = max(ice_water_path[k], 0.0) * scale
    end

    mapped = mapped_cloud_scattering_properties()
    if env_bool("RH_CLOUD_SCATTERING_LAYER_EFFECTIVE_RADIUS", true) &&
       all(name in variables for name in ("re_liquid", "re_ice"))
        add_layer_mapped_cloud_scattering!(
            cloudy, mapped, cloudy_liquid_water_path, cloudy_ice_water_path,
            full_cloud_fraction, Array(dataset["re_liquid"][:, column]),
            Array(dataset["re_ice"][:, column]);
            cloud_fraction_exponent = 1.0,
            liquid_extinction_scale = env_float("RH_LIQUID_CLOUD_SCATTERING_SW_SCALE", 1.0),
            ice_extinction_scale = env_float("RH_ICE_CLOUD_SCATTERING_SW_SCALE", 1.0),
            shortwave_scattering_scale = env_float("RH_CLOUD_SCATTERING_SW_SSA_SCALE", 1.0),
            delta_eddington_scale = env_bool("RH_CLOUD_SCATTERING_DELTA_EDDINGTON_SCALE", true),
        )
    else
        add_mapped_cloud_scattering!(cloudy, mapped.liquid_sw, mapped.ice_sw,
                                     cloudy_liquid_water_path, cloudy_ice_water_path,
                                     full_cloud_fraction;
                                     cloud_fraction_exponent = 1.0,
                                     liquid_extinction_scale =
                                         env_float("RH_LIQUID_CLOUD_SCATTERING_SW_SCALE", 1.0),
                                     ice_extinction_scale =
                                         env_float("RH_ICE_CLOUD_SCATTERING_SW_SCALE", 1.0),
                                     shortwave_scattering_scale =
                                         env_float("RH_CLOUD_SCATTERING_SW_SSA_SCALE", 1.0),
                                     delta_eddington_scale =
                                         env_bool("RH_CLOUD_SCATTERING_DELTA_EDDINGTON_SCALE", true))
    end
    return cloudy
end

function build_cloudy_region_longwave(clear_longwave, dataset, column)
    variables = String.(collect(keys(dataset)))
    all(name in variables for name in ("cloud_fraction", "liquid_water_path", "ice_water_path")) ||
        return clear_longwave
    env_bool("RH_CLOUD_SCATTERING_TABLE_OPTICS", false) || return clear_longwave

    cloudy = deepcopy(clear_longwave)
    cloud_fraction = Array(dataset["cloud_fraction"][:, column])
    liquid_water_path = Array(dataset["liquid_water_path"][:, column])
    ice_water_path = Array(dataset["ice_water_path"][:, column])
    cloudy_liquid_water_path = similar(liquid_water_path)
    cloudy_ice_water_path = similar(ice_water_path)
    full_cloud_fraction = ones(Float64, length(cloud_fraction))
    for k in eachindex(cloud_fraction)
        fraction = clamp(cloud_fraction[k], 0.0, 1.0)
        scale = fraction > 0 ? inv(fraction) : 0.0
        cloudy_liquid_water_path[k] = max(liquid_water_path[k], 0.0) * scale
        cloudy_ice_water_path[k] = max(ice_water_path[k], 0.0) * scale
    end

    mapped = mapped_cloud_scattering_properties()
    use_layer_radius = env_bool("RH_CLOUD_SCATTERING_LAYER_EFFECTIVE_RADIUS", true) &&
        all(name in variables for name in ("re_liquid", "re_ice"))
    if env_bool("RH_LW_CLOUD_SCATTERING", false)
        if use_layer_radius
            add_layer_mapped_cloud_longwave!(
                cloudy, mapped, cloudy_liquid_water_path, cloudy_ice_water_path,
                full_cloud_fraction, Array(dataset["re_liquid"][:, column]),
                Array(dataset["re_ice"][:, column]);
                cloud_fraction_exponent = 1.0,
                liquid_extinction_scale = env_float("RH_LIQUID_CLOUD_SCATTERING_LW_SCALE", 1.0),
                ice_extinction_scale = env_float("RH_ICE_CLOUD_SCATTERING_LW_SCALE", 1.0),
                scattering = true,
            )
        else
            add_mapped_cloud_longwave_scattering!(
                cloudy, mapped.liquid_lw, mapped.ice_lw,
                cloudy_liquid_water_path, cloudy_ice_water_path, full_cloud_fraction;
                cloud_fraction_exponent = 1.0,
                liquid_extinction_scale = env_float("RH_LIQUID_CLOUD_SCATTERING_LW_SCALE", 1.0),
                ice_extinction_scale = env_float("RH_ICE_CLOUD_SCATTERING_LW_SCALE", 1.0),
            )
        end
    else
        if use_layer_radius
            add_layer_mapped_cloud_longwave!(
                cloudy, mapped, cloudy_liquid_water_path, cloudy_ice_water_path,
                full_cloud_fraction, Array(dataset["re_liquid"][:, column]),
                Array(dataset["re_ice"][:, column]);
                cloud_fraction_exponent = 1.0,
                liquid_extinction_scale = env_float("RH_LIQUID_CLOUD_SCATTERING_LW_SCALE", 1.0),
                ice_extinction_scale = env_float("RH_ICE_CLOUD_SCATTERING_LW_SCALE", 1.0),
                extinction_as_absorption =
                    env_bool("RH_CLOUD_LW_TABLE_EXTINCTION_AS_ABSORPTION", false),
            )
        else
            add_mapped_cloud_longwave_absorption!(
                cloudy, mapped.liquid_lw, mapped.ice_lw,
                cloudy_liquid_water_path, cloudy_ice_water_path, full_cloud_fraction;
                cloud_fraction_exponent = 1.0,
                liquid_extinction_scale = env_float("RH_LIQUID_CLOUD_SCATTERING_LW_SCALE", 1.0),
                ice_extinction_scale = env_float("RH_ICE_CLOUD_SCATTERING_LW_SCALE", 1.0),
                extinction_as_absorption =
                    env_bool("RH_CLOUD_LW_TABLE_EXTINCTION_AS_ABSORPTION", false),
            )
        end
    end
    return cloudy
end

function aerosol_path_from_mmr(dataset, pressure_interfaces, column)
    mmr = Array(dataset["aerosol_mmr"][:, :, column])
    nlayers = size(mmr, 1)
    path = zeros(Float64, nlayers)
    for k in 1:nlayers
        Δp = pressure_interfaces[k + 1, column] - pressure_interfaces[k, column]
        path[k] = sum(mmr[k, :]) * Δp / GRAVITY
    end
    return path
end

function _nearest_index(values, target)
    best = firstindex(values)
    best_distance = abs(Float64(values[best]) - Float64(target))
    for index in eachindex(values)
        distance = abs(Float64(values[index]) - Float64(target))
        if distance < best_distance
            best = index
            best_distance = distance
        end
    end
    return best
end

function _ecrad_rh_index(rh_lower, rh)
    if rh > rh_lower[end]
        return lastindex(rh_lower)
    end
    index = firstindex(rh_lower)
    while index < lastindex(rh_lower) && rh > rh_lower[index + 1]
        index += 1
    end
    return index
end

function _saturation_vapor_pressure_liquid(temperature)
    celsius = Float64(temperature) - 273.15
    return 611.2 * exp(17.67 * celsius / (celsius + 243.5))
end

function _saturation_specific_humidity_liquid(temperature, pressure)
    vapor_pressure = min(_saturation_vapor_pressure_liquid(temperature),
                         0.99 * Float64(pressure))
    return WATER_VAPOR_EPSILON * vapor_pressure /
        (Float64(pressure) - (1.0 - WATER_VAPOR_EPSILON) * vapor_pressure)
end

function _layer_relative_humidity(dataset, k, column)
    variables = String.(collect(keys(dataset)))
    if all(name in variables for name in ("h2o", "h2o_sat_liq"))
        q = max(Float64(dataset["h2o"][k, column]), 0.0)
        qsat = max(Float64(dataset["h2o_sat_liq"][k, column]), 0.0)
        return qsat > 0 ? clamp(q / qsat, 0.0, 1.0) :
            env_float("RH_IFS_AEROSOL_RELATIVE_HUMIDITY", 0.8)
    end
    all(name in variables for name in ("h2o", "temperature_layer", "pressure_layer")) ||
        return env_float("RH_IFS_AEROSOL_RELATIVE_HUMIDITY", 0.8)
    q = max(Float64(dataset["h2o"][k, column]), 0.0)
    qsat = _saturation_specific_humidity_liquid(dataset["temperature_layer"][k, column],
                                                dataset["pressure_layer"][k, column])
    return qsat > 0 ? clamp(q / qsat, 0.0, 1.0) :
        env_float("RH_IFS_AEROSOL_RELATIVE_HUMIDITY", 0.8)
end

function _mapped_ifs_aerosol_properties()
    cached = IFS_AEROSOL_CACHE[]
    cached === nothing || return cached
    nc = require_ncdatasets()
    path = reference_path(get(ENV, "RH_IFS_AEROSOL_OPTICS_PATH",
        "validation/external/ecrad/data/aerosol_ifs_49R1_20230119.nc"))
    sw_mapping_path = reference_path(get(ENV, "RH_CLOUD_SW_MAPPING_PATH",
        "validation/external/ecrad/data/ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc"))
    lw_mapping_path = reference_path(get(ENV, "RH_CLOUD_LW_MAPPING_PATH",
        "validation/external/ecrad/data/ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc"))
    sw_mapping = read_ecckd_spectral_mapping(sw_mapping_path)
    lw_mapping = read_ecckd_spectral_mapping(lw_mapping_path)
    aerosol_mapping = CloudScatteringTable(
        wavenumber = [0.0, 1.0],
        effective_radius = [1.0],
        mass_extinction_coefficient = zeros(2, 1),
        single_scattering_albedo = zeros(2, 1),
        asymmetry_factor = zeros(2, 1),
    )
    nc.NCDataset(path) do dataset
        wavenumber = Float64.(Array(dataset["wavenumber"]))
        table = CloudScatteringTable(
            wavenumber = wavenumber,
            effective_radius = [1.0],
            mass_extinction_coefficient = zeros(length(wavenumber), 1),
            single_scattering_albedo = zeros(length(wavenumber), 1),
            asymmetry_factor = zeros(length(wavenumber), 1),
        )
        sw_weights = Lightflux._ecrad_cloud_mapping_matrix(table, sw_mapping)
        lw_weights = Lightflux._ecrad_cloud_mapping_matrix(table, lw_mapping)
        rh_lower = Float64.(Array(dataset["relative_humidity1"]))
        ng_sw = size(sw_weights, 1)
        ng_lw = size(lw_weights, 1)
        ntypes = length(IFS_AEROSOL_TYPE_MAP)
        nrh = length(rh_lower)
        sw_ext = zeros(Float64, ng_sw, nrh, ntypes)
        sw_ssa = zeros(Float64, ng_sw, nrh, ntypes)
        sw_asymmetry = zeros(Float64, ng_sw, nrh, ntypes)
        lw_absorption = zeros(Float64, nrh, ng_lw, ntypes)
        for (itype, mapped_type) in enumerate(IFS_AEROSOL_TYPE_MAP)
            rh_indices = mapped_type < 0 ? (1:nrh) : (1:1)
            for irh in rh_indices
                if mapped_type < 0
                    index = abs(mapped_type)
                    raw_ext_sw = Float64.(Array(dataset["mass_ext_hydrophilic"][:, irh, index]))
                    raw_ssa_sw = clamp.(Float64.(Array(dataset["ssa_hydrophilic"][:, irh, index])), 0.0, 1.0)
                    raw_asym_sw = clamp.(Float64.(Array(dataset["asymmetry_hydrophilic"][:, irh, index])), -1.0, 1.0)
                    raw_ext_lw = raw_ext_sw
                    raw_ssa_lw = raw_ssa_sw
                else
                    index = mapped_type
                    raw_ext_sw = Float64.(Array(dataset["mass_ext_hydrophobic"][:, index]))
                    raw_ssa_sw = clamp.(Float64.(Array(dataset["ssa_hydrophobic"][:, index])), 0.0, 1.0)
                    raw_asym_sw = clamp.(Float64.(Array(dataset["asymmetry_hydrophobic"][:, index])), -1.0, 1.0)
                    raw_ext_lw = raw_ext_sw
                    raw_ssa_lw = raw_ssa_sw
                end
                raw_scat_sw = raw_ext_sw .* raw_ssa_sw
                for ig in 1:ng_sw
                    ext = sum(view(sw_weights, ig, :) .* raw_ext_sw)
                    scat = sum(view(sw_weights, ig, :) .* raw_scat_sw)
                    sw_ext[ig, irh, itype] = ext
                    sw_ssa[ig, irh, itype] = ext > 0 ? clamp(scat / ext, 0.0, 1.0) : 0.0
                    sw_asymmetry[ig, irh, itype] = scat > 0 ?
                        clamp(sum(view(sw_weights, ig, :) .* raw_scat_sw .* raw_asym_sw) / scat,
                              -1.0, 1.0) : 0.0
                end
                for ig in 1:ng_lw
                    lw_absorption[irh, ig, itype] =
                        sum(view(lw_weights, ig, :) .* raw_ext_lw .* (1.0 .- raw_ssa_lw))
                end
            end
            if mapped_type > 0 && nrh > 1
                for irh in 2:nrh
                    sw_ext[:, irh, itype] .= sw_ext[:, 1, itype]
                    sw_ssa[:, irh, itype] .= sw_ssa[:, 1, itype]
                    sw_asymmetry[:, irh, itype] .= sw_asymmetry[:, 1, itype]
                    lw_absorption[irh, :, itype] .= lw_absorption[1, :, itype]
                end
            end
        end
        cached = (
            sw_extinction = sw_ext,
            sw_single_scattering_albedo = sw_ssa,
            sw_asymmetry = sw_asymmetry,
            lw_absorption = lw_absorption,
            relative_humidity_lower_bounds = rh_lower,
            path = path,
        )
    end
    IFS_AEROSOL_CACHE[] = cached
    return cached
end

function add_ifs_aerosol_optics!(longwave, shortwave, dataset, pressure_interfaces, column)
    variables = String.(collect(keys(dataset)))
    "aerosol_mmr" in variables || return false
    props = _mapped_ifs_aerosol_properties()
    mmr = Array(dataset["aerosol_mmr"][:, :, column])
    nlayers = size(mmr, 1)
    ng_lw = size(longwave.optical_depth, 1)
    ng_sw = size(shortwave.optical_depth, 1)
    if size(props.lw_absorption, 2) != ng_lw || size(props.sw_extinction, 1) != ng_sw
        return false
    end
    ntypes = min(size(mmr, 2), size(props.sw_extinction, 3))
    use_layer_rh = env_bool("RH_IFS_AEROSOL_LAYER_RELATIVE_HUMIDITY", true)
    delta_scale_sw = env_bool("RH_AEROSOL_DELTA_EDDINGTON_SCALE", true)
    for k in 1:nlayers
        Δp = pressure_interfaces[k + 1, column] - pressure_interfaces[k, column]
        paths = max.(Float64.(view(mmr, k, 1:ntypes)), 0.0) .* Δp ./ GRAVITY
        rh = use_layer_rh ? _layer_relative_humidity(dataset, k, column) :
            env_float("RH_IFS_AEROSOL_RELATIVE_HUMIDITY", 0.8)
        irh = _ecrad_rh_index(props.relative_humidity_lower_bounds, rh)
        for ig in 1:ng_lw
            longwave.optical_depth[ig, k] +=
                sum(view(props.lw_absorption, irh, ig, 1:ntypes) .* paths)
        end
        for ig in 1:ng_sw
            ext_by_type = view(props.sw_extinction, ig, irh, 1:ntypes) .* paths
            ext = sum(ext_by_type)
            scat_by_type =
                ext_by_type .* view(props.sw_single_scattering_albedo, ig, irh, 1:ntypes)
            scat = sum(scat_by_type)
            aerosol_asym = scat > 0 ?
                sum(scat_by_type .* view(props.sw_asymmetry, ig, irh, 1:ntypes)) / scat :
                0.0
            if delta_scale_sw && scat > 0
                forward_fraction = aerosol_asym^2
                ext -= scat * forward_fraction
                scat *= 1.0 - forward_fraction
                aerosol_asym /= 1.0 + aerosol_asym
            end
            shortwave.optical_depth[ig, k] += max(ext - scat, 0.0)
            shortwave.rayleigh_optical_depth[ig, k] += scat
            if scat > 0
                old_scat = shortwave.rayleigh_optical_depth[ig, k] - scat
                old_asym = shortwave.scattering_asymmetry[ig, k]
                shortwave.scattering_asymmetry[ig, k] =
                    (old_scat * old_asym + scat * aerosol_asym) / (old_scat + scat)
            end
        end
    end
    return true
end

function maybe_add_aerosol_optics!(longwave, shortwave, dataset, pressure_interfaces, column)
    variables = String.(collect(keys(dataset)))
    env_bool("RH_AEROSOL_OPTICS", false) || return false
    "aerosol_mmr" in variables || return false
    env_bool("RH_IFS_AEROSOL_TABLE_OPTICS", false) &&
        return add_ifs_aerosol_optics!(longwave, shortwave, dataset,
                                       pressure_interfaces, column)

    aerosol_path = aerosol_path_from_mmr(dataset, pressure_interfaces, column)
    aerosol = AerosolOpticalProperties(zeros(Float64, length(aerosol_path)),
                                       zeros(Float64, length(aerosol_path)))
    model = LayerAerosolOpticsModel(
        aerosol_path = aerosol_path,
        longwave_mass_absorption = env_float("RH_AEROSOL_LW_MASS_ABSORPTION", 0.0),
        shortwave_mass_extinction = env_float("RH_AEROSOL_SW_MASS_EXTINCTION", 0.0),
        shortwave_single_scattering_albedo =
            env_float("RH_AEROSOL_SW_SINGLE_SCATTERING_ALBEDO", 0.0),
        shortwave_scattering_asymmetry =
            env_float("RH_AEROSOL_SW_SCATTERING_ASYMMETRY", 0.0),
    )
    aerosol_optical_properties!(aerosol, model, (;))
    add_aerosol_optical_depths!(longwave, shortwave, aerosol)
    return true
end

function longwave_surface_boundary(gas_optics, surface_temperature, target_broadband_flux)
    if gas_optics isa EcCKDTabulatedGasOpticsModel &&
       gas_optics.longwave_source_table !== nothing
        spectral = [
            Lightflux._longwave_source(gas_optics, ig, surface_temperature)
            for ig in axes(gas_optics.longwave_absorption, 1)
        ]
        broadband = sum(gas_optics.longwave_weights .* spectral)
        if isfinite(broadband) && broadband > 0
            return spectral .* (target_broadband_flux / broadband)
        end
    end
    return target_broadband_flux
end

function run_candidate_for_file!(path)
    nc = require_ncdatasets()
    gas_optics = candidate_gas_optics(Float64)
    nc.NCDataset(reference_path(path), "a") do dataset
        pressure_layers = Array(dataset["pressure_layer"])
        pressure_interfaces = Array(dataset["pressure_interface"])
        temperature_layers = Array(dataset["temperature_layer"])
        temperature_interfaces = Array(dataset["temperature_interface"])
        gas_amounts = gas_column_amounts(dataset, pressure_interfaces)
        surface_temperature = Array(dataset["surface_temperature"])
        surface_albedo = Array(dataset["surface_albedo"])
        variables = String.(collect(keys(dataset)))
        surface_albedo_spectral = "surface_albedo_spectral" in variables ?
            Array(dataset["surface_albedo_spectral"]) : nothing
        surface_albedo_direct_spectral =
            "surface_albedo_direct_gpoint" in variables ?
            Array(dataset["surface_albedo_direct_gpoint"]) :
            "surface_albedo_direct_spectral" in variables ?
            Array(dataset["surface_albedo_direct_spectral"]) : nothing
        surface_longwave_up_spectral = "surface_longwave_up_spectral" in variables ?
            Array(dataset["surface_longwave_up_spectral"]) : nothing
        toa_shortwave_down_spectral = "toa_shortwave_down_spectral" in variables ?
            Array(dataset["toa_shortwave_down_spectral"]) : nothing
        if surface_albedo_spectral !== nothing &&
           size(surface_albedo_spectral, 1) != size(gas_optics.shortwave_absorption, 1)
            surface_albedo_spectral = nothing
        end
        if surface_albedo_direct_spectral !== nothing &&
           size(surface_albedo_direct_spectral, 1) != size(gas_optics.shortwave_absorption, 1)
            surface_albedo_direct_spectral = nothing
        end
        if toa_shortwave_down_spectral !== nothing &&
           size(toa_shortwave_down_spectral, 1) != size(gas_optics.shortwave_absorption, 1)
            toa_shortwave_down_spectral = nothing
        end
        if surface_longwave_up_spectral !== nothing &&
           size(surface_longwave_up_spectral, 1) != size(gas_optics.longwave_absorption, 1)
            surface_longwave_up_spectral = nothing
        end
        reference_lw_up = Array(dataset["lw_up"])
        reference_sw_up = Array(dataset["sw_up"])
        reference_sw_down = Array(dataset["sw_down"])
        toa_shortwave_down = Array(dataset["sw_down"])[1, :]
        solar_irradiance = "solar_irradiance" in variables ?
            Array(dataset["solar_irradiance"]) :
            fill(DEFAULT_SOLAR_IRRADIANCE, length(toa_shortwave_down))
        cos_zenith = "cos_solar_zenith_angle" in variables ?
            Array(dataset["cos_solar_zenith_angle"]) :
            clamp.(toa_shortwave_down ./ solar_irradiance, 0.0, 1.0)

        nlayers, ncolumns = size(pressure_layers)
        lw_up = zeros(Float64, nlayers + 1, ncolumns)
        lw_down = zeros(Float64, nlayers + 1, ncolumns)
        sw_up = zeros(Float64, nlayers + 1, ncolumns)
        sw_down = zeros(Float64, nlayers + 1, ncolumns)
        heating = zeros(Float64, nlayers, ncolumns)
        clear_lw_up = zeros(Float64, nlayers + 1, ncolumns)
        clear_lw_down = zeros(Float64, nlayers + 1, ncolumns)
        clear_sw_up = zeros(Float64, nlayers + 1, ncolumns)
        clear_sw_down = zeros(Float64, nlayers + 1, ncolumns)
        clear_heating = zeros(Float64, nlayers, ncolumns)

        use_interface_sources = env_bool("RH_LW_INTERFACE_SOURCES",
                                         candidate_mode() == "official_ecckd")
        use_lw_scattering = env_bool(
            "RH_LW_CLOUD_SCATTERING",
            false,
        )
        ng_lw = size(gas_optics.longwave_absorption, 1)
        longwave = if use_interface_sources
            LongwaveOpticalProperties(
                zeros(Float64, ng_lw, nlayers),
                zeros(Float64, ng_lw, nlayers),
                source_top = zeros(Float64, ng_lw, nlayers),
                source_bottom = zeros(Float64, ng_lw, nlayers),
                single_scattering_albedo =
                    use_lw_scattering ? zeros(Float64, ng_lw, nlayers) : nothing,
                scattering_asymmetry =
                    use_lw_scattering ? zeros(Float64, ng_lw, nlayers) : nothing,
            )
        else
            LongwaveOpticalProperties(
                zeros(Float64, ng_lw, nlayers),
                zeros(Float64, ng_lw, nlayers),
                single_scattering_albedo =
                    use_lw_scattering ? zeros(Float64, ng_lw, nlayers) : nothing,
                scattering_asymmetry =
                    use_lw_scattering ? zeros(Float64, ng_lw, nlayers) : nothing,
            )
        end
        shortwave = ShortwaveOpticalProperties(
            zeros(Float64, size(gas_optics.shortwave_absorption, 1), nlayers),
        )
        fluxes = RadiativeFluxes(
            longwave_up = zeros(Float64, nlayers + 1),
            longwave_down = zeros(Float64, nlayers + 1),
            shortwave_up = zeros(Float64, nlayers + 1),
            shortwave_down = zeros(Float64, nlayers + 1),
        )
        column_heating = zeros(Float64, nlayers)
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
            if longwave.single_scattering_albedo !== nothing
                fill!(longwave.single_scattering_albedo, 0.0)
                fill!(longwave.scattering_asymmetry, 0.0)
            end
            optical_properties!(longwave, shortwave, gas_optics, atmosphere)
            if toa_shortwave_down_spectral !== nothing
                incoming = max.(toa_shortwave_down_spectral[:, j], 0.0)
                incoming_sum = sum(incoming)
                if incoming_sum > 0
                    shortwave.weights .= incoming ./ incoming_sum
                end
            end

            radiative_fluxes!(
                fluxes,
                CloudlessLongwave(),
                longwave,
                atmosphere,
                LongwaveBoundaryConditions(
                    surface_longwave_up = surface_longwave_up_spectral === nothing ?
                        longwave_surface_boundary(
                            gas_optics,
                            surface_temperature[j],
                            reference_lw_up[end, j],
                        ) :
                        surface_longwave_up_spectral[:, j] ./ gas_optics.longwave_weights,
                ),
            )
            effective_surface_albedo = reference_sw_down[end, j] <= 0 ?
                surface_albedo[j] :
                clamp(reference_sw_up[end, j] / reference_sw_down[end, j], 0.0, 1.0)
            shortwave_surface_albedo = surface_albedo_spectral === nothing ?
                effective_surface_albedo :
                surface_albedo_spectral[:, j]
            shortwave_surface_albedo_direct = surface_albedo_direct_spectral === nothing ?
                shortwave_surface_albedo :
                surface_albedo_direct_spectral[:, j]
            radiative_fluxes!(
                fluxes,
                CloudlessShortwave(
                    rayleigh_backscatter_fraction =
                        env_float("RH_SW_RAYLEIGH_BACKSCATTER_FRACTION", 0.5),
                ),
                shortwave,
                atmosphere,
                ShortwaveBoundaryConditions(
                    toa_shortwave_down = max(toa_shortwave_down[j], 0.0),
                    surface_albedo = shortwave_surface_albedo,
                    surface_albedo_direct = shortwave_surface_albedo_direct,
                ),
            )
            heating_rates!(
                column_heating,
                fluxes,
                atmosphere;
                gravity = 9.80665,
                heat_capacity = 1004.0,
            )
            clear_lw_up[:, j] = fluxes.longwave_up
            clear_lw_down[:, j] = fluxes.longwave_down
            clear_sw_up[:, j] = fluxes.shortwave_up
            clear_sw_down[:, j] = fluxes.shortwave_down
            clear_heating[:, j] = 86400.0 .* column_heating

            clear_longwave_for_overlap = deepcopy(longwave)
            clear_shortwave_for_overlap = deepcopy(shortwave)
            maybe_add_aerosol_optics!(clear_longwave_for_overlap,
                                      clear_shortwave_for_overlap,
                                      dataset, pressure_interfaces, j)
            fill!(fluxes.longwave_up, 0.0)
            fill!(fluxes.longwave_down, 0.0)
            fill!(fluxes.shortwave_up, 0.0)
            fill!(fluxes.shortwave_down, 0.0)
            radiative_fluxes!(
                fluxes,
                CloudlessLongwave(),
                clear_longwave_for_overlap,
                atmosphere,
                LongwaveBoundaryConditions(
                    surface_longwave_up = surface_longwave_up_spectral === nothing ?
                        longwave_surface_boundary(
                            gas_optics,
                            surface_temperature[j],
                            reference_lw_up[end, j],
                        ) :
                        surface_longwave_up_spectral[:, j] ./ gas_optics.longwave_weights,
                ),
            )
            radiative_fluxes!(
                fluxes,
                CloudlessShortwave(
                    rayleigh_backscatter_fraction =
                        env_float("RH_SW_RAYLEIGH_BACKSCATTER_FRACTION", 0.5),
                ),
                clear_shortwave_for_overlap,
                atmosphere,
                ShortwaveBoundaryConditions(
                    toa_shortwave_down = max(toa_shortwave_down[j], 0.0),
                    surface_albedo = shortwave_surface_albedo,
                    surface_albedo_direct = shortwave_surface_albedo_direct,
                ),
            )
            heating_rates!(
                column_heating,
                fluxes,
                atmosphere;
                gravity = 9.80665,
                heat_capacity = 1004.0,
            )
            clear_lw_up[:, j] = fluxes.longwave_up
            clear_lw_down[:, j] = fluxes.longwave_down
            clear_sw_up[:, j] = fluxes.shortwave_up
            clear_sw_down[:, j] = fluxes.shortwave_down
            clear_heating[:, j] = 86400.0 .* column_heating

            maybe_add_cloud_optics!(longwave, shortwave, dataset, j)
            maybe_add_aerosol_optics!(longwave, shortwave, dataset, pressure_interfaces, j)

            longwave_boundary = LongwaveBoundaryConditions(
                surface_longwave_up = surface_longwave_up_spectral === nothing ?
                    longwave_surface_boundary(
                        gas_optics,
                        surface_temperature[j],
                        reference_lw_up[end, j],
                    ) :
                    surface_longwave_up_spectral[:, j] ./ gas_optics.longwave_weights,
            )
            if env_bool("RH_CLOUD_OVERLAP_LONGWAVE", false) &&
               "cloud_fraction" in String.(collect(keys(dataset)))
                cloudy_longwave_for_overlap =
                    env_bool("RH_CLOUD_OVERLAP_USE_CLOUDY_REGION_OPTICS", false) ?
                    build_cloudy_region_longwave(clear_longwave_for_overlap, dataset, j) :
                    longwave
                longwave_overlap_optics = LongwaveCloudOverlapOpticalProperties(
                    clear_longwave_for_overlap,
                    cloudy_longwave_for_overlap,
                    Array(dataset["cloud_fraction"][:, j]);
                    overlap_parameter = "overlap_param" in String.(collect(keys(dataset))) ?
                        Array(dataset["overlap_param"][:, j]) : nothing,
                    fractional_std = "fractional_std" in String.(collect(keys(dataset))) ?
                        Array(dataset["fractional_std"][:, j]) : nothing,
                )
                radiative_fluxes!(
                    fluxes,
                    CloudOverlapLongwave(
                        overlap = Symbol(get(ENV, "RH_CLOUD_OVERLAP_LONGWAVE_RULE", "adding")),
                        cloud_fraction_exponent =
                            env_float("RH_CLOUD_FRACTION_EXPONENT", 1.0),
                        inhomogeneity_overlap_exponent =
                            env_float("RH_LW_CLOUD_INHOM_OVERLAP_EXPONENT",
                                      env_float("RH_CLOUD_INHOM_OVERLAP_EXPONENT", 2.0)),
                    ),
                    longwave_overlap_optics,
                    atmosphere,
                    longwave_boundary,
                )
            else
                radiative_fluxes!(
                    fluxes,
                    CloudlessLongwave(),
                    longwave,
                    atmosphere,
                    longwave_boundary,
                )
            end
            shortwave_boundary = ShortwaveBoundaryConditions(
                toa_shortwave_down = max(toa_shortwave_down[j], 0.0),
                surface_albedo = shortwave_surface_albedo,
                surface_albedo_direct = shortwave_surface_albedo_direct,
            )
            if env_bool("RH_CLOUD_OVERLAP_SHORTWAVE", false) &&
               "cloud_fraction" in String.(collect(keys(dataset)))
                cloudy_shortwave_for_overlap =
                    env_bool("RH_CLOUD_OVERLAP_USE_CLOUDY_REGION_OPTICS", false) ?
                    build_cloudy_region_shortwave(clear_shortwave_for_overlap, dataset, j) :
                    shortwave
                overlap_optics = ShortwaveCloudOverlapOpticalProperties(
                    clear_shortwave_for_overlap,
                    cloudy_shortwave_for_overlap,
                    Array(dataset["cloud_fraction"][:, j]);
                    overlap_parameter = "overlap_param" in String.(collect(keys(dataset))) ?
                        Array(dataset["overlap_param"][:, j]) : nothing,
                    fractional_std = "fractional_std" in String.(collect(keys(dataset))) ?
                        Array(dataset["fractional_std"][:, j]) : nothing,
                )
                radiative_fluxes!(
                    fluxes,
                    CloudOverlapShortwave(
                        clear_solver = CloudlessShortwave(
                            rayleigh_backscatter_fraction =
                                env_float("RH_SW_RAYLEIGH_BACKSCATTER_FRACTION", 0.5),
                        ),
                        overlap = Symbol(get(ENV, "RH_CLOUD_OVERLAP_RULE", "maximum")),
                        cloud_fraction_exponent =
                            env_float("RH_CLOUD_FRACTION_EXPONENT", 1.0),
                        inhomogeneity_overlap_exponent =
                            env_float("RH_SW_CLOUD_INHOM_OVERLAP_EXPONENT",
                                      env_float("RH_CLOUD_INHOM_OVERLAP_EXPONENT", 2.0)),
                    ),
                    overlap_optics,
                    atmosphere,
                    shortwave_boundary,
                )
            else
                radiative_fluxes!(
                    fluxes,
                    CloudlessShortwave(
                        rayleigh_backscatter_fraction =
                            env_float("RH_SW_RAYLEIGH_BACKSCATTER_FRACTION", 0.5),
                    ),
                    shortwave,
                    atmosphere,
                    shortwave_boundary,
                )
            end
            heating_rates!(
                column_heating,
                fluxes,
                atmosphere;
                gravity = 9.80665,
                heat_capacity = 1004.0,
            )

            lw_up[:, j] = fluxes.longwave_up
            lw_down[:, j] = fluxes.longwave_down
            sw_up[:, j] = fluxes.shortwave_up
            sw_down[:, j] = fluxes.shortwave_down
            heating[:, j] = 86400.0 .* column_heating
        end

        write_variable!(dataset, CANDIDATE_PREFIX * "lw_up", lw_up, ("interface", "column"))
        write_variable!(dataset, CANDIDATE_PREFIX * "lw_down", lw_down, ("interface", "column"))
        write_variable!(dataset, CANDIDATE_PREFIX * "sw_up", sw_up, ("interface", "column"))
        write_variable!(dataset, CANDIDATE_PREFIX * "sw_down", sw_down, ("interface", "column"))
        write_variable!(dataset, CANDIDATE_PREFIX * "heating_rate", heating, ("layer", "column"))
        write_variable!(dataset, CANDIDATE_PREFIX * "lw_up_clear", clear_lw_up, ("interface", "column"))
        write_variable!(dataset, CANDIDATE_PREFIX * "lw_down_clear", clear_lw_down, ("interface", "column"))
        write_variable!(dataset, CANDIDATE_PREFIX * "sw_up_clear", clear_sw_up, ("interface", "column"))
        write_variable!(dataset, CANDIDATE_PREFIX * "sw_down_clear", clear_sw_down, ("interface", "column"))
        write_variable!(dataset, CANDIDATE_PREFIX * "heating_rate_clear", clear_heating, ("layer", "column"))
        dataset.attrib["radiative_heating_candidate"] =
            "Lightflux staged EcCKDGasOpticsModel + CloudlessLongwave/CloudlessShortwave"
        dataset.attrib["radiative_heating_candidate_gas_optics_mode"] = candidate_mode()
        dataset.attrib["radiative_heating_candidate_gas_amount_units"] = gas_amounts.units
        dataset.attrib["radiative_heating_candidate_gases"] =
            join(string.(Lightflux.gas_names(gas_optics)), ",")
        dataset.attrib["radiative_heating_longwave_interface_sources"] =
            string(use_interface_sources)
        dataset.attrib["radiative_heating_candidate_cos_zenith"] =
            "read from cos_solar_zenith_angle when available, otherwise inferred from top-interface sw_down / solar_irradiance"
        dataset.attrib["radiative_heating_cloud_optics"] =
            "liquid/ice layer cloud optical depths enabled by RH_LIQUID_CLOUD_LW_MASS_ABSORPTION, RH_ICE_CLOUD_LW_MASS_ABSORPTION, RH_LIQUID_CLOUD_SW_MASS_EXTINCTION, RH_ICE_CLOUD_SW_MASS_EXTINCTION, RH_LIQUID_CLOUD_SW_SINGLE_SCATTERING_ALBEDO, RH_ICE_CLOUD_SW_SINGLE_SCATTERING_ALBEDO, RH_LIQUID_CLOUD_SW_SCATTERING_ASYMMETRY, and RH_ICE_CLOUD_SW_SCATTERING_ASYMMETRY"
        dataset.attrib["radiative_heating_cloud_effective_radius_optics"] =
            string(env_bool("RH_CLOUD_EFFECTIVE_RADIUS_OPTICS", false))
        dataset.attrib["radiative_heating_cloud_scattering_table_optics"] =
            string(env_bool("RH_CLOUD_SCATTERING_TABLE_OPTICS", false))
        dataset.attrib["radiative_heating_cloud_scattering_mapping_method"] =
            get(ENV, "RH_CLOUD_SCATTERING_MAPPING_METHOD", "ecrad")
        dataset.attrib["radiative_heating_cloud_scattering_delta_eddington_average"] =
            string(env_bool("RH_CLOUD_SCATTERING_DELTA_EDDINGTON_AVERAGE", true))
        dataset.attrib["radiative_heating_cloud_scattering_thick_averaging"] =
            string(env_bool("RH_CLOUD_SCATTERING_THICK_AVERAGING", true))
        dataset.attrib["radiative_heating_cloud_scattering_delta_eddington_scale"] =
            string(env_bool("RH_CLOUD_SCATTERING_DELTA_EDDINGTON_SCALE", true))
        dataset.attrib["radiative_heating_cloud_scattering_sw_ssa_scale"] =
            string(env_float("RH_CLOUD_SCATTERING_SW_SSA_SCALE", 1.0))
        dataset.attrib["radiative_heating_cloud_overlap_longwave"] =
            string(env_bool("RH_CLOUD_OVERLAP_LONGWAVE", false))
        dataset.attrib["radiative_heating_cloud_overlap_longwave_rule"] =
            get(ENV, "RH_CLOUD_OVERLAP_LONGWAVE_RULE", "adding")
        dataset.attrib["radiative_heating_cloud_overlap_shortwave"] =
            string(env_bool("RH_CLOUD_OVERLAP_SHORTWAVE", false))
        dataset.attrib["radiative_heating_cloud_overlap_use_cloudy_region_optics"] =
            string(env_bool("RH_CLOUD_OVERLAP_USE_CLOUDY_REGION_OPTICS", false))
        dataset.attrib["radiative_heating_cloud_overlap_rule"] =
            get(ENV, "RH_CLOUD_OVERLAP_RULE", "maximum")
        dataset.attrib["radiative_heating_cloud_fraction_exponent"] =
            string(env_float("RH_CLOUD_FRACTION_EXPONENT", 1.0))
        dataset.attrib["radiative_heating_sw_cloud_inhom_overlap_exponent"] =
            string(env_float("RH_SW_CLOUD_INHOM_OVERLAP_EXPONENT",
                             env_float("RH_CLOUD_INHOM_OVERLAP_EXPONENT", 2.0)))
        dataset.attrib["radiative_heating_lw_cloud_inhom_overlap_exponent"] =
            string(env_float("RH_LW_CLOUD_INHOM_OVERLAP_EXPONENT",
                             env_float("RH_CLOUD_INHOM_OVERLAP_EXPONENT", 2.0)))
        dataset.attrib["radiative_heating_aerosol_optics"] =
            string(env_bool("RH_AEROSOL_OPTICS", false))
        dataset.attrib["radiative_heating_aerosol_delta_eddington_scale"] =
            string(env_bool("RH_AEROSOL_DELTA_EDDINGTON_SCALE", true))
        dataset.attrib["radiative_heating_candidate_written_at"] = string(Dates.now())
    end
    return path
end

function write_candidates()
    outputs = String[]
    for case in REQUIRED_CASES
        isfile(reference_path(case.path)) || continue
        push!(outputs, run_candidate_for_file!(case.path))
    end
    return outputs
end

function write_ecrad_candidates_main()
    outputs = write_candidates()
    println("# Wrote ecRad Candidate Variables")
    for output in outputs
        println("- ", output)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    write_ecrad_candidates_main()
end
