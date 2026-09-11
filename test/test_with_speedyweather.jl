using SpeedyWeather, Statistics
const SpeedyExt = Base.get_extension(NumericalRadiation,
                                     :NumericalRadiationSpeedyWeatherExt)

default_spectral_grid() = SpectralGrid(truncation = 16, nlayers = 8)

# Longwave-only model: the analytic-band scheme as the longwave part of the
# Radiation bundle, no shortwave, and no other parameterizations.
function longwave_only_model(spectral_grid; kwargs...)
    longwave = SpeedyExt.SpeedyAnalyticBandLongwave(spectral_grid)
    radiation = Radiation(spectral_grid; shortwave = nothing, longwave)
    model = PrimitiveWetModel(spectral_grid; radiation, parameterizations = (:radiation,), kwargs...)
    initialize!(model.radiation, model)
    return model
end

# Idealised column state: temperature increasing towards the surface
# (k = 1 is the top layer), moist, warm ocean and land surfaces.
function set_test_state!(vars, model)
    nlayers = model.spectral_grid.nlayers
    vars.parameterizations.surface_pressure .= 100000        # [Pa]
    for k in 1:nlayers
        vars.grid.temperature[:, k, :] .= 220 + 9 * (k - 1)   # all time steps
        vars.grid.humidity[:, k, :] .= 0.005
    end
    vars.prognostic.ocean.sea_surface_temperature .= 295
    vars.prognostic.land.soil_temperature .= 285
    vars.tendencies.grid.temperature .= 0
    return nothing
end

@testset "Model initializes and runs with SpeedyWeather" begin
    spectral_grid = default_spectral_grid()
    model = longwave_only_model(spectral_grid)
    @test model.radiation.longwave isa SpeedyExt.SpeedyAnalyticBandLongwave

    vars = Variables(model)
    set_test_state!(vars, model)

    # run the parameterizations, here just longwave
    SpeedyWeather.column_parameterizations!(vars, model)

    # After one call, temperature tendency should be non-zero (atmosphere cools)
    @test any(!=(zero(spectral_grid.NF)), vars.tendencies.grid.temperature)
    @test all(isfinite, vars.tendencies.grid.temperature)
    @test all(>(zero(spectral_grid.NF)), vars.parameterizations.outgoing_longwave)
end

@testset "SpeedyWeather runs CO2 forcing" begin
    # Two prescribed CO₂ concentrations: more CO₂ → less OLR, different tendency.
    spectral_grid = default_spectral_grid()
    NF = spectral_grid.NF

    co2 = CO2(spectral_grid, 280)
    model = longwave_only_model(spectral_grid; greenhouse_gases = (; co2 = co2))
    vars = Variables(model)

    # --- Run 1: 280 ppm CO₂ ---
    set_test_state!(vars, model)
    vars.prognostic.greenhouse_gases.co2[] = 280
    SpeedyWeather.column_parameterizations!(vars, model)
    dT1  = copy(vars.tendencies.grid.temperature)
    olr1 = copy(vars.parameterizations.outgoing_longwave)

    # --- Run 2: 600 ppm CO₂ ---
    set_test_state!(vars, model)
    vars.prognostic.greenhouse_gases.co2[] = 600
    SpeedyWeather.column_parameterizations!(vars, model)
    dT2  = copy(vars.tendencies.grid.temperature)
    olr2 = copy(vars.parameterizations.outgoing_longwave)

    @test all(isfinite, dT1)
    @test all(isfinite, dT2)
    @test any(!=(zero(NF)), dT1)
    @test any(!=(zero(NF)), dT2)
    @test dT1 != dT2

    @test all(isfinite, olr1)
    @test all(>(zero(NF)), olr1)
    @test all(>(zero(NF)), olr2)
    @test all(olr2 .< olr1)
    @test mean(olr2) < mean(olr1)
end

@testset "Full model time steps with the analytic-band longwave" begin
    # Default wet model with only the longwave scheme swapped; a few steps run through.
    spectral_grid = default_spectral_grid()
    longwave = SpeedyExt.SpeedyAnalyticBandLongwave(spectral_grid)
    model = PrimitiveWetModel(spectral_grid; radiation = Radiation(spectral_grid; longwave))
    simulation = initialize!(model)
    run!(simulation, steps = 4)
    @test all(isfinite, simulation.variables.parameterizations.outgoing_longwave)
    @test all(isfinite, simulation.variables.prognostic.temperature)
end

# -----------------------------------------------------------------------------
# EcCKDRadiation: clear-sky ecCKD as one SpeedyWeather radiation component
# -----------------------------------------------------------------------------
using NCDatasets   # activates the NetCDF loader for the reference ecCKD tables

const σ_SB = 5.670374419e-8

@testset "EcCKDRadiation construction and variables" begin
    spectral_grid = default_spectral_grid()
    NF = spectral_grid.NF
    radiation = SpeedyExt.EcCKDRadiation(spectral_grid, "32x32")
    @test radiation isa SpeedyWeather.AbstractRadiation
    @test eltype(radiation.gas_optics) === NF
    @test NumericalRadiation.gas_names(radiation.gas_optics) == (:composite, :h2o, :o3, :co2)

    model = PrimitiveWetModel(spectral_grid; radiation)
    @test model.radiation === radiation
    vars = Variables(model)
    W = vars.parameterizations.ecckd
    nlayers, npoints = spectral_grid.nlayers, spectral_grid.npoints
    @test size(W.longwave_optical_depth) == (npoints, nlayers, 32)
    @test size(W.shortwave_optical_depth) == (npoints, nlayers, 32)
    @test size(W.gas_amounts) == (npoints, nlayers, 4)
    @test size(W.longwave_up) == (npoints, nlayers + 1)
    @test haskey(vars.parameterizations, :outgoing_longwave)
    @test haskey(vars.parameterizations, :outgoing_shortwave)

    # a model gas without a mole fraction is rejected at construction
    gas_optics = read_reference_ecckd_gas_optics("32x32"; names = (:composite, :h2o, :co2, :ch4))
    @test_throws ArgumentError SpeedyExt.EcCKDRadiation(spectral_grid, gas_optics)
    @test SpeedyExt.EcCKDRadiation(spectral_grid, gas_optics; mole_fractions = (; ch4 = 1.8e-6)) isa
          SpeedyExt.EcCKDRadiation
end

@testset "EcCKDRadiation column physics" begin
    spectral_grid = default_spectral_grid()
    NF = spectral_grid.NF
    radiation = SpeedyExt.EcCKDRadiation(spectral_grid, "32x32")
    model = PrimitiveWetModel(spectral_grid; radiation, parameterizations = (:radiation,))
    initialize!(model.radiation, model)
    vars = Variables(model)
    nlayers = spectral_grid.nlayers

    # a realistic column: tropospheric lapse rate, exponentially decaying humidity
    vars.parameterizations.surface_pressure .= 100_000
    for k in 1:nlayers
        σ = model.geometry.σ_levels_full[k]
        vars.grid.temperature[:, k, :] .= max(210, 290 * σ^0.28)
        vars.grid.humidity[:, k, :] .= 0.012 * σ^3
    end
    vars.prognostic.ocean.sea_surface_temperature .= 293
    vars.prognostic.land.soil_temperature .= 288
    vars.parameterizations.ocean.albedo .= 0.06
    vars.parameterizations.land.albedo .= 0.3
    vars.parameterizations.cos_zenith .= 0.5
    vars.tendencies.grid.temperature .= 0

    SpeedyWeather.column_parameterizations!(vars, model)

    P = vars.parameterizations
    olr = P.outgoing_longwave
    @test all(isfinite, olr)
    @test all(200 .< olr .< 320)                                  # W/m², clear sky, ~290 K surface
    @test all(P.surface_longwave_down .< P.surface_longwave_up)    # surface loses longwave
    @test all(P.surface_longwave_up .< σ_SB * 293^4)               # emissivity < 1 and land colder
    @test all(P.outgoing_shortwave .> 0)                            # Rayleigh + surface reflection
    toa_down = model.planet.solar_constant * 0.5
    @test all(0.5 * toa_down .< P.surface_shortwave_down .< toa_down)
    @test all(P.outgoing_shortwave .< P.surface_shortwave_down)

    # heating rates are finite and moderate
    dTdt = vars.tendencies.grid.temperature[:, :, 1]
    @test all(isfinite, dTdt)
    @test maximum(abs, dTdt) * 86400 < 20                          # K/day

    # energy conservation: column-integrated heating equals the net flux convergence
    W = P.ecckd
    ij = 1
    net(k) = W.longwave_down[ij, k] - W.longwave_up[ij, k] + W.shortwave_down[ij, k] - W.shortwave_up[ij, k]
    convergence = net(1) - net(nlayers + 1)
    integrated = sum(dTdt[ij, k] * model.atmosphere.heat_capacity *
                     (W.pressure_interfaces[ij, k + 1] - W.pressure_interfaces[ij, k]) / model.planet.gravity
                     for k in 1:nlayers)
    @test integrated ≈ convergence rtol = 1e-3

    # the same column through the staged API directly gives the same fluxes
    p = W.pressure_layers[ij, :]
    p_half = W.pressure_interfaces[ij, :]
    T = Float64.(vars.grid.temperature[ij, :, 1])
    q = Float64.(vars.grid.humidity[ij, :, 1])
    T_surface = 293.0   # no land at ij = 1? use the blended value from the land fraction instead
    f = model.land_sea_mask.land_fraction[ij]
    T_surface = (1 - f) * 293.0 + f * 288.0
    T_half = zeros(nlayers + 1)
    SpeedyExt.interface_temperatures!(T_half, T, Float64.(p), Float64.(p_half), T_surface)
    gas_optics = read_reference_ecckd_gas_optics("32x32"; names = (:composite, :h2o, :o3, :co2))
    amounts = zeros(nlayers, 4)
    SpeedyExt.gas_amounts!(amounts, Val((:composite, :h2o, :o3, :co2)),
                           (mole_fractions = (; o3 = SpeedyExt.default_ozone_profile), molar_mass_dry_air = 0.0289647, molar_mass_water = 0.018015),
                           q, Float64.(p), Float64.(p_half), 280e-6, Float64(model.planet.gravity))
    atmosphere = ColumnAtmosphere(pressure_layers = Float64.(p), pressure_interfaces = Float64.(p_half),
                                  temperature_layers = T, temperature_interfaces = T_half,
                                  gases = (composite = amounts[:, 1], h2o = amounts[:, 2], o3 = amounts[:, 3], co2 = amounts[:, 4]),
                                  surface = (; temperature = T_surface), geometry = (; cos_zenith = 0.5))
    longwave = LongwaveOptics(zeros(32, nlayers), zeros(32, nlayers); source_top = zeros(32, nlayers),
                              source_bottom = zeros(32, nlayers), weights = zeros(32))
    shortwave = ShortwaveOptics(zeros(32, nlayers); rayleigh_optical_depth = zeros(32, nlayers),
                                scattering_asymmetry = zeros(32, nlayers), weights = zeros(32))
    optical_properties!(longwave, shortwave, gas_optics, atmosphere)
    fluxes = RadiativeFluxes(longwave_up = zeros(nlayers + 1), longwave_down = zeros(nlayers + 1),
                             shortwave_up = zeros(nlayers + 1), shortwave_down = zeros(nlayers + 1))
    emission = (1 - f) .* surface_longwave_emission(gas_optics, 293.0; emissivity = 0.98) .+
               f .* surface_longwave_emission(gas_optics, 288.0; emissivity = 0.98)
    radiative_fluxes!(fluxes, CloudlessLongwave(), longwave, atmosphere,
                      LongwaveBoundaryConditions(surface_longwave_up = emission))
    albedo = (1 - f) * 0.06 + f * 0.3
    radiative_fluxes!(fluxes, CloudlessShortwave(), shortwave, atmosphere,
                      ShortwaveBoundaryConditions(toa_shortwave_down = Float64(model.planet.solar_constant) * 0.5,
                                                  surface_albedo = albedo))
    @test W.longwave_up[ij, :] ≈ fluxes.longwave_up rtol = 2e-3
    @test W.longwave_down[ij, :] ≈ fluxes.longwave_down rtol = 2e-3
    @test W.shortwave_down[ij, :] ≈ fluxes.shortwave_down rtol = 2e-3
    @test W.shortwave_up[ij, :] ≈ fluxes.shortwave_up rtol = 2e-3

    # night: no shortwave, longwave unchanged
    olr_day = copy(olr)
    vars.parameterizations.cos_zenith .= 0
    vars.tendencies.grid.temperature .= 0
    SpeedyWeather.column_parameterizations!(vars, model)
    @test all(==(0), P.outgoing_shortwave)
    @test all(==(0), P.surface_shortwave_down)
    @test P.outgoing_longwave == olr_day
    # a clear-sky atmosphere emits more longwave than it absorbs: the mass-weighted
    # column heating (longwave only at night) is negative everywhere
    dTdt_night = vars.tendencies.grid.temperature[:, :, 1]
    column_heating = [sum(dTdt_night[ij, k] * (W.pressure_interfaces[ij, k + 1] - W.pressure_interfaces[ij, k])
                          for k in 1:nlayers) for ij in 1:spectral_grid.npoints]
    @test all(column_heating .< 0)

    # more CO₂, less OLR
    model_co2 = PrimitiveWetModel(spectral_grid; radiation, parameterizations = (:radiation,),
                                  greenhouse_gases = (; co2 = CO2(spectral_grid, 280)))
    vars_co2 = Variables(model_co2)
    for name in (:surface_pressure, :cos_zenith)
        getproperty(vars_co2.parameterizations, name) .= getproperty(vars.parameterizations, name)
    end
    vars_co2.grid.temperature .= vars.grid.temperature
    vars_co2.grid.humidity .= vars.grid.humidity
    vars_co2.prognostic.ocean.sea_surface_temperature .= 293
    vars_co2.prognostic.land.soil_temperature .= 288
    vars_co2.prognostic.greenhouse_gases.co2[] = 280
    SpeedyWeather.column_parameterizations!(vars_co2, model_co2)
    olr_280 = copy(vars_co2.parameterizations.outgoing_longwave)
    vars_co2.prognostic.greenhouse_gases.co2[] = 1120
    SpeedyWeather.column_parameterizations!(vars_co2, model_co2)
    olr_1120 = vars_co2.parameterizations.outgoing_longwave
    @test all(olr_1120 .< olr_280)
    @test 2 < mean(olr_280 .- olr_1120) < 12                       # ~2 doublings ≈ 7 W/m²
end

@testset "EcCKDRadiation runs in the full model" begin
    spectral_grid = default_spectral_grid()
    radiation = SpeedyExt.EcCKDRadiation(spectral_grid, "32x32")
    model = PrimitiveWetModel(spectral_grid; radiation)
    simulation = initialize!(model)
    run!(simulation, steps = 4)
    P = simulation.variables.parameterizations
    @test all(isfinite, P.outgoing_longwave)
    @test all(150 .< P.outgoing_longwave .< 350)
    @test all(isfinite, simulation.variables.prognostic.temperature)
end
