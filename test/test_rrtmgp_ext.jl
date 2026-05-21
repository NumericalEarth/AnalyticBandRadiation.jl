using AnalyticBandRadiation
using ClimaComms
using NCDatasets
using RRTMGP

@testset "RRTMGP extension direct ColumnAtmosphere comparison" begin
    ext = Base.get_extension(AnalyticBandRadiation, :AnalyticBandRadiationRRTMGPExt)
    @test ext !== nothing

    pressure_interfaces = [1.0e3, 2.0e4, 5.0e4, 8.0e4, 1.0e5]
    pressure_layers = (pressure_interfaces[1:end-1] .+ pressure_interfaces[2:end]) ./ 2
    temperature_interfaces = [210.0, 235.0, 260.0, 285.0, 300.0]
    temperature_layers = (temperature_interfaces[1:end-1] .+ temperature_interfaces[2:end]) ./ 2
    atmosphere = ColumnAtmosphere(
        pressure_layers = pressure_layers,
        pressure_interfaces = pressure_interfaces,
        temperature_layers = temperature_layers,
        temperature_interfaces = temperature_interfaces,
        gases = (h2o = [1e-4, 8e-4, 3e-3, 1e-2],
                 o3 = fill(1e-8, 4),
                 co2 = 400e-6),
        surface = (;),
        geometry = (;),
    )
    rrtmgp = ext.RRTMGPClearSkyModel(Float64)
    boundary = ext.RRTMGPBoundaryConditions(surface_temperature = 300.0,
                                            surface_emissivity = 0.98,
                                            surface_albedo = 0.1,
                                            toa_shortwave_down = 500.0,
                                            cos_zenith = 1.0)
    rrtmgp_fluxes = RadiativeFluxes(
        longwave_up = zeros(5),
        longwave_down = zeros(5),
        shortwave_up = zeros(5),
        shortwave_down = zeros(5),
    )
    workspace = radiation_workspace(rrtmgp, atmosphere)
    radiative_fluxes!(rrtmgp_fluxes, rrtmgp, atmosphere, boundary, workspace)

    @test all(isfinite, rrtmgp_fluxes.longwave_up)
    @test all(isfinite, rrtmgp_fluxes.longwave_down)
    @test all(isfinite, rrtmgp_fluxes.shortwave_up)
    @test all(isfinite, rrtmgp_fluxes.shortwave_down)

    gas = EcCKDGasOpticsModel(
        gas_names = (:h2o, :co2),
        longwave_absorption = fill(0.01, 2, 2),
        shortwave_absorption = fill(0.005, 2, 2),
        longwave_weights = fill(0.5, 2),
        shortwave_weights = fill(0.5, 2),
        longwave_source_scale = fill(1.0, 2),
    )
    longwave = LongwaveOpticalProperties(zeros(2, 4), zeros(2, 4); weights = fill(0.5, 2))
    shortwave = ShortwaveOpticalProperties(zeros(2, 4); weights = fill(0.5, 2))
    optical_properties!(longwave, shortwave, gas, atmosphere)
    rh_fluxes = RadiativeFluxes(
        longwave_up = zeros(5),
        longwave_down = zeros(5),
        shortwave_up = zeros(5),
        shortwave_down = zeros(5),
    )
    radiative_fluxes!(rh_fluxes,
                      CloudlessLongwave(),
                      longwave,
                      atmosphere,
                      LongwaveBoundaryConditions(surface_longwave_up = 0.98 * 5.670374419e-8 * 300.0^4))
    radiative_fluxes!(rh_fluxes,
                      CloudlessShortwave(),
                      shortwave,
                      atmosphere,
                      ShortwaveBoundaryConditions(toa_shortwave_down = 500.0,
                                                  surface_albedo = 0.1))
    metrics = radiative_flux_error_metrics(rh_fluxes, rrtmgp_fluxes, atmosphere;
                                           gravity = 9.80665,
                                           heat_capacity = 1004.0)
    @test isfinite(metrics.flux_rmse)
    @test isfinite(metrics.heating_rate_rmse)
end
