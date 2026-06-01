@testset "Radiation validation metrics" begin
    candidate_flux = [1.0, 3.0, 6.0]
    reference_flux = [1.0, 1.0, 2.0]
    candidate_heating = [0.10, 0.20, 0.40]
    reference_heating = [0.05, 0.25, 0.20]

    metrics = radiation_error_metrics(
        candidate_flux = candidate_flux,
        reference_flux = reference_flux,
        candidate_heating_rate = candidate_heating,
        reference_heating_rate = reference_heating,
        candidate_toa_flux = 240.5,
        reference_toa_flux = 240.0,
        candidate_surface_flux = 110.0,
        reference_surface_flux = 112.0,
    )

    @test metrics isa RadiationErrorMetrics
    @test metrics.flux_rmse ≈ sqrt((0.0^2 + 2.0^2 + 4.0^2) / 3)
    @test metrics.flux_max_abs == 4.0
    @test metrics.flux_bias == 2.0
    @test metrics.heating_rate_rmse ≈ sqrt((0.05^2 + (-0.05)^2 + 0.20^2) / 3)
    @test metrics.heating_rate_max_abs ≈ 0.20
    @test metrics.heating_rate_bias ≈ (0.05 - 0.05 + 0.20) / 3
    @test metrics.toa_forcing_error == 0.5
    @test metrics.surface_forcing_error == -2.0

    loose = RadiationThresholds(
        flux_rmse = 3.0,
        flux_max_abs = 4.0,
        flux_abs_bias = 2.0,
        heating_rate_rmse = 0.2,
        heating_rate_max_abs = 0.2,
        heating_rate_abs_bias = 0.1,
        toa_forcing_abs_error = 0.5,
        surface_forcing_abs_error = 2.0,
    )
    valid, errors = passes_thresholds(metrics, loose)
    @test valid
    @test isempty(errors)
    @test passes_thresholds(metrics, loose; throw_on_error = true)

    strict = RadiationThresholds(flux_rmse = 1.0, surface_forcing_abs_error = 1.0)
    valid_strict, strict_errors = passes_thresholds(metrics, strict)
    @test !valid_strict
    @test any(contains("flux_rmse"), strict_errors)
    @test any(contains("surface_forcing_abs_error"), strict_errors)
    @test_throws ArgumentError passes_thresholds(metrics, strict; throw_on_error = true)

    @test_throws DimensionMismatch radiation_error_metrics(
        candidate_flux = [1.0],
        reference_flux = [1.0, 2.0],
        candidate_heating_rate = candidate_heating,
        reference_heating_rate = reference_heating,
        candidate_toa_flux = 0.0,
        reference_toa_flux = 0.0,
        candidate_surface_flux = 0.0,
        reference_surface_flux = 0.0,
    )
end

@testset "RadiativeFluxes comparison metrics" begin
    atmosphere = ColumnAtmosphere(
        pressure_layers = [25_000.0, 75_000.0],
        pressure_interfaces = [0.0, 50_000.0, 100_000.0],
        temperature_layers = [250.0, 280.0],
        temperature_interfaces = [240.0, 265.0, 290.0],
        gases = (;),
        surface = (;),
        geometry = (;),
    )
    candidate = RadiativeFluxes(
        longwave_up = [10.0, 12.0, 14.0],
        longwave_down = [0.0, 2.0, 4.0],
        shortwave_up = [1.0, 2.0, 3.0],
        shortwave_down = [100.0, 80.0, 60.0],
    )
    reference = RadiativeFluxes(
        longwave_up = [9.0, 11.0, 13.0],
        longwave_down = [0.0, 1.0, 3.0],
        shortwave_up = [1.0, 1.0, 2.0],
        shortwave_down = [99.0, 79.0, 59.0],
    )

    metrics = radiative_flux_error_metrics(candidate, reference, atmosphere;
                                           gravity = 10.0,
                                           heat_capacity = 1000.0)
    @test metrics isa RadiationErrorMetrics
    @test metrics.flux_rmse ≈ sqrt(10 / 12)
    @test metrics.flux_max_abs == 1.0
    @test metrics.toa_forcing_error == 0.0
    @test metrics.surface_forcing_error == 0.0
    @test isfinite(metrics.heating_rate_rmse)
end
