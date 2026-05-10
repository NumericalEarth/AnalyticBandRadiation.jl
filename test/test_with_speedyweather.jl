using SpeedyWeather
const SpeedyExt = Base.get_extension(AnalyticBandRadiation,
                                     :AnalyticBandRadiationSpeedyWeatherExt)

default_spectral_grid() = SpectralGrid(trunc=15, nlayers=8)

@testset "Model initializes and runs with SpeedyWeather" begin
    # Basic smoke test: construct, initialize, run one step in low resolution
    spectral_grid = default_spectral_grid()

    @testset for do_CO₂ in (false, true)
        rad   = SpeedyExt.SpeedyAnalyticBandLongwave(spectral_grid; do_CO₂)

        # use only longwave radation as the only parameterization
        model = PrimitiveWetModel(spectral_grid; longwave_radiation=rad, parameterizations=(:longwave_radiation,))

        initialize!(model.longwave_radiation, model)
        vars = Variables(model)
        vars.grid.pressure_prev .= 100000
        # k=1 is top layer, k=nlayers is bottom layer → temperature increases with k
        for k in 1:spectral_grid.nlayers
            vars.grid.temperature_prev[:, k] .= 220 + 9 * (k-1)
            vars.grid.humidity_prev[:, k]    .= 0.005
        end
        vars.prognostic.ocean.sea_surface_temperature .= 295
        vars.prognostic.land.soil_temperature[:, 1]   .= 285

        # run the parameterizations, here just longwave 
        SpeedyWeather.column_parameterizations!(vars, model)
 
        # After one call, temperature tendency should be non-zero (atmosphere cools)
        @test any(!=(zero(spectral_grid.NF)), vars.tendencies.grid.temperature)
    end
end

@testset "SpeedyWeather runs CO2 forcing" begin
    # Construct, initialize, and run the column parameterization with two
    # different prescribed CO₂ concentrations. Verify that the tendency and
    # OLR respond in a physically sensible way (more CO₂ → less OLR, less
    # net longwave cooling of the column).
    spectral_grid = default_spectral_grid()
    NF = spectral_grid.NF

    rad = SpeedyExt.SpeedyAnalyticBandLongwave(spectral_grid; do_CO₂=true, CO₂_forcing=true)

    co2 = CO2(spectral_grid, 280)

    model = PrimitiveWetModel(spectral_grid; longwave_radiation=rad, parameterizations=(:longwave_radiation,),
            greenhouse_gases = (; co2 = co2))

    initialize!(model.longwave_radiation, model)
    vars = Variables(model)
    vars.grid.pressure_prev .= 100000
    # k=1 is top layer, k=nlayers is bottom layer → temperature increases with k
    for k in 1:spectral_grid.nlayers
        vars.grid.temperature_prev[:, k] .= 220 + 9 * (k-1)
        vars.grid.humidity_prev[:, k]    .= 0.005
    end
    vars.prognostic.ocean.sea_surface_temperature .= 295
    vars.prognostic.land.soil_temperature[:, 1]   .= 285

    # --- Run 1: 280 ppm CO₂ ---
    vars.prognostic.greenhouse_gases.co2[] = 280
    vars.tendencies.grid.temperature .= 0
    SpeedyWeather.column_parameterizations!(vars, model)
    dT1  = copy(vars.tendencies.grid.temperature)
    olr1 = copy(vars.parameterizations.outgoing_longwave)

    # --- Run 2: 600 ppm CO₂ (much higher) ---
    vars.prognostic.greenhouse_gases.co2[] = 600
    vars.tendencies.grid.temperature .= 0
    SpeedyWeather.column_parameterizations!(vars, model)
    dT2  = copy(vars.tendencies.grid.temperature)
    olr2 = copy(vars.parameterizations.outgoing_longwave)

    # Both runs produced finite, non-trivial tendencies.
    @test all(isfinite, dT1)
    @test all(isfinite, dT2)
    @test any(!=(zero(NF)), dT1)
    @test any(!=(zero(NF)), dT2)

    # Tendencies should differ between the two CO₂ levels (the forcing was applied).
    @test dT1 != dT2

    # OLR is positive and finite, and increasing CO₂ reduces OLR (greenhouse effect).
    @test all(isfinite, olr1)
    @test all(>(zero(NF)), olr1)
    @test all(>(zero(NF)), olr2)
    @test all(olr2 .< olr1)

    # The column's net longwave cooling weakens at higher CO₂, i.e. the
    # column-integrated temperature tendency is less negative (warmer) for
    # 600 ppm than for 280 ppm at every horizontal grid point.
    column_dT1 = sum(dT1, dims=2)
    column_dT2 = sum(dT2, dims=2)
    @test all(column_dT2 .> column_dT1)
end 