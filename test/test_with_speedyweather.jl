using SpeedyWeather
const SpeedyExt = Base.get_extension(AnalyticBandRadiation,
                                     :AnalyticBandRadiationSpeedyWeatherExt)

@testset "Model initializes and runs with SpeedyWeather" begin
    # Basic smoke test: construct, initialize, run one column.
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    @testset for do_CO₂ in (false, true)
        rad   = SpeedyExt.SpeedyAnalyticBandLongwave(spectral_grid; do_CO₂)
        model = PrimitiveWetModel(spectral_grid; longwave_radiation=rad)

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

        ij = rand(1:spectral_grid.npoints)
        SpeedyWeather.parameterization!(ij, vars, model.longwave_radiation, model)

        # After one call, temperature tendency should be non-zero (atmosphere cools)
        @test any(!=(zero(spectral_grid.NF)), vars.tendencies.grid.temperature[ij, :])
    end
end