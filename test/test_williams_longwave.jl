"Build a σ-coordinate column of `nlayers` layers with lapse-rate temperature,
constant humidity and a 100 kPa surface pressure."
function _test_column(::Type{NF}, nlayers; T_surface = 295, T_top = 220, q = 0.005) where NF
    σ_half = collect(NF.(range(0, 1, length = nlayers + 1)))
    geom   = ColumnGrid(σ_half)
    T      = NF.(collect(T_top .+ (T_surface - T_top) .* range(0, 1, length = nlayers)))
    q      = fill(NF(q), nlayers)
    Φ      = zeros(NF, nlayers)
    profile = AtmosphereProfile(temperature = T, humidity = q, geopotential = Φ,
                            surface_pressure = NF(100_000))
    return profile, geom
end

@testset "AnalyticBandLongwave: parameterization smoke test" begin
    NF = Float32
    nlayers = 8
    for do_CO₂ in (false, true)
        lw = AnalyticBandLongwave(NF; do_CO₂)
        profile, geom = _test_column(NF, nlayers)
        surface = SurfaceState{NF}(sea_surface_temperature = NF(295),
                                    land_surface_temperature = NF(285),
                                    land_fraction = NF(0.3))
        constants = PhysicalConstants{NF}()
        dTdt = zeros(NF, nlayers)
        diag = LongwaveDiagnostics{NF}()
        solve_longwave!(dTdt, diag, lw, profile, geom, surface, constants)
        @test any(!=(0), dTdt)
        @test isfinite(diag.outgoing_longwave)
        @test diag.outgoing_longwave > 0
    end
end

@testset "AnalyticBandLongwave: energy conservation sanity" begin
    NF = Float32
    nlayers = 8
    lw = AnalyticBandLongwave(NF; do_CO₂ = true)
    profile, geom = _test_column(NF, nlayers)
    surface = SurfaceState{NF}(sea_surface_temperature = NF(295),
                                land_surface_temperature = NF(285),
                                land_fraction = NF(0.3))
    constants = PhysicalConstants{NF}()
    dTdt = zeros(NF, nlayers)
    diag = LongwaveDiagnostics{NF}()
    solve_longwave!(dTdt, diag, lw, profile, geom, surface, constants)

    @test all(isfinite, dTdt)
    @test diag.outgoing_longwave > 0
    @test diag.surface_longwave_down >= 0
    @test diag.surface_longwave_down <= diag.surface_longwave_up
end

@testset "AnalyticBandLongwave: CO₂ forcing" begin
    NF = Float32
    nlayers = 8
    function _olr(CO₂_ppmv)
        lw = AnalyticBandLongwave(NF; do_CO₂ = true, CO₂_ppmv = NF(CO₂_ppmv))
        profile, geom = _test_column(NF, nlayers)
        surface = SurfaceState{NF}(sea_surface_temperature = NF(295),
                                    land_surface_temperature = NF(285),
                                    land_fraction = NF(0.3))
        constants = PhysicalConstants{NF}()
        dTdt = zeros(NF, nlayers)
        diag = LongwaveDiagnostics{NF}()
        solve_longwave!(dTdt, diag, lw, profile, geom, surface, constants)
        return diag.outgoing_longwave
    end
    olr_280 = _olr(280)
    olr_560 = _olr(560)
    @test olr_560 < olr_280
    forcing = olr_280 - olr_560
    @test 1 < forcing < 10
end

@testset "AnalyticBandLongwave: Float32 vs Float64 compatibility" begin
    for NF in (Float32, Float64)
        nlayers = 4
        lw = AnalyticBandLongwave(NF)
        profile, geom = _test_column(NF, nlayers)
        surface = SurfaceState{NF}(sea_surface_temperature = NF(295),
                                    land_surface_temperature = NF(285),
                                    land_fraction = NF(0.3))
        constants = PhysicalConstants{NF}()
        dTdt = zeros(NF, nlayers)
        diag = LongwaveDiagnostics{NF}()
        solve_longwave!(dTdt, diag, lw, profile, geom, surface, constants)
        @test all(isfinite, dTdt)
    end
end
