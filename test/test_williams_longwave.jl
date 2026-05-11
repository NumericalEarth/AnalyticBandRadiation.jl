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
    lw = AnalyticBandLongwave(NF)
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

@testset "AnalyticBandLongwave: energy conservation sanity" begin
    NF = Float32
    nlayers = 8
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
    @test diag.outgoing_longwave > 0
    @test diag.surface_longwave_down >= 0
    @test diag.surface_longwave_down <= diag.surface_longwave_up
end

@testset "AnalyticBandLongwave: CO₂ forcing" begin
    NF = Float32
    nlayers = 8

    # CO₂ as a compile-time scheme parameter (CO₂_forcing = false).
    function _olr_param(CO₂_concentration)
        lw = AnalyticBandLongwave(NF; CO₂_concentration = NF(CO₂_concentration))
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
    olr_280 = _olr_param(280)
    olr_560 = _olr_param(560)
    @test olr_560 < olr_280
    forcing = olr_280 - olr_560
    @test 1 < forcing < 10

    # CO₂ as an external input via the profile (CO₂_forcing = true).
    # The scheme's CO₂_concentration field is ignored and profile.CO₂_concentration is used instead.
    function _olr_forcing(CO₂_concentration)
        lw = AnalyticBandLongwave(NF; CO₂_forcing = true,
                                   CO₂_concentration = NF(0))   # set to 0 to confirm it is unused
        profile, geom = _test_column(NF, nlayers)
        profile = AtmosphereProfile(temperature = profile.temperature,
                                    humidity = profile.humidity,
                                    geopotential = profile.geopotential,
                                    surface_pressure = profile.surface_pressure,
                                    CO₂_concentration = NF(CO₂_concentration))
        surface = SurfaceState{NF}(sea_surface_temperature = NF(295),
                                    land_surface_temperature = NF(285),
                                    land_fraction = NF(0.3))
        constants = PhysicalConstants{NF}()
        dTdt = zeros(NF, nlayers)
        diag = LongwaveDiagnostics{NF}()
        solve_longwave!(dTdt, diag, lw, profile, geom, surface, constants)
        return diag.outgoing_longwave
    end
    olr_280_f = _olr_forcing(280)
    olr_560_f = _olr_forcing(560)
    @test olr_560_f < olr_280_f
    forcing_f = olr_280_f - olr_560_f
    @test 1 < forcing_f < 10
    # Both code paths should yield the same OLR for the same CO₂ concentration.
    @test olr_280 ≈ olr_280_f
    @test olr_560 ≈ olr_560_f
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
