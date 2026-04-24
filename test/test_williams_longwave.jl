"Build a σ-coordinate column of `nlayers` layers with lapse-rate temperature,
constant humidity and a 100 kPa surface pressure."
function _test_column(::Type{NF}, nlayers; T_surface = 295, T_top = 220, q = 0.005) where NF
    σ_half = collect(NF.(range(0, 1, length = nlayers + 1)))
    geom   = ColumnGeometry(σ_half)
    T      = NF.(collect(T_top .+ (T_surface - T_top) .* range(0, 1, length = nlayers)))
    q      = fill(NF(q), nlayers)
    Φ      = zeros(NF, nlayers)
    profile = ColumnProfile(temperature = T, humidity = q, geopotential = Φ,
                            surface_pressure = NF(100_000))
    return profile, geom
end

@testset "WilliamsLongwave: parameterization smoke test" begin
    NF = Float32
    nlayers = 8
    for do_co2 in (false, true)
        lw = WilliamsLongwave(NF; do_co2)
        profile, geom = _test_column(NF, nlayers)
        surface = ColumnSurface{NF}(sea_surface_temperature = NF(295),
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

@testset "WilliamsLongwave: energy conservation sanity" begin
    NF = Float32
    nlayers = 8
    lw = WilliamsLongwave(NF; do_co2 = true)
    profile, geom = _test_column(NF, nlayers)
    surface = ColumnSurface{NF}(sea_surface_temperature = NF(295),
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

@testset "WilliamsLongwave: CO₂ forcing" begin
    NF = Float32
    nlayers = 8
    function _olr(co2_ppmv)
        lw = WilliamsLongwave(NF; do_co2 = true, co2_ppmv = NF(co2_ppmv))
        profile, geom = _test_column(NF, nlayers)
        surface = ColumnSurface{NF}(sea_surface_temperature = NF(295),
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

@testset "WilliamsLongwave: Float32 vs Float64 compatibility" begin
    for NF in (Float32, Float64)
        nlayers = 4
        lw = WilliamsLongwave(NF)
        profile, geom = _test_column(NF, nlayers)
        surface = ColumnSurface{NF}(sea_surface_temperature = NF(295),
                                    land_surface_temperature = NF(285),
                                    land_fraction = NF(0.3))
        constants = PhysicalConstants{NF}()
        dTdt = zeros(NF, nlayers)
        diag = LongwaveDiagnostics{NF}()
        solve_longwave!(dTdt, diag, lw, profile, geom, surface, constants)
        @test all(isfinite, dTdt)
    end
end
