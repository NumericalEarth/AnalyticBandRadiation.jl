function _test_sw_column(::Type{NF}, nlayers; q = 0.005) where NF
    σ_half = collect(NF.(range(0, 1, length = nlayers + 1)))
    geom   = ColumnGeometry(σ_half)
    T      = NF.(collect(220 .+ 9 .* (0:(nlayers - 1))))
    qv     = fill(NF(q), nlayers)
    Φ      = zeros(NF, nlayers)
    profile = ColumnProfile(temperature = T, humidity = qv, geopotential = Φ,
                            surface_pressure = NF(100_000))
    return profile, geom
end

@testset "TransparentShortwave: surface energy budget" begin
    NF = Float32
    nlayers = 4
    profile, geom = _test_sw_column(NF, nlayers)
    surface = ColumnSurface{NF}(sea_surface_temperature = NF(295),
                                land_surface_temperature = NF(NaN),
                                land_fraction = NF(0),
                                ocean_albedo = NF(0.2),
                                land_albedo = NF(0.15),
                                cos_zenith = NF(0.5))
    constants = PhysicalConstants{NF}()
    thermo = ThermodynamicConstants{NF}()
    dTdt = zeros(NF, nlayers)
    diag = ShortwaveDiagnostics{NF}(nlayers)
    solve_shortwave!(dTdt, diag, TransparentShortwave(), profile, geom, surface, constants, thermo)

    @test diag.surface_shortwave_down ≈ constants.solar_constant * 0.5f0
    @test diag.outgoing_shortwave ≈ 0.2f0 * constants.solar_constant * 0.5f0
    @test all(==(0), dTdt)
end

@testset "OneBandGreyShortwave: absorption in the atmosphere" begin
    NF = Float32
    nlayers = 8
    profile, geom = _test_sw_column(NF, nlayers)
    surface = ColumnSurface{NF}(sea_surface_temperature = NF(295),
                                land_surface_temperature = NF(NaN),
                                land_fraction = NF(0),
                                ocean_albedo = NF(0.1),
                                land_albedo = NF(0.1),
                                cos_zenith = NF(0.5))
    constants = PhysicalConstants{NF}()
    thermo = ThermodynamicConstants{NF}()
    dTdt = zeros(NF, nlayers)
    t_scratch = similar(profile.temperature)
    diag = ShortwaveDiagnostics{NF}(nlayers)
    scheme = OneBandGreyShortwave(NF)
    solve_shortwave!(dTdt, diag, scheme, profile, geom, surface, constants, thermo;
                     transmissivity_scratch = t_scratch)

    D_toa = constants.solar_constant * NF(0.5)
    @test 0 < diag.surface_shortwave_down < D_toa
    @test all(dTdt .> 0)              # heating everywhere
    @test diag.outgoing_shortwave > 0
    @test 0 < diag.albedo < 1
end

@testset "OneBandShortwave (full SPEEDY): runs with diagnostic clouds" begin
    NF = Float32
    nlayers = 8
    profile, geom = _test_sw_column(NF, nlayers; q = 0.01)
    profile = ColumnProfile(temperature = profile.temperature,
                            humidity = profile.humidity,
                            geopotential = profile.geopotential,
                            surface_pressure = profile.surface_pressure,
                            rain_rate = NF(1e-6))                 # ~0.09 mm/day
    surface = ColumnSurface{NF}(sea_surface_temperature = NF(295),
                                land_surface_temperature = NF(285),
                                land_fraction = NF(0.3),
                                ocean_albedo = NF(0.07),
                                land_albedo = NF(0.25),
                                cos_zenith = NF(0.6))
    constants = PhysicalConstants{NF}()
    thermo = ThermodynamicConstants{NF}()
    dTdt = zeros(NF, nlayers)
    t_scratch = similar(profile.temperature)
    diag = ShortwaveDiagnostics{NF}(nlayers)
    scheme = OneBandShortwave(NF)
    solve_shortwave!(dTdt, diag, scheme, profile, geom, surface, constants, thermo;
                     transmissivity_scratch = t_scratch)

    @test all(isfinite, dTdt)
    @test diag.surface_shortwave_down > 0
    @test diag.outgoing_shortwave > 0
    @test 0 <= diag.cloud_cover <= 1
end
