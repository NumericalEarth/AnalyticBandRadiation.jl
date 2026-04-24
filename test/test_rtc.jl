@testset "RadiativeTransferColumn: umbrella API" begin
    nlayers = 8
    σ_half = collect(range(0.0, 1.0, length = nlayers + 1))
    grid    = ColumnGrid(σ_half)
    profile = AtmosphereProfile(
        temperature      = collect(range(220.0, 295.0, length = nlayers)),
        humidity         = fill(0.005, nlayers),
        geopotential     = zeros(nlayers),
        surface_pressure = 100_000.0,
    )
    surface = SurfaceState(
        sea_surface_temperature  = 295.0,
        land_surface_temperature = 285.0,
        land_fraction            = 0.3,
        ocean_albedo             = 0.07,
        land_albedo              = 0.25,
        cos_zenith               = 0.5,
    )

    rtm = RadiativeTransferColumn(; grid, profile, surface)

    solve_longwave!(rtm)
    @test isfinite(rtm.longwave_diagnostics.outgoing_longwave)
    @test rtm.longwave_diagnostics.outgoing_longwave > 0
    @test rtm.longwave_diagnostics.surface_longwave_down >= 0

    solve_shortwave!(rtm)
    @test rtm.shortwave_diagnostics.surface_shortwave_down > 0
    @test rtm.shortwave_diagnostics.outgoing_shortwave > 0
    @test 0 < rtm.shortwave_diagnostics.albedo < 1

    # reset! zeroes tendency + diagnostics
    reset!(rtm)
    @test all(==(0), rtm.temperature_tendency)
    @test rtm.longwave_diagnostics.outgoing_longwave == 0
end

@testset "RadiativeTransferColumn: Float32 propagation" begin
    NF = Float32
    nlayers = 4
    σ_half = collect(NF.(range(0, 1, length = nlayers + 1)))
    grid   = ColumnGrid(σ_half)
    profile = AtmosphereProfile(
        temperature = collect(NF.(range(220, 295, length = nlayers))),
        humidity    = fill(NF(0.005), nlayers),
        geopotential = zeros(NF, nlayers),
        surface_pressure = NF(100_000),
    )
    surface = SurfaceState{NF}(
        sea_surface_temperature  = NF(295),
        land_surface_temperature = NF(285),
        land_fraction            = NF(0.3),
    )
    rtm = RadiativeTransferColumn(; grid, profile, surface,
        longwave_scheme = AnalyticBandLongwave(NF),
        shortwave_scheme = OneBandShortwave(NF),
        physical_constants = PhysicalConstants{NF}(),
        thermodynamic_constants = ThermodynamicConstants{NF}())
    solve_longwave!(rtm)
    @test eltype(rtm.temperature_tendency) === NF
    @test isa(rtm.longwave_diagnostics.outgoing_longwave, NF)
end

@testset "solve_longwave!: duck-typed constants" begin
    nlayers = 4
    σ_half = collect(range(0.0, 1.0, length = nlayers + 1))
    grid   = ColumnGrid(σ_half)
    profile = AtmosphereProfile(
        temperature = collect(range(220.0, 295.0, length = nlayers)),
        humidity = fill(0.005, nlayers),
        geopotential = zeros(nlayers),
        surface_pressure = 100_000.0,
    )
    surface = SurfaceState(sea_surface_temperature = 295.0,
                           land_surface_temperature = 285.0,
                           land_fraction = 0.3)

    # Plain NamedTuple with the expected field names — should work.
    foreign = (gravity = 9.80665, heat_capacity = 1004.64,
               stefan_boltzmann = 5.670374419e-8, solar_constant = 1361.0)
    dTdt = zeros(nlayers)
    diag = LongwaveDiagnostics()
    solve_longwave!(dTdt, diag, AnalyticBandLongwave(), profile, grid, surface, foreign)
    @test diag.outgoing_longwave > 0
end
