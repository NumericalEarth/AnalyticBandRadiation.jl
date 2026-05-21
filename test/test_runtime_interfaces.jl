@testset "Runtime interface scaffold" begin
    @test AbstractAtmosphericState isa DataType
    @test AbstractGasOpticsModel isa DataType
    @test AbstractCloudOpticsModel isa DataType
    @test AbstractAerosolOpticsModel isa DataType
    @test AbstractRadiativeTransferSolver isa DataType
    @test AbstractRadiationBackend isa DataType

    nlayers = 4
    pressure_interfaces = collect(range(1_000.0, 100_000.0, length = nlayers + 1))
    pressure_layers = @views (pressure_interfaces[1:end-1] .+ pressure_interfaces[2:end]) ./ 2
    temperature_interfaces = collect(range(210.0, 290.0, length = nlayers + 1))
    temperature_layers = @views (temperature_interfaces[1:end-1] .+ temperature_interfaces[2:end]) ./ 2

    atmosphere = ColumnAtmosphere(
        pressure_layers = pressure_layers,
        pressure_interfaces = pressure_interfaces,
        temperature_layers = temperature_layers,
        temperature_interfaces = temperature_interfaces,
        gases = (; h2o = fill(0.005, nlayers), co2 = 420.0),
        surface = (; temperature = 290.0, emissivity = 1.0),
        geometry = (; cos_zenith = 0.5),
    )

    @test atmosphere isa AbstractAtmosphericState
    @test eltype(atmosphere) === Float64

    fluxes = RadiativeFluxes(
        longwave_up = zeros(nlayers + 1),
        longwave_down = zeros(nlayers + 1),
        shortwave_up = zeros(nlayers + 1),
        shortwave_down = zeros(nlayers + 1),
    )

    @test eltype(fluxes) === Float64
    @test length(fluxes.longwave_up) == nlayers + 1
end

@testset "RadiativeFluxes to heating rates" begin
    nlayers = 2
    pressure_interfaces = [0.0, 50_000.0, 100_000.0]
    pressure_layers = [25_000.0, 75_000.0]
    temperature_interfaces = fill(280.0, nlayers + 1)
    temperature_layers = fill(280.0, nlayers)
    atmosphere = ColumnAtmosphere(
        pressure_layers = pressure_layers,
        pressure_interfaces = pressure_interfaces,
        temperature_layers = temperature_layers,
        temperature_interfaces = temperature_interfaces,
        gases = (;),
        surface = (;),
        geometry = (;),
    )

    fluxes = RadiativeFluxes(
        longwave_up = zeros(nlayers + 1),
        longwave_down = zeros(nlayers + 1),
        shortwave_up = zeros(nlayers + 1),
        shortwave_down = [200.0, 150.0, 150.0],
    )
    heating = zeros(nlayers)

    heating_rates!(heating, fluxes, atmosphere; gravity = 10.0, heat_capacity = 1000.0)

    @test heating[1] ≈ 1.0e-5
    @test heating[2] == 0.0
    column_energy = sum(heating .* diff(pressure_interfaces)) * 1000.0 / 10.0
    @test column_energy ≈ 50.0

    fluxes.longwave_up .= [100.0, 110.0, 110.0]
    fluxes.shortwave_down .= 0.0
    heating_rates!(heating, fluxes, atmosphere; gravity = 10.0, heat_capacity = 1000.0)

    @test heating[1] ≈ 2.0e-6
    @test heating[2] == 0.0

    fluxes.longwave_down .= [10.0, 15.0, 20.0]
    fluxes.longwave_up .= [80.0, 70.0, 65.0]
    fluxes.shortwave_down .= [500.0, 450.0, 410.0]
    fluxes.shortwave_up .= [50.0, 55.0, 60.0]
    heating_rates!(heating, fluxes, atmosphere; gravity = 10.0, heat_capacity = 1000.0)

    net_flux = fluxes.longwave_down .- fluxes.longwave_up .+
        fluxes.shortwave_down .- fluxes.shortwave_up
    layer_energy = heating .* diff(pressure_interfaces) .* 1000.0 ./ 10.0
    @test layer_energy ≈ net_flux[1:end-1] .- net_flux[2:end]
    @test maximum(abs.(layer_energy .- (net_flux[1:end-1] .- net_flux[2:end]))) < 1.0e-12
end

@testset "RadiativeTransferColumn staged wrapper" begin
    nlayers = 6
    grid = ColumnGrid(collect(range(0.0, 1.0, length = nlayers + 1)))
    profile = AtmosphereProfile(
        temperature = collect(range(220.0, 295.0, length = nlayers)),
        humidity = fill(0.005, nlayers),
        geopotential = zeros(nlayers),
        surface_pressure = 100_000.0,
    )
    surface = SurfaceState(
        sea_surface_temperature = 295.0,
        land_surface_temperature = 285.0,
        land_fraction = 0.3,
        ocean_albedo = 0.07,
        land_albedo = 0.25,
        cos_zenith = 0.5,
    )

    rtm = RadiativeTransferColumn(; grid, profile, surface)
    @test radiation_workspace(rtm) === rtm

    returned = radiative_heating!(rtm)
    @test returned === rtm
    @test any(!iszero, rtm.temperature_tendency)
    @test rtm.longwave_diagnostics.outgoing_longwave > 0
    @test rtm.shortwave_diagnostics.surface_shortwave_down > 0

    heating = similar(rtm.temperature_tendency)
    heating_rates!(heating, rtm)
    @test heating == rtm.temperature_tendency

    previous = copy(rtm.temperature_tendency)
    radiative_heating!(rtm; reset = false, shortwave = false)
    @test rtm.temperature_tendency != previous
end
