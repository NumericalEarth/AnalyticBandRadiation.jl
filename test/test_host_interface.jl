using Test
using NumericalRadiation

# Host-model facing preparations (plan Phase 1): in-place surface emission,
# element-type conversion of tabulated models, and ColumnAtmosphere built from
# views of differently shaped host arrays.

# Small tabulated model with every optional table populated, in Float64.
function host_fixture_model()
    np, nt, nh2o = 4, 3, 3
    ng_lw, ng_sw, ngas = 3, 2, 3
    pressure_grid = exp.(range(log(5_000.0), log(100_000.0), length = np))
    temperature_grid = [200.0, 250.0, 300.0]
    source_temperature_grid = [180.0, 240.0, 300.0]
    return EcCKDTabulatedGasOpticsModel(
        names = (:h2o, :co2, :composite),
        pressure_grid = pressure_grid,
        temperature_grid = temperature_grid,
        h2o_mole_fraction_grid = [1e-6, 1e-4, 1e-2],
        gas_reference_mole_fractions = [0.0, 4e-4, 0.0],
        longwave_absorption =
            [1e-4 * (7ig + 3j) * (1 + 1e-5 * pressure_grid[ip]) * (1 + 1e-3 * temperature_grid[it])
             for ig in 1:ng_lw, j in 1:ngas, ip in 1:np, it in 1:nt],
        shortwave_absorption =
            [1e-5 * (5ig + 2j) * (1 + 2e-5 * pressure_grid[ip]) * (1 + 2e-3 * temperature_grid[it])
             for ig in 1:ng_sw, j in 1:ngas, ip in 1:np, it in 1:nt],
        longwave_h2o_absorption =
            [1e-3 * ig * (1 + 10ih) for ig in 1:ng_lw, ip in 1:np, it in 1:nt, ih in 1:nh2o],
        shortwave_h2o_absorption =
            [1e-4 * ig * (1 + 5ih) for ig in 1:ng_sw, ip in 1:np, it in 1:nt, ih in 1:nh2o],
        shortwave_rayleigh_molar_scattering = [1.1e-6, 3.7e-6],
        longwave_source_temperature_grid = source_temperature_grid,
        longwave_source_table = [(ig + 2) * st^2 for ig in 1:ng_lw, st in source_temperature_grid],
        longwave_weights = [0.2, 0.3, 0.5],
        shortwave_weights = [0.45, 0.55],
    )
end

function host_fixture_optics(FT, model, nlayers)
    ng_lw = length(model.longwave_weights)
    ng_sw = length(model.shortwave_weights)
    longwave = LongwaveOptics(zeros(FT, ng_lw, nlayers), zeros(FT, ng_lw, nlayers);
                              source_top = zeros(FT, ng_lw, nlayers),
                              source_bottom = zeros(FT, ng_lw, nlayers),
                              weights = zeros(FT, ng_lw))
    shortwave = ShortwaveOptics(zeros(FT, ng_sw, nlayers);
                                rayleigh_optical_depth = zeros(FT, ng_sw, nlayers),
                                scattering_asymmetry = zeros(FT, ng_sw, nlayers),
                                weights = zeros(FT, ng_sw))
    return longwave, shortwave
end

@testset "surface_longwave_emission! (in place)" begin
    model = host_fixture_model()
    for T in (250.0, 300.0), emissivity in (1.0, 0.95)
        out = zeros(length(model.longwave_weights))
        returned = surface_longwave_emission!(out, model, T; emissivity)
        @test returned === out
        @test out == surface_longwave_emission(model, T; emissivity)
    end
    @test_throws DimensionMismatch surface_longwave_emission!(zeros(2), model, 300.0)

    # Float32 model, Float32 output, no allocation in the hot call
    model32 = EcCKDTabulatedGasOpticsModel{Float32}(model)
    out32 = zeros(Float32, length(model.longwave_weights))
    surface_longwave_emission!(out32, model32, 300.0f0)   # warm up
    @test (@allocated surface_longwave_emission!(out32, model32, 300.0f0)) == 0
    @test out32 ≈ Float32.(surface_longwave_emission(model, 300.0)) rtol = 1e-5
end

@testset "EcCKDTabulatedGasOpticsModel{FT} element-type conversion" begin
    model = host_fixture_model()
    model32 = EcCKDTabulatedGasOpticsModel{Float32}(model)

    @test eltype(model) === Float64
    @test eltype(model32) === Float32
    @test NumericalRadiation.gas_names(model32) == NumericalRadiation.gas_names(model)
    for name in fieldnames(typeof(model))
        a, b = getfield(model, name), getfield(model32, name)
        a === nothing && (@test b === nothing; continue)
        @test eltype(b) === Float32
        @test size(b) == size(a)
        @test b ≈ a rtol = 1e-6
    end

    # converting to the same type reuses the arrays
    same = EcCKDTabulatedGasOpticsModel{Float64}(model)
    @test same.longwave_absorption === model.longwave_absorption
    @test same.pressure_grid === model.pressure_grid

    # the Float32 model produces Float32 optics close to the Float64 ones
    nlayers = 3
    atmosphere(FT) = ColumnAtmosphere(
        pressure_layers = FT[7_500, 33_000, 88_000],
        pressure_interfaces = FT[1_000, 14_000, 52_000, 101_000],
        temperature_layers = FT[217.3, 263.9, 291.4],
        temperature_interfaces = FT[205.1, 231.7, 279.2, 297.8],
        gases = (h2o = FT[3.1, 12.7, 41.9], co2 = FT(8.3), composite = FT[4.2e2, 1.1e3, 2.6e3]),
        surface = (;), geometry = (;))
    lw64, sw64 = host_fixture_optics(Float64, model, nlayers)
    lw32, sw32 = host_fixture_optics(Float32, model32, nlayers)
    optical_properties!(lw64, sw64, model, atmosphere(Float64))
    optical_properties!(lw32, sw32, model32, atmosphere(Float32))
    @test eltype(lw32.optical_depth) === Float32
    @test lw32.optical_depth ≈ lw64.optical_depth rtol = 1e-4
    @test lw32.source ≈ lw64.source rtol = 1e-4
    @test sw32.optical_depth ≈ sw64.optical_depth rtol = 1e-4
    @test sw32.rayleigh_optical_depth ≈ sw64.rayleigh_optical_depth rtol = 1e-4

    # absent optional tables stay absent
    plain = EcCKDTabulatedGasOpticsModel(
        names = (:h2o, :co2),
        pressure_grid = [10_000.0, 100_000.0],
        temperature_grid = [220.0, 300.0],
        longwave_absorption = ones(2, 2, 2, 2),
        shortwave_absorption = ones(2, 2, 2, 2),
    )
    plain32 = EcCKDTabulatedGasOpticsModel{Float32}(plain)
    @test isempty(plain32.longwave_h2o_absorption)
    @test isempty(plain32.h2o_mole_fraction_grid)
    @test plain32.longwave_source_table === nothing
    @test eltype(plain32.longwave_weights) === Float32
end

@testset "ColumnAtmosphere from differently shaped host views" begin
    model = host_fixture_model()
    nlayers = 3
    npoints = 5
    ij = 2

    # host-style storage: layers in a (npoints, nlayers) matrix, interfaces in a
    # (npoints, nlayers + 1) matrix, temperatures with a trailing step dimension
    pressure_layers = zeros(npoints, nlayers)
    pressure_interfaces = zeros(npoints, nlayers + 1)
    temperature_layers = zeros(npoints, nlayers, 2)
    temperature_interfaces = zeros(npoints, nlayers + 1)
    pressure_layers[ij, :] .= [7_500, 33_000, 88_000]
    pressure_interfaces[ij, :] .= [1_000, 14_000, 52_000, 101_000]
    temperature_layers[ij, :, 1] .= [217.3, 263.9, 291.4]
    temperature_interfaces[ij, :] .= [205.1, 231.7, 279.2, 297.8]
    gases = (h2o = [3.1, 12.7, 41.9], co2 = 8.3, composite = [4.2e2, 1.1e3, 2.6e3])

    view_atmosphere = ColumnAtmosphere(
        pressure_layers = view(pressure_layers, ij, :),
        pressure_interfaces = view(pressure_interfaces, ij, :),
        temperature_layers = view(temperature_layers, ij, :, 1),
        temperature_interfaces = view(temperature_interfaces, ij, :),
        gases = gases, surface = (;), geometry = (;))
    vector_atmosphere = ColumnAtmosphere(
        pressure_layers = pressure_layers[ij, :],
        pressure_interfaces = pressure_interfaces[ij, :],
        temperature_layers = temperature_layers[ij, :, 1],
        temperature_interfaces = temperature_interfaces[ij, :],
        gases = gases, surface = (;), geometry = (;))
    @test eltype(view_atmosphere) === Float64

    results = map((view_atmosphere, vector_atmosphere)) do atmosphere
        longwave, shortwave = host_fixture_optics(Float64, model, nlayers)
        optical_properties!(longwave, shortwave, model, atmosphere)
        fluxes = RadiativeFluxes(longwave_up = zeros(nlayers + 1), longwave_down = zeros(nlayers + 1),
                                 shortwave_up = zeros(nlayers + 1), shortwave_down = zeros(nlayers + 1))
        radiative_fluxes!(fluxes, CloudlessLongwave(), longwave, atmosphere,
                          LongwaveBoundaryConditions(surface_longwave_up = surface_longwave_emission(model, 297.8)))
        radiative_fluxes!(fluxes, CloudlessShortwave(), shortwave, atmosphere,
                          ShortwaveBoundaryConditions(toa_shortwave_down = 700.0, surface_albedo = 0.1))
        heating = zeros(nlayers)
        heating_rates!(heating, fluxes, atmosphere; gravity = 9.80665, heat_capacity = 1004.0)
        (longwave.optical_depth, fluxes.longwave_up, fluxes.shortwave_down, heating)
    end
    for (a, b) in zip(results[1], results[2])
        @test a == b
    end
    @test all(isfinite, results[1][4])
end
