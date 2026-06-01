@testset "ecCKD-style forward gas optics" begin
    nlayers = 3
    atmosphere = ColumnAtmosphere(
        pressure_layers = [20_000.0, 50_000.0, 80_000.0],
        pressure_interfaces = [10_000.0, 35_000.0, 65_000.0, 95_000.0],
        temperature_layers = [220.0, 260.0, 300.0],
        temperature_interfaces = [210.0, 240.0, 280.0, 305.0],
        gases = (
            h2o = [1.0, 2.0, 3.0],
            co2 = 4.0,
        ),
        surface = (;),
        geometry = (;),
    )

    model = EcCKDGasOpticsModel(
        gas_names = (:h2o, :co2),
        longwave_absorption = [0.1 0.01;
                               0.2 0.02],
        shortwave_absorption = [0.03 0.003],
        longwave_source_scale = [0.5, 1.0],
        longwave_weights = [0.25, 0.75],
        shortwave_weights = [1.0],
    )

    longwave = LongwaveOpticalProperties(zeros(2, nlayers), zeros(2, nlayers);
                                         weights = zeros(2))
    shortwave = ShortwaveOpticalProperties(zeros(1, nlayers); weights = zeros(1))

    returned = optical_properties!(longwave, shortwave, model, atmosphere)
    @test returned == (longwave, shortwave)

    @test longwave.optical_depth ≈ [0.14 0.24 0.34;
                                    0.28 0.48 0.68]
    @test shortwave.optical_depth ≈ [0.042 0.072 0.102]
    @test longwave.weights == [0.25, 0.75]
    @test shortwave.weights == [1.0]

    stefan_boltzmann = 5.670374419e-8
    @test longwave.source[1, :] ≈ 0.5 .* stefan_boltzmann .* atmosphere.temperature_layers .^ 4
    @test longwave.source[2, :] ≈ stefan_boltzmann .* atmosphere.temperature_layers .^ 4

    optical_properties!(longwave, shortwave, model, atmosphere)
    @test (@allocated optical_properties!(longwave, shortwave, model, atmosphere)) == 0

    bad_longwave = LongwaveOpticalProperties(zeros(1, nlayers), zeros(1, nlayers);
                                             weights = zeros(1))
    @test_throws DimensionMismatch optical_properties!(bad_longwave, shortwave, model, atmosphere)

    @test_throws DimensionMismatch EcCKDGasOpticsModel(
        gas_names = (:h2o,),
        longwave_absorption = [0.1 0.01],
        shortwave_absorption = reshape([0.03], 1, 1),
    )
end

@testset "ecCKD-style tabulated longwave source table" begin
    atmosphere = ColumnAtmosphere(
        pressure_layers = [15_000.0],
        pressure_interfaces = [10_000.0, 20_000.0],
        temperature_layers = [275.0],
        temperature_interfaces = [250.0, 300.0],
        gases = (h2o = [1.0], co2 = [1.0]),
        surface = (;),
        geometry = (;),
    )
    model = EcCKDTabulatedGasOpticsModel(
        gas_names = (:h2o, :co2),
        pressure_grid = [10_000.0, 20_000.0],
        temperature_grid = [250.0, 300.0],
        longwave_absorption = fill(0.1, 2, 2, 2, 2),
        shortwave_absorption = fill(0.01, 1, 2, 2, 2),
        longwave_source_temperature_grid = [250.0, 300.0],
        longwave_source_table = [10.0 20.0;
                                 100.0 200.0],
        longwave_weights = [0.5, 0.5],
        shortwave_weights = [1.0],
    )
    longwave = LongwaveOpticalProperties(zeros(2, 1), zeros(2, 1);
                                         weights = zeros(2))
    shortwave = ShortwaveOpticalProperties(zeros(1, 1); weights = zeros(1))

    optical_properties!(longwave, shortwave, model, atmosphere)

    @test longwave.source[:, 1] ≈ [15.0, 150.0]
end

@testset "ecCKD-style tabulated shortwave Rayleigh channel" begin
    atmosphere = ColumnAtmosphere(
        pressure_layers = [15_000.0],
        pressure_interfaces = [10_000.0, 20_000.0],
        temperature_layers = [275.0],
        temperature_interfaces = [250.0, 300.0],
        gases = (h2o = [1.0], co2 = [1.0]),
        surface = (;),
        geometry = (;),
    )
    model = EcCKDTabulatedGasOpticsModel(
        gas_names = (:h2o, :co2),
        pressure_grid = [10_000.0, 20_000.0],
        temperature_grid = [250.0, 300.0],
        longwave_absorption = fill(0.1, 1, 2, 2, 2),
        shortwave_absorption = fill(0.01, 2, 2, 2, 2),
        shortwave_rayleigh_molar_scattering = [1.0e-8, 2.0e-8],
        longwave_weights = [1.0],
        shortwave_weights = [0.5, 0.5],
    )
    longwave = LongwaveOpticalProperties(zeros(1, 1), zeros(1, 1);
                                         weights = zeros(1))
    shortwave = ShortwaveOpticalProperties(zeros(2, 1); weights = zeros(2))

    optical_properties!(longwave, shortwave, model, atmosphere)

    @test shortwave.rayleigh_optical_depth[1, 1] > 0
    @test shortwave.rayleigh_optical_depth[2, 1] ≈ 2shortwave.rayleigh_optical_depth[1, 1]
end

@testset "ecCKD-style tabulated gas optics" begin
    nlayers = 2
    pressure_grid = [10_000.0, 20_000.0]
    temperature_grid = [250.0, 300.0]

    longwave_table = zeros(2, 2, 2, 2)
    shortwave_table = zeros(1, 2, 2, 2)
    for ig in axes(longwave_table, 1), j in axes(longwave_table, 2),
        ip in axes(longwave_table, 3), it in axes(longwave_table, 4)
        longwave_table[ig, j, ip, it] =
            100ig + 10j + 0.001pressure_grid[ip] + 0.01temperature_grid[it]
    end
    for ig in axes(shortwave_table, 1), j in axes(shortwave_table, 2),
        ip in axes(shortwave_table, 3), it in axes(shortwave_table, 4)
        shortwave_table[ig, j, ip, it] =
            10ig + j + 0.0001pressure_grid[ip] + 0.001temperature_grid[it]
    end

    atmosphere = ColumnAtmosphere(
        pressure_layers = [15_000.0, 20_000.0],
        pressure_interfaces = [10_000.0, 17_500.0, 25_000.0],
        temperature_layers = [275.0, 250.0],
        temperature_interfaces = [260.0, 280.0, 245.0],
        gases = (
            h2o = [2.0, 3.0],
            co2 = 4.0,
        ),
        surface = (;),
        geometry = (;),
    )

    model = EcCKDTabulatedGasOpticsModel(
        gas_names = (:h2o, :co2),
        pressure_grid = pressure_grid,
        temperature_grid = temperature_grid,
        longwave_absorption = longwave_table,
        shortwave_absorption = shortwave_table,
        longwave_source_scale = [1.0, 2.0],
        longwave_weights = [0.4, 0.6],
        shortwave_weights = [1.0],
    )

    longwave = LongwaveOpticalProperties(zeros(2, nlayers), zeros(2, nlayers);
                                         weights = zeros(2))
    shortwave = ShortwaveOpticalProperties(zeros(1, nlayers); weights = zeros(1))

    optical_properties!(longwave, shortwave, model, atmosphere)

    lw_coeff(ig, j, p, t) = 100ig + 10j + 0.001p + 0.01t
    sw_coeff(ig, j, p, t) = 10ig + j + 0.0001p + 0.001t
    @test longwave.optical_depth[1, 1] ≈
          lw_coeff(1, 1, 15_000.0, 275.0) * 2.0 +
          lw_coeff(1, 2, 15_000.0, 275.0) * 4.0
    @test longwave.optical_depth[2, 2] ≈
          lw_coeff(2, 1, 20_000.0, 250.0) * 3.0 +
          lw_coeff(2, 2, 20_000.0, 250.0) * 4.0
    @test shortwave.optical_depth[1, 1] ≈
          sw_coeff(1, 1, 15_000.0, 275.0) * 2.0 +
          sw_coeff(1, 2, 15_000.0, 275.0) * 4.0
    @test longwave.weights == [0.4, 0.6]
    @test shortwave.weights == [1.0]

    stefan_boltzmann = 5.670374419e-8
    @test longwave.source[1, :] ≈ stefan_boltzmann .* atmosphere.temperature_layers .^ 4
    @test longwave.source[2, :] ≈ 2 .* stefan_boltzmann .* atmosphere.temperature_layers .^ 4

    optical_properties!(longwave, shortwave, model, atmosphere)
    @test (@allocated optical_properties!(longwave, shortwave, model, atmosphere)) == 0

    @test_throws DimensionMismatch EcCKDTabulatedGasOpticsModel(
        gas_names = (:h2o, :co2),
        pressure_grid = pressure_grid,
        temperature_grid = temperature_grid,
        longwave_absorption = zeros(2, 2, 1, 2),
        shortwave_absorption = shortwave_table,
    )
end
