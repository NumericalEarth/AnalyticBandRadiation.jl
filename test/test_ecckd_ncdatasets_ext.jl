using NCDatasets

@testset "NCDatasets ecCKD reader extension" begin
    path = tempname() * ".nc"

    NCDataset(path, "c") do ds
        defDim(ds, "lw_bands", 2)
        defDim(ds, "sw_bands", 3)
        defDim(ds, "lw_gpoints", 4)
        defDim(ds, "sw_gpoints", 5)
        defDim(ds, "gases", 2)
        defDim(ds, "pressure", 6)
        defDim(ds, "temperature", 7)

        defVar(ds, "lw_absorption", Float64,
               ("gases", "lw_gpoints", "pressure", "temperature"))
        defVar(ds, "sw_absorption", Float64,
               ("gases", "sw_gpoints", "pressure", "temperature"))
        defVar(ds, "lw_source", Float64, ("lw_gpoints", "temperature"))
        defVar(ds, "sw_rayleigh", Float64, ("sw_gpoints", "pressure"))

        ds.attrib["model_name"] = "toy-netcdf-ecCKD"
        ds.attrib["version"] = "0.2"
        ds.attrib["gas_names"] = ["h2o", "co2"]
    end

    definition = read_ecckd_definition(path)
    @test definition isa EcCKDDefinition
    @test validate_ecckd_definition(definition)

    summary = summarize_ecckd_definition(definition)
    @test summary.model_name == "toy-netcdf-ecCKD"
    @test summary.version == "0.2"
    @test summary.lw_bands == 2
    @test summary.sw_bands == 3
    @test summary.lw_gpoints == 4
    @test summary.sw_gpoints == 5
    @test summary.gases == ["h2o", "co2"]
    @test summary.pressure_grid_size == 6
    @test summary.temperature_grid_size == 7
    @test summary.source_tables_present
    @test summary.rayleigh_tables_present
end

@testset "official ecCKD runtime LUT ingestion" begin
    paths = official_ecckd_definition_paths(require = false)
    lw_path = paths.longwave
    sw_path = paths.shortwave

    if lw_path !== nothing && sw_path !== nothing && isfile(lw_path) && isfile(sw_path)
        model = read_ecckd_tabulated_gas_optics(lw_path, sw_path;
                                                gas_names = (:h2o, :co2),
                                                h2o_mole_fraction = 0.005)
        @test model isa EcCKDTabulatedGasOpticsModel
        @test size(model.longwave_absorption) == (64, 2, 53, 6)
        @test size(model.shortwave_absorption) == (32, 2, 53, 6)
        @test length(model.h2o_mole_fraction_grid) == 12
        @test size(model.longwave_h2o_absorption) == (64, 53, 6, 12)
        @test size(model.shortwave_h2o_absorption) == (32, 53, 6, 12)
        @test all(iszero, model.longwave_absorption[:, 1, :, :])
        @test all(iszero, model.shortwave_absorption[:, 1, :, :])
        @test size(model.temperature_grid) == (53, 6)
        @test model.pressure_grid[begin] < model.pressure_grid[end]
        @test all(isfinite, model.longwave_absorption)
        @test all(isfinite, model.shortwave_absorption)
        @test length(model.shortwave_rayleigh_molar_scattering) == 32
        @test maximum(model.shortwave_rayleigh_molar_scattering) > 0
        @test sum(model.longwave_weights) ≈ 1.0
        @test sum(model.shortwave_weights) ≈ 1.0

        atmosphere = ColumnAtmosphere(
            pressure_layers = [20_000.0, 80_000.0],
            pressure_interfaces = [10_000.0, 50_000.0, 100_000.0],
            temperature_layers = [240.0, 290.0],
            temperature_interfaces = [230.0, 265.0, 300.0],
            gases = (
                h2o = [1.0e-3, 5.0e-3],
                co2 = [4.0e-4, 4.0e-4],
            ),
            surface = (;),
            geometry = (;),
        )
        longwave = LongwaveOpticalProperties(zeros(64, 2), zeros(64, 2);
                                             weights = zeros(64))
        shortwave = ShortwaveOpticalProperties(zeros(32, 2);
                                               weights = zeros(32))
        optical_properties!(longwave, shortwave, model, atmosphere)

        @test all(isfinite, longwave.optical_depth)
        @test all(isfinite, shortwave.optical_depth)
        @test all(longwave.optical_depth .>= 0)
        @test all(shortwave.optical_depth .>= 0)
        @test all(isfinite, shortwave.rayleigh_optical_depth)
        @test all(shortwave.rayleigh_optical_depth .>= 0)
        @test maximum(shortwave.rayleigh_optical_depth) > 0
    else
        @info "Skipping official ecCKD runtime LUT ingestion check; ecRad data files are not present" lw_path sw_path
        @test_skip "official ecCKD runtime LUT files are not present"
    end
end

@testset "official ecCKD definition files" begin
    paths = official_ecckd_definition_paths(require = false)
    lw_path = paths.longwave
    sw_path = paths.shortwave

    if lw_path !== nothing && sw_path !== nothing && isfile(lw_path) && isfile(sw_path)
        lw = read_ecckd_definition(lw_path)
        sw = read_ecckd_definition(sw_path)

        @test validate_ecckd_definition(lw)
        @test validate_ecckd_definition(sw)

        lw_summary = summarize_ecckd_definition(lw)
        sw_summary = summarize_ecckd_definition(sw)

        @test lw_summary.lw_gpoints == 64
        @test lw_summary.lw_bands == 13
        @test lw_summary.sw_gpoints == 0
        @test lw_summary.source_tables_present
        @test "h2o" in lw_summary.gases
        @test "co2" in lw_summary.gases

        @test sw_summary.sw_gpoints == 32
        @test sw_summary.sw_bands == 5
        @test sw_summary.lw_gpoints == 0
        @test sw_summary.rayleigh_tables_present
        @test "h2o" in sw_summary.gases
        @test "co2" in sw_summary.gases
    else
        @info "Skipping official ecCKD definition file checks; ecRad data files are not present" lw_path sw_path
        @test_skip "official ecCKD definition files are not present"
    end
end
