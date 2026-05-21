const ABR_ROOT_FOR_ARTIFACT_TEST = normpath(joinpath(@__DIR__, ".."))
if Base.find_package("Lightflux") === nothing
    push!(LOAD_PATH, ABR_ROOT_FOR_ARTIFACT_TEST)
end

using Lightflux

@testset "official ecCKD artifact resolver" begin
    inventory = official_ecckd_model_inventory()
    @test "ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc" in inventory
    @test "ecckd-1.4_sw_climate_vfine-96b_ckd-definition.nc" in inventory

    root = ecrad_data_path(require = true)
    @test isdir(root)

    mktempdir() do temp_root
        withenv("RH_ECRAD_DATA_PATH" => temp_root) do
            @test ecrad_data_path(require = true) == normpath(temp_root)
        end
    end

    mktempdir() do temp_root
        filename = "ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc"
        nested_data = joinpath(temp_root, "ecrad-archive", "data")
        mkpath(nested_data)
        write(joinpath(nested_data, filename), "")
        withenv("RH_ECRAD_DATA_PATH" => temp_root) do
            @test official_ecckd_definition_path(filename; require = true) ==
                  normpath(joinpath(nested_data, filename))
        end
    end

    source_root = ecckd_source_path(require = true)
    @test isfile(joinpath(source_root, "src", "ecckd", "optimize_lut.cpp"))

    paths = official_ecckd_definition_paths(require = true)
    @test isfile(paths.longwave)
    @test isfile(paths.shortwave)
    @test basename(paths.longwave) == "ecckd-1.2_lw_climate_narrow-64b_ckd-definition.nc"
    @test basename(paths.shortwave) == "ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc"

    @test_throws ArgumentError official_ecckd_definition_path(:not_a_model; require = false)
end
