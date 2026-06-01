using JSON

module EcckdModelInventoryValidation
include(joinpath(@__DIR__, "..", "validation", "ecckd_model_inventory.jl"))
end

@testset "official ecCKD model inventory artifact" begin
    root = normpath(joinpath(@__DIR__, ".."))
    json_path = joinpath(root, "validation", "results", "ecckd_model_inventory.json")
    md_path = joinpath(root, "validation", "results", "ecckd_model_inventory.md")
    redirect_stdout(devnull) do
        EcckdModelInventoryValidation.ecckd_model_inventory_main()
    end

    @test isfile(json_path)
    @test isfile(md_path)
    output = read(md_path, String)
    @test occursin("ecCKD Model Inventory", output)
    @test occursin("Status: **passed**", output)
    @test occursin("ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc", output)
    @test occursin("ecckd-1.4_sw_climate_vfine-96b_ckd-definition.nc", output)

    result = JSON.parsefile(json_path)
    @test result["case"] == "ecckd_model_inventory"
    @test result["status"] == "passed"
    @test length(result["entries"]) == 6

    entries = Dict(entry["filename"] => entry for entry in result["entries"])
    @test sort([entry["gpoints"] for entry in result["entries"]]) == [32, 32, 32, 64, 64, 96]
    @test entries["ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc"]["source_tables_present"]
    @test !entries["ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc"]["rayleigh_tables_present"]
    @test entries["ecckd-1.2_lw_climate_narrow-64b_ckd-definition.nc"]["source_tables_present"]
    @test !entries["ecckd-1.2_lw_climate_narrow-64b_ckd-definition.nc"]["rayleigh_tables_present"]
    for filename in (
        "ecckd-1.0_sw_climate_rgb-32b_ckd-definition.nc",
        "ecckd-1.2_sw_climate_window-64b_ckd-definition.nc",
        "ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc",
        "ecckd-1.4_sw_climate_vfine-96b_ckd-definition.nc",
    )
        @test !entries[filename]["source_tables_present"]
        @test entries[filename]["rayleigh_tables_present"]
    end
end
