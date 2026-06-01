using JSON

module ECCKDDerivedFluxGenerationPlanValidation
include(joinpath(@__DIR__, "..", "validation", "ecckd_derived_flux_generation_plan.jl"))
end

@testset "ecCKD derived flux generation plan" begin
    root = normpath(joinpath(@__DIR__, ".."))
    json_path = joinpath(root, "validation", "results", "ecckd_derived_flux_generation_plan.json")
    md_path = joinpath(root, "validation", "results", "ecckd_derived_flux_generation_plan.md")

    redirect_stdout(devnull) do
        ECCKDDerivedFluxGenerationPlanValidation.ecckd_derived_flux_generation_plan_main()
    end

    @test isfile(json_path)
    @test isfile(md_path)

    result = JSON.parsefile(json_path)
    @test result["case"] == "ecckd_derived_flux_generation_plan"
    @test result["expected_derived_flux_count"] == 18
    @test 0 <= result["present_derived_flux_count"] <= 18
    @test 0 <= result["missing_derived_flux_count"] <= 18
    @test result["present_derived_flux_count"] + result["missing_derived_flux_count"] == 18
    @test result["expected_raw_chunk_count"] >= result["present_raw_chunk_count"]
    @test result["expected_raw_chunk_count"] >= 0
    @test haskey(result, "completed_equivalent_raw_chunk_count")
    @test result["expected_raw_chunk_count"] >= result["completed_equivalent_raw_chunk_count"]
    @test result["completed_equivalent_raw_chunk_count"] >= result["present_raw_chunk_count"]
    @test haskey(result, "raw_chunk_rate")
    @test haskey(result["raw_chunk_rate"], "observed_raw_chunk_rate_per_hour")
    @test haskey(result["raw_chunk_rate"], "estimated_hours_remaining")
    @test haskey(result, "ncrcat")
    @test haskey(result["ncrcat"], "present")
    @test haskey(result["ncrcat"], "path")
    @test haskey(result["ncrcat"], "julia_concat_shim")
    @test haskey(result["ncrcat"], "uses_julia_concat")
    @test haskey(result, "derived_flux_progress")
    @test length(result["derived_flux_progress"]) == 18
    @test result["upstream_preflight_status"] in (
        "missing_ckdmip_data_root",
        "incomplete_ckdmip_upstream_data",
        "ready_for_derived_flux_generation",
        "ready_for_original_ecckd_objective",
    )
    @test "5gas-415" in result["lw_scenarios"]
    @test "rel-415" in result["lw_scenarios"]
    @test "rel-415" in result["sw_scenarios"]
    @test !("5gas-415" in result["sw_scenarios"])
    @test any(row -> row["path"] == "test/run_lw_lbl_evaluation.sh" && row["present"],
              result["required_ecckd_scripts"])
    @test any(row -> row["path"] == "test/run_sw_lbl_evaluation.sh" && row["present"],
              result["required_ecckd_scripts"])
    @test all(row -> haskey(row, "raw_chunk_count") &&
                     haskey(row, "expected_raw_chunk_count") &&
                     haskey(row, "missing_raw_chunk_count") &&
                     haskey(row, "raw_chunk_bytes") &&
                     haskey(row, "first_raw_chunk_unix_time") &&
                     haskey(row, "latest_raw_chunk_unix_time") &&
                     haskey(row, "progress_state"),
              result["derived_flux_progress"])

    md = read(md_path, String)
    @test occursin("ecCKD Derived Flux Generation Plan", md)
    @test occursin("Raw Chunk Progress", md)
    @test occursin("Concatenation Tool", md)
    @test occursin("Completed-equivalent raw chunks", md)
    @test occursin("Observed raw chunk rate", md)
    @test occursin("not public CKDMIP archive files", md)
    @test occursin("run_lw_lbl_evaluation.sh", md)
    @test occursin("run_sw_lbl_evaluation.sh", md)
    @test occursin("generate_ecckd_derived_fluxes.sh", md)
end
