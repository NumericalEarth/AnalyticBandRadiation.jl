using JSON

module RecoveryGoalAuditValidation
include(joinpath(@__DIR__, "..", "validation", "recovery_goal_audit.jl"))
end

@testset "recovery goal audit" begin
    root = normpath(joinpath(@__DIR__, ".."))
    json_path = joinpath(root, "validation", "results", "recovery_goal_audit.json")
    md_path = joinpath(root, "validation", "results", "recovery_goal_audit.md")

    redirect_stdout(devnull) do
        RecoveryGoalAuditValidation.recovery_goal_audit_main()
    end

    @test isfile(json_path)
    @test isfile(md_path)

    result = JSON.parsefile(json_path)
    @test result["case"] == "recovery_goal_audit"
    @test result["status"] in ("not_complete", "complete")
    @test result["blocked_count"] >= 0
    @test result["partial_count"] >= 0
    @test haskey(result, "prompt_to_artifact_checklist")
    @test length(result["prompt_to_artifact_checklist"]) == 4
    @test result["unmet_requirement_count"] ==
          count(item -> !item["covered"], result["prompt_to_artifact_checklist"])
    @test Set(result["unmet_requirement_ids"]) ==
          Set(item["requirement_id"] for item in result["prompt_to_artifact_checklist"]
              if !item["covered"])
    @test haskey(result, "original_objective_assets_ready")
    @test result["teacher_student_recovery_status"] == "passed"
    @test haskey(result, "derived_flux_progress")
    @test haskey(result, "derived_flux_plan_status")
    @test haskey(result, "derived_flux_plan_summary")
    @test haskey(result["derived_flux_plan_summary"], "ncrcat")
    @test result["derived_flux_progress"]["expected"] >= 0
    @test result["derived_flux_progress"]["final_present"] >= 0
    @test result["derived_flux_progress"]["products_with_raw_chunks"] >= 0
    @test result["derived_flux_progress"]["expected_raw_chunks"] >=
          result["derived_flux_progress"]["present_raw_chunks"]
    @test haskey(result["derived_flux_progress"], "completed_equivalent_raw_chunks")
    @test result["derived_flux_progress"]["expected_raw_chunks"] >=
          result["derived_flux_progress"]["completed_equivalent_raw_chunks"]
    @test result["derived_flux_progress"]["completed_equivalent_raw_chunks"] >=
          result["derived_flux_progress"]["present_raw_chunks"]
    @test haskey(result["derived_flux_progress"], "observed_raw_chunk_rate_per_hour")
    @test haskey(result["derived_flux_progress"], "estimated_raw_chunk_hours_remaining")
    @test result["objective_reconstruction_status"] in (
        "blocked_missing_original_training_assets",
        "ready_to_reconstruct_original_objective",
        "passed",
    )
    @test result["ckdmip_preflight_status"] in (
        "missing_ckdmip_data_root",
        "incomplete_ckdmip_upstream_data",
        "ready_for_derived_flux_generation",
        "ready_for_original_ecckd_objective",
    )
    @test result["reduced_model_summary"]["full_32x32_passed"]
    @test result["reduced_model_summary"]["hard_boundary_forcing_threshold_w_m2"] == 0.3
    @test result["reduced_model_summary"]["official_32x32_worst_boundary_forcing_error_w_m2"] < 0.3
    @test haskey(result["reduced_model_summary"], "best_reduced_candidate")
    @test result["breeze_summary"]["runtime_supported"]
    @test result["breeze_summary"]["final_4x_claim_supported"]

    mktempdir() do dir
        missing_breeze = joinpath(dir, "missing_breeze.json")
        env_result = withenv("RH_BREEZE_RCEMIP_JSON" => missing_breeze) do
            RecoveryGoalAuditValidation.run_recovery_goal_audit()
        end
        @test env_result.breeze_summary.artifact == missing_breeze
        @test !env_result.breeze_summary.present
        @test any(row -> row.id == "breeze_dynamic_integration" && row.status == "blocked",
                  env_result.requirements)
    end

    statuses = Dict(row["id"] => row["status"] for row in result["requirements"])
    @test statuses["ecrad_full_and_reduced_parity"] in ("blocked", "partial", "passed")
    @test statuses["reduced_vs_rrtmgp_representative_states"] in ("blocked", "partial", "passed")
    @test statuses["breeze_dynamic_integration"] == "passed"
    @test statuses["reactant_enzyme_ecckd_training_recovery"] in ("blocked", "partial", "passed")

    markdown = read(md_path, String)
    @test occursin("Recovery Goal Audit", markdown)
    @test occursin("Prompt-to-Artifact Checklist", markdown)
    @test occursin("Derived flux products:", markdown)
    @test occursin("Derived flux generation plan status", markdown)
    @test occursin("Derived raw chunks:", markdown)
    @test occursin("Completed-equivalent derived raw chunks", markdown)
    @test occursin("ncrcat", markdown)
    @test occursin("Observed derived raw chunk rate", markdown)
    @test occursin("Quantitative Reduced-Model Status", markdown)
    @test occursin("RH_CKDMIP_DATA_PATH", markdown) ||
          occursin("derived ecCKD", markdown)
end
