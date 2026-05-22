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
    @test haskey(result, "official_training_summary")
    @test result["official_training_summary"]["present"]
    @test result["official_training_summary"]["reactant_status"] == "passed"
    @test result["official_training_summary"]["enzyme_status"] == "passed"
    @test result["official_training_summary"]["final_objective_target_ratio"] > 1
    @test !result["official_training_summary"]["hard_accuracy_target_met"]
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
    @test haskey(result, "reduced_near_miss_summary")
    @test result["reduced_near_miss_summary"]["present"]
    @test result["reduced_near_miss_summary"]["ng_lw"] == 32
    @test result["reduced_near_miss_summary"]["ng_sw"] == 31
    @test result["reduced_near_miss_summary"]["limiting_metric"] == "heating_rate_rmse"
    @test result["reduced_near_miss_summary"]["limiting_metric_ratio"] > 12
    @test result["reduced_near_miss_summary"]["objective_best_omitted_gpoint"] == 23
    @test result["reduced_near_miss_summary"]["objective_best_limiting_metric"] ==
          "heating_rate_rmse"
    @test 1.6 < result["reduced_near_miss_summary"]["objective_best_objective"] < 1.7
    @test haskey(result, "reduced_weight_coordinate_summary")
    @test result["reduced_weight_coordinate_summary"]["present"]
    @test result["reduced_weight_coordinate_summary"]["accepted"]
    @test result["reduced_weight_coordinate_summary"]["omitted_gpoint"] == 23
    @test 1.48 < result["reduced_weight_coordinate_summary"]["accepted_objective"] < 1.49
    @test result["reduced_weight_coordinate_summary"]["accepted_worst_boundary_forcing_error_w_m2"] < 0.3
    @test haskey(result, "reduced_weight_coordinate_descent_summary")
    @test result["reduced_weight_coordinate_descent_summary"]["present"]
    @test result["reduced_weight_coordinate_descent_summary"]["accepted"]
    @test result["reduced_weight_coordinate_descent_summary"]["omitted_gpoint"] == 23
    @test result["reduced_weight_coordinate_descent_summary"]["accepted_move_count"] == 6
    @test 1.26 < result["reduced_weight_coordinate_descent_summary"]["final_objective"] < 1.27
    @test result["reduced_weight_coordinate_descent_summary"]["final_worst_boundary_forcing_error_w_m2"] < 0.3
    @test haskey(result, "reduced_weight_coordinate_descent_continuation_summary")
    @test result["reduced_weight_coordinate_descent_continuation_summary"]["present"]
    @test result["reduced_weight_coordinate_descent_continuation_summary"]["accepted"]
    @test result["reduced_weight_coordinate_descent_continuation_summary"]["omitted_gpoint"] == 23
    @test result["reduced_weight_coordinate_descent_continuation_summary"]["accepted_move_count"] == 8
    @test 1.08 < result["reduced_weight_coordinate_descent_continuation_summary"]["final_objective"] < 1.09
    @test result["reduced_weight_coordinate_descent_continuation_summary"]["final_worst_boundary_forcing_error_w_m2"] < 0.3
    @test haskey(result, "reduced_weight_coordinate_boundary_polish_summary")
    @test result["reduced_weight_coordinate_boundary_polish_summary"]["present"]
    @test result["reduced_weight_coordinate_boundary_polish_summary"]["accepted"]
    @test result["reduced_weight_coordinate_boundary_polish_summary"]["passed"]
    @test result["reduced_weight_coordinate_boundary_polish_summary"]["omitted_gpoint"] == 23
    @test result["reduced_weight_coordinate_boundary_polish_summary"]["accepted_move_count"] == 17
    @test 0.99 < result["reduced_weight_coordinate_boundary_polish_summary"]["final_objective"] < 1.0
    @test result["reduced_weight_coordinate_boundary_polish_summary"]["final_worst_boundary_forcing_error_w_m2"] < 0.3
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
    @test occursin("Reduced near-miss limiter", markdown)
    @test occursin("Objective-best dense leave-one-out diagnostic", markdown)
    @test occursin("omitted SW g-point 23", markdown)
    @test occursin("Exact weight-coordinate improvement", markdown)
    @test occursin("Exact weight-coordinate descent", markdown)
    @test occursin("heating_rate_rmse", markdown)
    @test occursin("Quantitative Training-Recovery Status", markdown)
    @test occursin("final/target", markdown)
    @test occursin("RH_CKDMIP_DATA_PATH", markdown) ||
          occursin("derived ecCKD", markdown)
end
