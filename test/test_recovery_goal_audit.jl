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
    @test haskey(result, "original_objective_terms_summary")
    @test result["original_objective_terms_summary"]["present"]
    @test result["original_objective_terms_summary"]["status"] ==
          "objective_terms_captured"
    @test result["original_objective_terms_summary"]["implementation_status"] ==
          "terms_captured_not_yet_recovered"
    @test result["original_objective_terms_summary"]["longwave_term_count"] >= 8
    @test result["original_objective_terms_summary"]["shortwave_term_count"] >= 10
    @test result["original_objective_terms_summary"]["longwave_terms_present"]
    @test result["original_objective_terms_summary"]["shortwave_terms_present"]
    @test haskey(result, "ckdmip_objective_dataset_summary")
    @test result["ckdmip_objective_dataset_summary"]["present"]
    @test result["ckdmip_objective_dataset_summary"]["status"] ==
          "dataset_samples_ready"
    @test result["ckdmip_objective_dataset_summary"]["sample_count"] == 2
    @test result["ckdmip_objective_dataset_summary"]["training_flux_file_count"] >= 52
    @test result["ckdmip_objective_dataset_summary"]["training_flux_schema_ok_count"] ==
          result["ckdmip_objective_dataset_summary"]["training_flux_file_count"]
    @test result["ckdmip_objective_dataset_summary"]["longwave_sample_ready"]
    @test result["ckdmip_objective_dataset_summary"]["shortwave_sample_ready"]
    @test result["ckdmip_objective_dataset_summary"]["self_loss_zero"]
    @test haskey(result, "ckdmip_objective_ad_batch_summary")
    @test result["ckdmip_objective_ad_batch_summary"]["present"]
    @test result["ckdmip_objective_ad_batch_summary"]["status"] ==
          "optimizer_batch_ready"
    @test result["ckdmip_objective_ad_batch_summary"]["parameter_count"] > 0
    @test result["ckdmip_objective_ad_batch_summary"]["accepted_step"]
    @test result["ckdmip_objective_ad_batch_summary"]["loss_reduction_factor"] > 1
    @test result["ckdmip_objective_ad_batch_summary"]["gradient_method"] ==
          "central_finite_difference"
    @test haskey(result, "published_recovery_target_summary")
    @test result["published_recovery_target_summary"]["present"]
    @test result["published_recovery_target_summary"]["status"] ==
          "published_recovery_target_ready"
    @test result["published_recovery_target_summary"]["model_count"] == 6
    @test result["published_recovery_target_summary"]["shortwave_32_parameter_count"] > 0
    @test result["published_recovery_target_summary"]["longwave_32_parameter_count"] > 0
    @test result["published_recovery_target_summary"]["final_objective_target_ratio_max"] <= 1.05
    @test result["published_recovery_target_summary"]["optical_depth_log_rmse_max"] <= 0.02
    @test result["published_recovery_target_summary"]["optimizer_only_delta_rule_present"]
    @test haskey(result, "published_recovery_vector_summary")
    @test result["published_recovery_vector_summary"]["present"]
    @test result["published_recovery_vector_summary"]["status"] == "passed"
    @test result["published_recovery_vector_summary"]["array_count"] == 9
    @test result["published_recovery_vector_summary"]["parameter_count"] == 204896
    @test result["published_recovery_vector_summary"]["roundtrip_max_abs_error"] == 0.0
    @test result["published_recovery_vector_summary"]["roundtrip_l1_relative_error"] == 0.0
    @test result["published_recovery_vector_summary"]["recovery_metrics_status"] == "passed"
    @test haskey(result, "published_recovery_vector_training_summary")
    @test result["published_recovery_vector_training_summary"]["present"]
    @test result["published_recovery_vector_training_summary"]["status"] == "passed"
    @test result["published_recovery_vector_training_summary"]["parameter_count"] == 204896
    @test result["published_recovery_vector_training_summary"]["trained_parameter_count"] == 64
    @test result["published_recovery_vector_training_summary"]["loss_reduction_factor"] > 1.0e10
    @test !result["published_recovery_vector_training_summary"]["enzyme_requested"]
    @test !result["published_recovery_vector_training_summary"]["reactant_check_requested"]
    @test result["published_recovery_vector_training_summary"]["recovery_metrics_status"] == "passed"
    @test haskey(result, "candidate_objective_score_summary")
    @test result["candidate_objective_score_summary"]["present"]
    @test result["candidate_objective_score_summary"]["status"] ==
          "candidate_objective_score_ready"
    @test result["candidate_objective_score_summary"]["candidate_metrics_status"] == "passed"
    @test result["candidate_objective_score_summary"]["sample_count"] == 2
    @test result["candidate_objective_score_summary"]["zero_forward_losses_zero"]
    @test result["candidate_objective_score_summary"]["perturbation_increases_loss"]
    @test haskey(result, "candidate_transfer_smoke_summary")
    @test result["candidate_transfer_smoke_summary"]["present"]
    @test result["candidate_transfer_smoke_summary"]["status"] ==
          "candidate_transfer_smoke_ready"
    @test result["candidate_transfer_smoke_summary"]["layer_count"] == 54
    @test result["candidate_transfer_smoke_summary"]["longwave_gpoints"] == 32
    @test result["candidate_transfer_smoke_summary"]["shortwave_gpoints"] == 32
    @test result["candidate_transfer_smoke_summary"]["longwave_loss"] > 0
    @test result["candidate_transfer_smoke_summary"]["shortwave_loss"] > 0
    @test result["candidate_transfer_smoke_summary"]["longwave_projection_band_count"] > 0
    @test result["candidate_transfer_smoke_summary"]["shortwave_projection_band_count"] > 0
    @test result["candidate_transfer_smoke_summary"]["longwave_projection_loss"] > 0
    @test result["candidate_transfer_smoke_summary"]["shortwave_projection_loss"] > 0
    @test haskey(result, "candidate_transfer_optimizer_probe_summary")
    @test result["candidate_transfer_optimizer_probe_summary"]["present"]
    @test result["candidate_transfer_optimizer_probe_summary"]["status"] ==
          "optimizer_probe_passed"
    @test result["candidate_transfer_optimizer_probe_summary"]["parameter_count"] == 4
    @test result["candidate_transfer_optimizer_probe_summary"]["accepted_step"]
    @test result["candidate_transfer_optimizer_probe_summary"]["final_loss"] <
          result["candidate_transfer_optimizer_probe_summary"]["initial_loss"]
    @test result["candidate_transfer_optimizer_probe_summary"]["loss_reduction_factor"] > 1
    @test haskey(result, "candidate_table_parameter_probe_summary")
    @test result["candidate_table_parameter_probe_summary"]["present"]
    @test result["candidate_table_parameter_probe_summary"]["status"] ==
          "table_parameter_probe_passed"
    @test result["candidate_table_parameter_probe_summary"]["parameter_count"] == 4
    @test result["candidate_table_parameter_probe_summary"]["accepted_step"]
    @test result["candidate_table_parameter_probe_summary"]["final_loss"] <
          result["candidate_table_parameter_probe_summary"]["initial_loss"]
    @test result["candidate_table_parameter_probe_summary"]["loss_reduction_factor"] > 1
    @test haskey(result, "candidate_table_writeback_probe_summary")
    @test result["candidate_table_writeback_probe_summary"]["present"]
    @test result["candidate_table_writeback_probe_summary"]["status"] ==
          "table_writeback_probe_passed"
    @test result["candidate_table_writeback_probe_summary"]["accepted_writeback"]
    @test result["candidate_table_writeback_probe_summary"]["written_projected_loss"] <
          result["candidate_table_writeback_probe_summary"]["baseline_projected_loss"]
    @test result["candidate_table_writeback_probe_summary"]["loss_reduction_factor"] > 1
    @test haskey(result, "candidate_table_writeback_continuation_summary")
    @test result["candidate_table_writeback_continuation_summary"]["present"]
    @test result["candidate_table_writeback_continuation_summary"]["status"] ==
          "table_writeback_continuation_passed"
    @test result["candidate_table_writeback_continuation_summary"]["accepted_writeback"]
    @test result["candidate_table_writeback_continuation_summary"]["accepted_steps"] >= 1
    @test result["candidate_table_writeback_continuation_summary"]["final_in_memory_loss"] <=
          result["candidate_table_writeback_continuation_summary"]["initial_in_memory_loss"]
    @test result["candidate_table_writeback_continuation_summary"]["written_projected_loss"] <
          result["candidate_table_writeback_continuation_summary"]["baseline_projected_loss"]
    @test result["candidate_table_writeback_continuation_summary"]["writeback_loss_reduction_factor"] > 1
    @test haskey(result, "candidate_table_writeback_multisample_summary")
    @test result["candidate_table_writeback_multisample_summary"]["present"]
    @test result["candidate_table_writeback_multisample_summary"]["status"] in
          ("table_writeback_multisample_passed", "table_writeback_multisample_partial")
    @test result["candidate_table_writeback_multisample_summary"]["scenario_count"] >= 1
    @test result["candidate_table_writeback_multisample_summary"]["present_count"] >= 1
    @test result["candidate_table_writeback_multisample_summary"]["improved_count"] >= 0
    @test result["candidate_table_writeback_multisample_summary"]["worst_loss_ratio"] >= 0
    @test haskey(result, "candidate_table_multisample_optimizer_summary")
    @test result["candidate_table_multisample_optimizer_summary"]["present"]
    @test result["candidate_table_multisample_optimizer_summary"]["status"] in (
        "table_multisample_optimizer_passed",
        "table_multisample_optimizer_seed_passed",
        "table_multisample_optimizer_partial",
        "table_multisample_optimizer_failed",
    )
    @test result["candidate_table_multisample_optimizer_summary"]["present_count"] >= 1
    @test result["candidate_table_multisample_optimizer_summary"]["accepted_steps"] >= 0
    @test result["candidate_table_multisample_optimizer_summary"]["writeback_improved_count"] >= 0
    @test result["candidate_table_multisample_optimizer_summary"]["writeback_worst_loss_ratio"] >= 0
    @test haskey(result, "candidate_table_written_coordinate_scan_summary")
    @test result["candidate_table_written_coordinate_scan_summary"]["present"]
    @test result["candidate_table_written_coordinate_scan_summary"]["status"] in (
        "written_coordinate_scan_improved",
        "written_coordinate_scan_no_descent",
        "written_coordinate_scan_failed",
    )
    @test result["candidate_table_written_coordinate_scan_summary"]["present_count"] >= 1
    @test result["candidate_table_written_coordinate_scan_summary"]["tested_move_count"] >= 1
    @test result["candidate_table_written_coordinate_scan_summary"]["aggregate_loss_reduction_factor"] >= 0
    @test result["candidate_table_written_coordinate_scan_summary"]["best_improved_count"] >= 0
    @test result["candidate_table_written_coordinate_scan_summary"]["best_worst_loss_ratio"] >= 0
    @test haskey(result, "candidate_table_written_coordinate_descent_summary")
    @test result["candidate_table_written_coordinate_descent_summary"]["present"]
    @test result["candidate_table_written_coordinate_descent_summary"]["status"] in (
        "written_coordinate_descent_improved",
        "written_coordinate_descent_no_descent",
        "written_coordinate_descent_failed",
    )
    @test result["candidate_table_written_coordinate_descent_summary"]["present_count"] >= 1
    @test result["candidate_table_written_coordinate_descent_summary"]["accepted_move_count"] >= 0
    @test result["candidate_table_written_coordinate_descent_summary"]["aggregate_loss_reduction_factor"] >= 0
    @test result["candidate_table_written_coordinate_descent_summary"]["final_improved_count"] >= 0
    @test result["candidate_table_written_coordinate_descent_summary"]["final_worst_loss_ratio"] >= 0
    @test haskey(result, "candidate_table_written_minimax_descent_summary")
    @test result["candidate_table_written_minimax_descent_summary"]["present"]
    @test result["candidate_table_written_minimax_descent_summary"]["status"] in (
        "written_minimax_descent_improved",
        "written_minimax_descent_no_descent",
        "written_minimax_descent_failed",
    )
    @test result["candidate_table_written_minimax_descent_summary"]["present_count"] >= 1
    @test result["candidate_table_written_minimax_descent_summary"]["tested_move_count"] >= 1
    @test result["candidate_table_written_minimax_descent_summary"]["accepted_move_count"] >= 0
    @test result["candidate_table_written_minimax_descent_summary"]["aggregate_loss_reduction_factor"] >= 0
    @test result["candidate_table_written_minimax_descent_summary"]["worst_loss_ratio_reduction_factor"] >= 0
    @test result["candidate_table_written_minimax_descent_summary"]["final_improved_count"] >= 0
    @test result["candidate_table_written_minimax_descent_summary"]["final_worst_loss_ratio"] >= 0
    @test haskey(result, "candidate_table_humidity_split_probe_summary")
    @test result["candidate_table_humidity_split_probe_summary"]["present"]
    @test result["candidate_table_humidity_split_probe_summary"]["status"] in (
        "humidity_split_probe_improved",
        "humidity_split_probe_no_descent",
    )
    @test result["candidate_table_humidity_split_probe_summary"]["present_count"] >= 1
    @test result["candidate_table_humidity_split_probe_summary"]["tested_move_count"] >= 1
    @test result["candidate_table_humidity_split_probe_summary"]["aggregate_loss_reduction_factor"] >= 0
    @test result["candidate_table_humidity_split_probe_summary"]["worst_loss_ratio_reduction_factor"] >= 0
    @test result["candidate_table_humidity_split_probe_summary"]["best_aggregate_worst_loss_ratio"] >= 0
    @test result["candidate_table_humidity_split_probe_summary"]["best_minimax_worst_loss_ratio"] >= 0
    @test haskey(result, "candidate_table_humidity_split_descent_summary")
    @test result["candidate_table_humidity_split_descent_summary"]["present"]
    @test result["candidate_table_humidity_split_descent_summary"]["status"] in (
        "humidity_split_descent_improved",
        "humidity_split_descent_no_descent",
        "humidity_split_descent_failed",
    )
    @test result["candidate_table_humidity_split_descent_summary"]["present_count"] >= 1
    @test result["candidate_table_humidity_split_descent_summary"]["tested_move_count"] >= 1
    @test result["candidate_table_humidity_split_descent_summary"]["accepted_move_count"] >= 0
    @test result["candidate_table_humidity_split_descent_summary"]["aggregate_loss_reduction_factor"] >= 0
    @test result["candidate_table_humidity_split_descent_summary"]["worst_loss_ratio_reduction_factor"] >= 0
    @test result["candidate_table_humidity_split_descent_summary"]["final_improved_count"] >= 0
    @test result["candidate_table_humidity_split_descent_summary"]["final_worst_loss_ratio"] >= 0
    @test haskey(result, "candidate_table_humidity_tribin_probe_summary")
    @test result["candidate_table_humidity_tribin_probe_summary"]["present"]
    @test result["candidate_table_humidity_tribin_probe_summary"]["status"] in (
        "humidity_tribin_probe_improved",
        "humidity_tribin_probe_no_descent",
    )
    @test result["candidate_table_humidity_tribin_probe_summary"]["present_count"] >= 1
    @test result["candidate_table_humidity_tribin_probe_summary"]["tested_move_count"] >= 1
    @test result["candidate_table_humidity_tribin_probe_summary"]["aggregate_loss_reduction_factor"] >= 0
    @test result["candidate_table_humidity_tribin_probe_summary"]["worst_loss_ratio_reduction_factor"] >= 0
    @test result["candidate_table_humidity_tribin_probe_summary"]["best_aggregate_worst_loss_ratio"] >= 0
    @test result["candidate_table_humidity_tribin_probe_summary"]["best_minimax_worst_loss_ratio"] >= 0
    @test haskey(result, "candidate_table_humidity_tribin_descent_summary")
    @test result["candidate_table_humidity_tribin_descent_summary"]["present"]
    @test result["candidate_table_humidity_tribin_descent_summary"]["status"] in (
        "humidity_tribin_descent_improved",
        "humidity_tribin_descent_no_descent",
        "humidity_tribin_descent_failed",
    )
    @test result["candidate_table_humidity_tribin_descent_summary"]["present_count"] >= 1
    @test result["candidate_table_humidity_tribin_descent_summary"]["tested_move_count"] >= 1
    @test result["candidate_table_humidity_tribin_descent_summary"]["accepted_move_count"] >= 0
    @test result["candidate_table_humidity_tribin_descent_summary"]["aggregate_loss_reduction_factor"] >= 0
    @test result["candidate_table_humidity_tribin_descent_summary"]["worst_loss_ratio_reduction_factor"] >= 0
    @test result["candidate_table_humidity_tribin_descent_summary"]["final_improved_count"] >= 0
    @test result["candidate_table_humidity_tribin_descent_summary"]["final_worst_loss_ratio"] >= 0
    @test haskey(result, "candidate_table_humidity_tribin_weighted_descent_summary")
    @test result["candidate_table_humidity_tribin_weighted_descent_summary"]["present"]
    @test result["candidate_table_humidity_tribin_weighted_descent_summary"]["status"] in (
        "humidity_tribin_weighted_descent_improved",
        "humidity_tribin_weighted_descent_no_descent",
    )
    @test result["candidate_table_humidity_tribin_weighted_descent_summary"]["present_count"] >= 1
    @test result["candidate_table_humidity_tribin_weighted_descent_summary"]["tested_weight_set_count"] >= 1
    @test result["candidate_table_humidity_tribin_weighted_descent_summary"]["best_accepted_move_count"] >= 0
    @test result["candidate_table_humidity_tribin_weighted_descent_summary"]["best_weighted_loss_ratio"] >= 0
    @test result["candidate_table_humidity_tribin_weighted_descent_summary"]["best_weighted_loss_ratio_reduction_factor"] >= 0
    @test result["candidate_table_humidity_tribin_weighted_descent_summary"]["best_aggregate_loss_reduction_factor"] >= 0
    @test result["candidate_table_humidity_tribin_weighted_descent_summary"]["best_worst_loss_ratio"] >= 0
    @test result["candidate_table_humidity_tribin_weighted_descent_summary"]["best_worst_loss_ratio_reduction_factor"] >= 0
    @test result["candidate_table_humidity_tribin_weighted_descent_summary"]["best_final_improved_count"] >= 0
    @test haskey(result, "candidate_table_humidity_tribin_constrained_descent_summary")
    @test result["candidate_table_humidity_tribin_constrained_descent_summary"]["present"]
    @test result["candidate_table_humidity_tribin_constrained_descent_summary"]["status"] in (
        "humidity_tribin_constrained_descent_improved",
        "humidity_tribin_constrained_descent_no_descent",
    )
    @test result["candidate_table_humidity_tribin_constrained_descent_summary"]["present_count"] >= 1
    @test result["candidate_table_humidity_tribin_constrained_descent_summary"]["tested_tolerance_count"] >= 1
    @test result["candidate_table_humidity_tribin_constrained_descent_summary"]["best_accepted_move_count"] >= 0
    @test result["candidate_table_humidity_tribin_constrained_descent_summary"]["best_aggregate_loss_reduction_factor"] >= 0
    @test result["candidate_table_humidity_tribin_constrained_descent_summary"]["best_worst_loss_ratio"] >= 0
    @test result["candidate_table_humidity_tribin_constrained_descent_summary"]["best_worst_loss_ratio_reduction_factor"] >= 0
    @test result["candidate_table_humidity_tribin_constrained_descent_summary"]["best_final_improved_count"] >= 0
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
    @test haskey(result, "published_model_accuracy_summary")
    @test result["published_model_accuracy_summary"]["present"]
    @test result["published_model_accuracy_summary"]["status"] == "passed"
    @test result["published_model_accuracy_summary"]["model_count"] == 6
    @test result["published_model_accuracy_summary"]["passed_count"] == 6
    @test result["published_model_accuracy_summary"]["boundary_compatible_count"] == 6
    @test result["published_model_accuracy_summary"]["isolation_diagnostic_count"] == 3
    @test result["published_model_accuracy_summary"]["boundary_projection_diagnostic_count"] == 5
    @test haskey(result, "matched_reference_summary")
    @test result["matched_reference_summary"]["present"]
    @test result["matched_reference_summary"]["status"] ==
          "ready_for_published_parity_validation"
    @test result["matched_reference_summary"]["required_case_count"] == 16
    @test result["matched_reference_summary"]["missing_case_count"] == 0
    @test result["matched_reference_summary"]["all_sky_case_count"] == 6
    @test result["matched_reference_summary"]["all_sky_ready_count"] == 6
    @test haskey(result, "published_all_sky_accuracy_summary")
    @test result["published_all_sky_accuracy_summary"]["present"]
    @test result["published_all_sky_accuracy_summary"]["status"] in
          ("passed", "failed_threshold")
    @test result["published_all_sky_accuracy_summary"]["model_count"] == 6
    @test result["published_all_sky_accuracy_summary"]["passed_count"] >= 1
    @test result["published_all_sky_accuracy_summary"]["passed_count"] <= 6
    @test result["published_all_sky_accuracy_summary"]["worst_hard_objective"] > 0
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
    @test occursin("64x96", markdown) || occursin("64x64", markdown)
    @test occursin("boundary-projection diagnostics", markdown)
    @test occursin("Reduced near-miss limiter", markdown)
    @test occursin("Objective-best dense leave-one-out diagnostic", markdown)
    @test occursin("omitted SW g-point 23", markdown)
    @test occursin("Exact weight-coordinate improvement", markdown)
    @test occursin("Exact weight-coordinate descent", markdown)
    @test occursin("heating_rate_rmse", markdown)
    @test occursin("Quantitative Training-Recovery Status", markdown)
    @test occursin("Original objective terms:", markdown)
    @test occursin("CKDMIP objective dataset:", markdown)
    @test occursin("CKDMIP objective optimizer batch:", markdown)
    @test occursin("Published recovery target:", markdown)
    @test occursin("Published recovery vector:", markdown)
    @test occursin("Published recovery vector training:", markdown)
    @test occursin("Candidate original-objective score:", markdown)
    @test occursin("Candidate transfer smoke:", markdown)
    @test occursin("Candidate transfer optimizer probe:", markdown)
    @test occursin("Candidate table-parameter probe:", markdown)
    @test occursin("Candidate table-writeback probe:", markdown)
    @test occursin("Candidate table-writeback continuation:", markdown)
    @test occursin("Candidate table-writeback multisample:", markdown)
    @test occursin("Candidate table multisample optimizer:", markdown)
    @test occursin("Candidate table written-coordinate scan:", markdown)
    @test occursin("Candidate table written-coordinate descent:", markdown)
    @test occursin("Candidate table written minimax descent:", markdown)
    @test occursin("Candidate table humidity-split probe:", markdown)
    @test occursin("Candidate table humidity-split descent:", markdown)
    @test occursin("Candidate table humidity-tribin probe:", markdown)
    @test occursin("Candidate table humidity-tribin descent:", markdown)
    @test occursin("Candidate table humidity-tribin weighted descent:",
                   markdown)
    @test occursin("Candidate table humidity-tribin constrained descent:",
                   markdown)
    @test occursin("projected loss", markdown)
    @test occursin("final/target", markdown)
    @test occursin("RH_CKDMIP_DATA_PATH", markdown) ||
          occursin("derived ecCKD", markdown)
end
