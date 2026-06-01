using Test

include(joinpath(@__DIR__, "..", "validation", "reduced_ecckd_optimization_preflight.jl"))

optional_ad_skipped() = get(ENV, "RH_SKIP_OPTIONAL_AD_CHECKS", "false") == "true"

@testset "greedy checkpoint iteration budget" begin
    mktemp() do path, io
        close(io)
        parameters = initial_parameters()
        checkpoint_parameters = join(fill("0.0", length(parameters)), ", ")
        write(path, """
        {
          "greedy_coordinate_descent": {
            "iterations_completed": 5,
            "final_objective": 2.0,
            "final_parameters": [$checkpoint_parameters]
          }
        }
        """)
        f(x) = sum(abs2, x)
        result, used_checkpoint = greedy_coordinate_descent_with_checkpoint(
            f,
            parameters,
            f(parameters);
            total_iterations = 2,
            checkpoint_path = path,
        )
        @test used_checkpoint
        @test result.iterations_requested == 5
        @test result.iterations_completed == 5
        @test result.final_objective == 2.0
    end
end

@testset "reduced ecCKD optimization preflight artifact" begin
    result = reduced_optimization_preflight()
    @test result.case == "reduced_ecckd_optimization_preflight"
    @test result.status == "preflight_ready"
    @test result.ng_sw == 16
    @test result.parameter_count == 48
    @test result.objective_target == 1.0
    @test result.initial_objective > result.objective_target
    @test isfinite(result.directional_derivative)
    @test isfinite(result.best_trial_objective)
    @test result.best_trial_objective <= result.initial_objective
    @test length(result.trial_line_search) == 5
    @test result.optimization_smoke.iterations == 1
    @test result.optimization_smoke.improved
    @test result.optimization_smoke.final_objective == result.best_trial_objective
    @test result.optimization_smoke.objective_reduction > 0
    @test result.directional_search.iterations == 2
    @test result.directional_search.improved
    @test result.directional_search.final_objective <= result.optimization_smoke.final_objective
    @test result.directional_search.objective_reduction >= result.optimization_smoke.objective_reduction
    @test length(result.block_trial_scan.rows) == 3
    @test result.block_trial_scan.best_block in ("weights", "absorption", "rayleigh")
    @test isfinite(result.block_trial_scan.best_objective)
    @test result.coordinate_coefficient_scan.candidate_count == 64
    @test result.coordinate_coefficient_scan.best_block in ("absorption", "rayleigh")
    @test isfinite(result.coordinate_coefficient_scan.best_objective)
    @test result.coordinate_coefficient_scan.best_objective < result.initial_objective
    @test result.greedy_coordinate_descent.iterations_completed >= 1
    @test result.greedy_coordinate_descent.improved
    @test result.greedy_coordinate_descent.final_objective <= result.coordinate_coefficient_scan.best_objective
    @test result.greedy_coordinate_descent.objective_reduction > 0
    @test result.coefficient_joint_direction_scan.candidate_count >= 1
    @test isfinite(result.coefficient_joint_direction_scan.best_objective)
    @test result.coefficient_joint_direction_scan.best_objective <=
          result.greedy_coordinate_descent.final_objective ||
          result.coefficient_joint_direction_scan.best_improvement <= 0
    @test result.post_joint_coordinate_refinement.iterations_completed >= 1
    @test isfinite(result.post_joint_coordinate_refinement.final_objective)
    @test result.post_joint_coordinate_refinement.final_objective <=
          min(result.greedy_coordinate_descent.final_objective,
              result.coefficient_joint_direction_scan.best_objective)
    @test result.post_coefficient_weight_refinement.iterations_completed >= 1
    @test isfinite(result.post_coefficient_weight_refinement.final_objective)
    @test result.post_coefficient_weight_refinement.final_objective <=
          result.post_coefficient_weight_refinement.initial_objective
    @test result.finite_difference_coefficient_direction_refinement.active_parameter_count == 32
    @test result.finite_difference_coefficient_direction_refinement.candidate_count >= 32
    @test isfinite(result.finite_difference_coefficient_direction_refinement.gradient_norm)
    @test isfinite(result.finite_difference_coefficient_direction_refinement.best_objective)
    @test result.finite_difference_coefficient_direction_refinement.best_objective <=
          result.finite_difference_coefficient_direction_refinement.initial_objective ||
          !result.finite_difference_coefficient_direction_refinement.improved
    @test result.smooth_objective_coefficient_refinement.beta > 0
    @test result.smooth_objective_coefficient_refinement.active_parameter_count == 32
    @test result.smooth_objective_coefficient_refinement.candidate_count >= 32
    @test isfinite(result.smooth_objective_coefficient_refinement.gradient_norm)
    @test result.smooth_objective_coefficient_refinement.best_smooth_objective <=
          result.smooth_objective_coefficient_refinement.initial_smooth_objective ||
          !result.smooth_objective_coefficient_refinement.smooth_improved
    @test result.smooth_objective_coefficient_refinement.final_objective <=
          result.smooth_objective_coefficient_refinement.initial_hard_objective
    @test result.final_objective_breakdown.objective >= result.objective_target
    @test result.final_objective_breakdown.worst_case in String.(REDUCED_CASE_NAMES)
    @test !isempty(result.final_objective_breakdown.worst_metric)
    @test result.final_objective_breakdown.worst_value > 0
    @test result.final_objective_breakdown.worst_threshold > 0
    @test length(result.final_objective_breakdown.rows) ==
          8 * length(REDUCED_CASES)
    @test result.targeted_worst_metric_refinement.target_case in String.(REDUCED_CASE_NAMES)
    @test !isempty(result.targeted_worst_metric_refinement.target_metric)
    @test result.targeted_worst_metric_refinement.scan.candidate_count >= 1
    @test result.targeted_worst_metric_refinement.target_best_objective <=
          result.targeted_worst_metric_refinement.target_initial_objective ||
          !result.targeted_worst_metric_refinement.target_improved
    @test result.targeted_worst_metric_refinement.final_objective <=
          result.targeted_worst_metric_refinement.initial_full_objective
    @test result.separated_component_refinement.parameter_count == 64
    @test result.separated_component_refinement.target_case in String.(REDUCED_CASE_NAMES)
    @test result.separated_component_refinement.target_metric ==
          result.final_objective_breakdown.worst_metric
    @test result.separated_component_refinement.scan.candidate_count >= 1
    @test result.separated_component_refinement.final_objective <=
          result.separated_component_refinement.initial_full_objective
    @test result.pressure_band_table_refinement.iterations_completed >= 1
    @test result.pressure_band_table_refinement.max_candidates_per_iteration >= 1
    @test result.pressure_band_table_refinement.final_objective <=
          result.pressure_band_table_refinement.initial_full_objective
    @test result.pressure_band_table_refinement.accepted_move_count ==
          length(result.pressure_band_table_refinement.accepted_moves)
    @test isapprox(result.pressure_band_table_refinement.refined_breakdown.objective,
                   result.pressure_band_table_refinement.final_objective;
                   rtol = 1.0e-12)
    @test result.active_table_entry_refinement.iterations_completed >= 0
    @test result.active_table_entry_refinement.max_candidates_per_iteration >= 1
    @test result.active_table_entry_refinement.targeted_candidates_enabled isa Bool
    @test result.active_table_entry_refinement.targeted_candidate_count >=
          result.active_table_entry_refinement.accepted_move_count
    @test result.active_table_entry_refinement.final_objective <=
          result.active_table_entry_refinement.initial_full_objective
    @test result.active_table_entry_refinement.accepted_move_count ==
          length(result.active_table_entry_refinement.accepted_moves)
    @test result.constrained_table_optimizer_target.target_case ==
          result.final_objective_breakdown.worst_case
    @test result.constrained_table_optimizer_target.target_metric ==
          result.final_objective_breakdown.worst_metric
    @test result.constrained_table_optimizer_target.required_absolute_reduction > 0
    @test result.constrained_table_optimizer_target.required_relative_reduction > 0
    @test occursin("nonnegative shortwave coefficient-table",
                   result.constrained_table_optimizer_target.recommended_next_parameterization)
    @test result.warm_started_topology_neighbor_refinement.candidate_count >= 1
    @test result.warm_started_topology_neighbor_refinement.radius >= 1
    @test result.warm_started_topology_neighbor_refinement.total_candidate_count >=
          result.warm_started_topology_neighbor_refinement.candidate_count
    @test isfinite(result.warm_started_topology_neighbor_refinement.best_objective)
    @test result.warm_started_topology_neighbor_refinement.best_objective <=
          result.warm_started_topology_neighbor_refinement.base_objective ||
          !result.warm_started_topology_neighbor_refinement.improved
    @test result.ranked_topology_neighbor_refinement.candidate_count >= 1
    @test result.ranked_topology_neighbor_refinement.radius >= 1
    @test result.ranked_topology_neighbor_refinement.total_candidate_count >=
          result.ranked_topology_neighbor_refinement.candidate_count
    @test result.ranked_topology_neighbor_refinement.ranked_candidate_count ==
          result.ranked_topology_neighbor_refinement.total_candidate_count
    @test isfinite(result.ranked_topology_neighbor_refinement.best_initial_objective)
    @test isfinite(result.ranked_topology_neighbor_refinement.best_objective)
    @test result.ranked_topology_neighbor_refinement.best_objective <=
          result.ranked_topology_neighbor_refinement.base_objective ||
          !result.ranked_topology_neighbor_refinement.improved
    @test result.topology_neighbor_scan.candidate_count >= 1
    @test result.topology_neighbor_scan.total_candidate_count >= result.topology_neighbor_scan.candidate_count
    @test isapprox(result.topology_neighbor_scan.base_objective,
                   min(result.greedy_coordinate_descent.final_objective,
                       result.coefficient_joint_direction_scan.best_objective,
                       result.post_joint_coordinate_refinement.final_objective,
                       result.post_coefficient_weight_refinement.final_objective,
                       result.targeted_worst_metric_refinement.final_objective,
                       result.separated_component_refinement.final_objective,
                       result.warm_started_topology_neighbor_refinement.best_objective,
                       result.ranked_topology_neighbor_refinement.best_objective);
                   rtol = 1.0e-10)
    @test isfinite(result.topology_neighbor_scan.best_objective)
    @test result.topology_neighbor_scan.best_objective <= result.topology_neighbor_scan.base_objective ||
          !result.topology_neighbor_scan.improved
    @test result.topology_candidate_scan.status == "all_32x16_topologies_fail_forcing_gate"
    @test result.topology_candidate_scan.candidate_count >= 1
    @test result.topology_candidate_scan.best_forcing_objective_lower_bound > result.objective_target
    @test !result.topology_candidate_scan.best_passed_hard_thresholds
    @test result.subset_search_topology_scan.status == "all_subset_topologies_fail_forcing_gate"
    @test result.subset_search_topology_scan.candidate_count == 4
    @test result.subset_search_topology_scan.best_forcing_objective_lower_bound > result.objective_target
    @test !result.subset_search_topology_scan.best_passed_hard_thresholds
    @test result.final_objective_target_ratio > 1
    @test result.acceptance_gap_status == "far_above_objective_target"
    @test occursin("coefficient", result.next_required_work)
    @test occursin("flux and heating residuals", result.next_required_work)
    @test result.finite_directional_derivative
    if optional_ad_skipped()
        @test result.reactant_check.status == "skipped_by_RH_SKIP_OPTIONAL_AD_CHECKS"
        @test result.enzyme_check.status == "skipped_by_RH_SKIP_OPTIONAL_AD_CHECKS"
    else
        @test result.reactant_check.status == "passed"
        @test result.enzyme_check.status == "passed"
        @test result.enzyme_check.relative_error <= result.enzyme_check.threshold
    end

    main(; result)
    @test isfile(PREFLIGHT_JSON)
    @test isfile(PREFLIGHT_MD)
    json = read(PREFLIGHT_JSON, String)
    @test occursin("\"status\": \"preflight_ready\"", json)
    @test occursin("\"parameter_count\": 48", json)
    @test occursin("\"trial_line_search\"", json)
    @test occursin("\"optimization_smoke\"", json)
    @test occursin("\"directional_search\"", json)
    @test occursin("\"block_trial_scan\"", json)
    @test occursin("\"coordinate_coefficient_scan\"", json)
    @test occursin("\"candidate_count\": 64", json)
    @test occursin("\"greedy_coordinate_descent\"", json)
    @test occursin("\"greedy_coordinate_descent_checkpoint_used\"", json)
    @test occursin("\"coefficient_joint_direction_scan\"", json)
    @test occursin("\"post_joint_coordinate_refinement\"", json)
    @test occursin("\"finite_difference_coefficient_direction_refinement\"", json)
    @test occursin("\"smooth_objective_coefficient_refinement\"", json)
    @test occursin("\"post_coefficient_weight_refinement\"", json)
    @test occursin("\"final_objective_breakdown\"", json)
    @test occursin("\"targeted_worst_metric_refinement\"", json)
    @test occursin("\"separated_component_refinement\"", json)
    @test occursin("\"pressure_band_table_refinement\"", json)
    @test occursin("\"active_table_entry_refinement\"", json)
    @test occursin("\"targeted_candidate_count\"", json)
    @test occursin("\"constrained_table_optimizer_target\"", json)
    @test occursin("\"warm_started_topology_neighbor_refinement\"", json)
    @test occursin("\"ranked_topology_neighbor_refinement\"", json)
    @test occursin("\"topology_neighbor_scan\"", json)
    @test occursin("\"topology_candidate_scan\"", json)
    @test occursin("\"subset_search_topology_scan\"", json)
    @test occursin("\"all_subset_topologies_fail_forcing_gate\"", json)
    @test occursin("\"all_32x16_topologies_fail_forcing_gate\"", json)
    @test occursin("\"acceptance_gap_status\": \"far_above_objective_target\"", json)
    @test occursin("\"next_required_work\"", json)
    if optional_ad_skipped()
        @test occursin("\"reactant_check\": {\n  \"status\": \"skipped_by_RH_SKIP_OPTIONAL_AD_CHECKS\"", json)
        @test occursin("\"enzyme_check\": {\n  \"status\": \"skipped_by_RH_SKIP_OPTIONAL_AD_CHECKS\"", json)
    else
        @test occursin("\"reactant_check\": {\n  \"status\": \"passed\"", json)
        @test occursin("\"enzyme_check\": {\n  \"status\": \"passed\"", json)
    end
end
