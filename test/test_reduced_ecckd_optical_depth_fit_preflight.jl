using Test

include(joinpath(@__DIR__, "..", "validation", "reduced_ecckd_optical_depth_fit_preflight.jl"))

@testset "reduced ecCKD optical-depth fit preflight artifact" begin
    result = optical_depth_fit_preflight()
    @test result.case == "reduced_ecckd_optical_depth_fit_preflight"
    @test result.status == "optical_depth_refit_target_ready"
    @test result.ng_sw == 16
    @test result.sample_count > 0
    @test result.baseline.rmse > 0
    @test result.fitted.rmse < result.baseline.rmse
    @test result.component_fitted.rmse < result.fitted.rmse
    @test result.relative_rmse_reduction > 0
    @test result.component_relative_rmse_reduction > result.relative_rmse_reduction
    @test length(result.fitted_scales) == 16
    @test size(result.component_scales) == (16, 3)
    @test result.scale_min > 0
    @test result.scale_max >= result.scale_min
    @test result.component_scale_max >= result.component_scale_min
    @test result.flux_baseline_objective > 1
    @test result.flux_scaled_objective > result.flux_baseline_objective
    @test !result.flux_scaled_improved
    @test !result.flux_scaled_passed_hard_thresholds
    @test result.flux_component_scaled_objective > result.flux_baseline_objective
    @test !result.flux_component_scaled_improved
    @test !result.flux_component_scaled_passed_hard_thresholds
    @test result.coefficient_table_fit.parameter_count_per_g > 0
    @test result.coefficient_table_fit.sample_count == result.sample_count
    @test result.coefficient_table_fit.physical_target_optical_depth.rmse < 1.0e-10
    @test result.coefficient_table_fit.physical_target_flux_objective > 1
    @test !result.coefficient_table_fit.physical_target_flux_passed_hard_thresholds
    @test result.coefficient_table_fit.raw_least_squares_optical_depth.rmse <
          result.component_fitted.rmse
    @test result.coefficient_table_fit.clipped_model_optical_depth.rmse >
          result.baseline.rmse
    @test result.coefficient_table_fit.clipped_parameter_count > 0
    @test result.coefficient_table_fit.flux_objective > result.flux_baseline_objective
    @test !result.coefficient_table_fit.flux_improved
    @test !result.coefficient_table_fit.flux_passed_hard_thresholds
    @test length(result.case_rows) == length(REDUCED_CASES)
    @test occursin("optical-depth targets", result.next_required_work)

    main()
    @test isfile(OPTICAL_DEPTH_FIT_JSON)
    @test isfile(OPTICAL_DEPTH_FIT_MD)
    json = read(OPTICAL_DEPTH_FIT_JSON, String)
    @test occursin("\"case\": \"reduced_ecckd_optical_depth_fit_preflight\"", json)
    @test occursin("\"status\": \"optical_depth_refit_target_ready\"", json)
    @test occursin("\"fitted_scales\"", json)
    @test occursin("\"component_scales\"", json)
    @test occursin("\"relative_rmse_reduction\"", json)
    @test occursin("\"flux_scaled_improved\": false", json)
    @test occursin("\"flux_component_scaled_improved\": false", json)
    @test occursin("\"coefficient_table_fit\"", json)
    markdown = read(OPTICAL_DEPTH_FIT_MD, String)
    @test occursin("Reduced ecCKD Optical-Depth Fit Preflight", markdown)
    @test occursin("Baseline optical-depth RMSE", markdown)
    @test occursin("Component-fitted optical-depth RMSE", markdown)
    @test occursin("Scaled flux objective", markdown)
    @test occursin("Physical projected table optical-depth RMSE", markdown)
    @test occursin("Coefficient-table raw LS optical-depth RMSE", markdown)
    @test occursin("Coefficient-table clipped-model optical-depth RMSE", markdown)
end
