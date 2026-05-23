using JSON

module EcckdPublishedRecoveryVectorTrainingValidation
include(joinpath(@__DIR__, "..", "validation", "ecckd_published_recovery_vector_training.jl"))
end

@testset "ecCKD published recovery vector training artifact" begin
    root = normpath(joinpath(@__DIR__, ".."))
    json_path = joinpath(root, "validation", "results",
                         "ecckd_published_recovery_vector_training.json")
    md_path = joinpath(root, "validation", "results",
                       "ecckd_published_recovery_vector_training.md")

    redirect_stdout(devnull) do
        EcckdPublishedRecoveryVectorTrainingValidation.ecckd_published_recovery_vector_training_main()
    end

    @test isfile(json_path)
    @test isfile(md_path)

    result = JSON.parsefile(json_path)
    @test result["case"] == "ecckd_published_recovery_vector_training"
    @test result["status"] == "passed"
    @test result["array_count"] == 9
    @test result["parameter_count"] == 204896
    @test result["trained_parameter_count"] == 64
    @test result["iterations"] == 36
    @test result["final_loss"] < result["initial_loss"]
    @test result["loss_reduction_factor"] > 1.0e10
    @test !result["enzyme_requested"]
    @test !result["enzyme_used_for_training"]
    @test !result["reactant_check_requested"]
    @test result["reactant_compile_check"]["status"] == "skipped"
    @test result["recovery_metrics"]["status"] == "passed"
    @test result["recovery_metrics"]["worst_log_coefficient_rmse"] < 1.0e-3
    @test result["recovery_metrics"]["worst_p99_relative_coefficient_error"] < 1.0e-2
    @test result["recovery_metrics"]["gpoint_weight_max_abs_error"] < 1.0e-6
    @test result["roundtrip_error"]["l1_relative_error"] < 1.0e-6
    @test occursin("deterministic slice", result["interpretation"])
    @test occursin("not the full original CKDMIP flux objective", result["interpretation"])

    markdown = read(md_path, String)
    @test occursin("ecCKD Published Recovery Vector Training", markdown)
    @test occursin("analytic quadratic gradient", markdown)
    @test occursin("Reactant compile check", markdown)
end
