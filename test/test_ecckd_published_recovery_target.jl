using JSON

module EcckdPublishedRecoveryTargetValidation
include(joinpath(@__DIR__, "..", "validation", "ecckd_published_recovery_target.jl"))
end

@testset "ecCKD published recovery target artifact" begin
    root = normpath(joinpath(@__DIR__, ".."))
    json_path = joinpath(root, "validation", "results", "ecckd_published_recovery_target.json")
    md_path = joinpath(root, "validation", "results", "ecckd_published_recovery_target.md")

    redirect_stdout(devnull) do
        EcckdPublishedRecoveryTargetValidation.ecckd_published_recovery_target_main()
    end

    @test isfile(json_path)
    @test isfile(md_path)

    result = JSON.parsefile(json_path)
    @test result["case"] == "ecckd_published_recovery_target"
    @test result["status"] == "published_recovery_target_ready"
    @test result["model_count"] == 6
    @test length(result["target_models"]) == 6
    @test haskey(result["primary_recovery_targets"], "shortwave_32")
    @test haskey(result["primary_recovery_targets"], "longwave_32")
    @test result["primary_recovery_targets"]["shortwave_32"]["filename"] ==
          "ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc"
    @test result["primary_recovery_targets"]["longwave_32"]["filename"] ==
          "ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc"
    @test result["primary_recovery_targets"]["shortwave_32"]["coefficient_array_count"] == 6
    @test result["primary_recovery_targets"]["longwave_32"]["coefficient_array_count"] == 8
    @test result["primary_recovery_targets"]["shortwave_32"]["coefficient_parameter_count"] > 0
    @test result["primary_recovery_targets"]["longwave_32"]["coefficient_parameter_count"] > 0
    @test result["primary_recovery_targets"]["shortwave_32"]["support_array_count"] >= 4
    @test result["primary_recovery_targets"]["longwave_32"]["support_array_count"] >= 3
    @test result["acceptance_metrics"]["final_objective_target_ratio_max"] <= 1.05
    @test result["acceptance_metrics"]["optical_depth_log_rmse_max"] <= 0.02
    @test occursin("optimizer settings", result["optimizer_only_delta_rule"])

    markdown = read(md_path, String)
    @test occursin("ecCKD Published Recovery Target", markdown)
    @test occursin("Published Targets", markdown)
    @test occursin("ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc", markdown)
    @test occursin("Primary 32-G Targets", markdown)
end
