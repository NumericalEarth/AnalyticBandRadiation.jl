using JSON

module EcckdRecoveryMetricsValidation
include(joinpath(@__DIR__, "..", "validation", "ecckd_recovery_metrics.jl"))
end

@testset "official ecCKD recovery metrics artifact" begin
    root = normpath(joinpath(@__DIR__, ".."))
    json_path = joinpath(root, "validation", "results", "ecckd_recovery_metrics.json")
    md_path = joinpath(root, "validation", "results", "ecckd_recovery_metrics.md")
    redirect_stdout(devnull) do
        EcckdRecoveryMetricsValidation.ecckd_recovery_metrics_main()
    end

    @test isfile(json_path)
    @test isfile(md_path)
    output = read(md_path, String)
    @test occursin("ecCKD Recovery Metrics", output)
    @test occursin("Status: **passed**", output)
    @test occursin("published_self_recovery_sanity", output)
    @test occursin("Worst log-coefficient RMSE", output)

    result = JSON.parsefile(json_path)
    @test result["case"] == "ecckd_recovery_metrics"
    @test result["status"] == "passed"
    @test result["recovery_mode"] == "published_self_recovery_sanity"

    metrics = result["metrics"]
    @test metrics["status"] == "passed"
    @test metrics["kind"] == "shortwave"
    @test metrics["coefficient_count"] == 6
    @test sort([coefficient["name"] for coefficient in metrics["coefficients"]]) == [
        "ch4_molar_absorption_coeff",
        "co2_molar_absorption_coeff",
        "composite_molar_absorption_coeff",
        "h2o_molar_absorption_coeff",
        "n2o_molar_absorption_coeff",
        "o3_molar_absorption_coeff",
    ]
    @test metrics["worst_log_coefficient_rmse"] == 0.0
    @test metrics["gpoint_weight_max_abs_error"] == 0.0
    @test metrics["band_weight_max_abs_error"] == 0.0
    @test metrics["thresholds"]["log_coefficient_rmse"] == 1.0e-3
end
