using JSON

module EcckdPublishedRecoveryVectorValidation
include(joinpath(@__DIR__, "..", "validation", "ecckd_published_recovery_vector.jl"))
end

@testset "ecCKD published recovery vector artifact" begin
    root = normpath(joinpath(@__DIR__, ".."))
    json_path = joinpath(root, "validation", "results", "ecckd_published_recovery_vector.json")
    md_path = joinpath(root, "validation", "results", "ecckd_published_recovery_vector.md")

    redirect_stdout(devnull) do
        EcckdPublishedRecoveryVectorValidation.ecckd_published_recovery_vector_main()
    end

    @test isfile(json_path)
    @test isfile(md_path)

    result = JSON.parsefile(json_path)
    @test result["case"] == "ecckd_published_recovery_vector"
    @test result["status"] == "passed"
    @test result["array_count"] == 9
    @test result["parameter_count"] == 204896
    @test result["roundtrip_error"]["max_abs_error"] == 0.0
    @test result["roundtrip_error"]["l1_relative_error"] == 0.0
    @test result["recovery_metrics"]["status"] == "passed"
    @test result["recovery_metrics"]["worst_log_coefficient_rmse"] == 0.0
    @test result["recovery_metrics"]["gpoint_weight_max_abs_error"] == 0.0
    @test any(row -> row["name"] == "gpoint_fraction", result["arrays"])
    @test any(row -> row["name"] == "solar_irradiance", result["arrays"])
    @test any(row -> row["name"] == "rayleigh_molar_scattering_coeff", result["arrays"])

    markdown = read(md_path, String)
    @test occursin("ecCKD Published Recovery Vector", markdown)
    @test occursin("Round-trip max abs error", markdown)
    @test occursin("optimizer handoff vector", markdown)
end
