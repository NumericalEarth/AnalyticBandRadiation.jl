using JSON

module EcckdPublishedTrainingManifestValidation
include(joinpath(@__DIR__, "..", "validation", "ecckd_published_training_manifest.jl"))
end

@testset "ecCKD published training manifest" begin
    root = normpath(joinpath(@__DIR__, ".."))
    json_path = joinpath(root, "validation", "results", "ecckd_published_training_manifest.json")
    md_path = joinpath(root, "validation", "results", "ecckd_published_training_manifest.md")
    redirect_stdout(devnull) do
        EcckdPublishedTrainingManifestValidation.ecckd_published_training_manifest_main()
    end

    @test isfile(json_path)
    @test isfile(md_path)
    @test occursin("ecCKD Published Training Manifest", read(md_path, String))

    result = JSON.parsefile(json_path)
    @test result["case"] == "ecckd_published_training_manifest"
    @test result["status"] == "passed"
    @test result["config"]["training_code"] == "evaluation1"
    @test result["config"]["evaluation_code"] == "evaluation2"
    @test result["config"]["training_both"] == "no"
    @test all(file -> file["exists"], result["source_files"])

    scripts = Dict(script["path"] => script for script in result["optimization_scripts"])
    @test haskey(scripts, "test/optimize_lut_lw.sh")
    @test haskey(scripts, "test/optimize_lut_sw.sh")
    @test occursin("prior_error=8.0", scripts["test/optimize_lut_lw.sh"]["selected_common_options"])
    @test occursin("bounded_optimization=0", scripts["test/optimize_lut_sw.sh"]["selected_common_options"])

    lw_modes = Set(vcat((mode["names"] for mode in scripts["test/optimize_lut_lw.sh"]["modes"])...))
    sw_modes = Set(vcat((mode["names"] for mode in scripts["test/optimize_lut_sw.sh"]["modes"])...))
    @test "climate" in lw_modes
    @test "all-in-one" in lw_modes
    @test "relative-base" in sw_modes
end
