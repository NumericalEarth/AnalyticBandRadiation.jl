using Test
using JSON

include(joinpath(@__DIR__, "..", "validation", "ecckd_published_model_accuracy.jl"))

@testset "published ecCKD model accuracy" begin
    result = main()
    @test isfile(PUBLISHED_MODEL_ACCURACY_JSON)
    @test isfile(PUBLISHED_MODEL_ACCURACY_MD)
    artifact = JSON.parsefile(PUBLISHED_MODEL_ACCURACY_JSON)
    @test artifact["case"] == "ecckd_published_model_accuracy"
    @test artifact["status"] == "failed_threshold"
    @test length(artifact["models"]) == 3
    @test length(artifact["isolation_diagnostics"]) == 3
    models = Dict(model["label"] => model for model in artifact["models"])
    model_32 = models["official ecCKD 1.0 32-LW x 32-SW climate model"]
    model_64 = models["official ecCKD 1.2 64-LW x 64-SW climate model"]
    model_96 = models["official ecCKD 1.2/1.4 64-LW x 96-SW climate/vfine model"]
    @test model_32["ng_lw"] == 32
    @test model_32["ng_sw"] == 32
    @test model_64["ng_lw"] == 64
    @test model_64["ng_sw"] == 64
    @test model_96["ng_lw"] == 64
    @test model_96["ng_sw"] == 96
    @test model_32["passed_hard_thresholds"]
    @test model_32["hard_objective"]["value"] <= 1.0
    @test !model_64["passed_hard_thresholds"]
    @test !model_96["passed_hard_thresholds"]
    @test model_64["hard_objective"]["value"] > 1.0
    @test model_96["hard_objective"]["value"] > 1.0
    diagnostics = Dict(model["label"] => model for model in artifact["isolation_diagnostics"])
    sw64 = diagnostics["component isolation: published 32-LW x 64-SW"]
    sw96 = diagnostics["component isolation: published 32-LW x 96-SW"]
    lw64 = diagnostics["component isolation: published 64-LW x 32-SW"]
    @test sw64["ng_lw"] == 32
    @test sw64["ng_sw"] == 64
    @test sw96["ng_lw"] == 32
    @test sw96["ng_sw"] == 96
    @test lw64["ng_lw"] == 64
    @test lw64["ng_sw"] == 32
    @test !sw64["passed_hard_thresholds"]
    @test !sw96["passed_hard_thresholds"]
    @test !lw64["passed_hard_thresholds"]
    @test sw64["hard_objective"]["metric"] == "surface_forcing"
    @test sw96["hard_objective"]["metric"] == "surface_forcing"
    @test lw64["hard_objective"]["metric"] == "heating_rate_max_abs"
    md = read(PUBLISHED_MODEL_ACCURACY_MD, String)
    @test occursin("64-LW x 96-SW", md)
    @test occursin("Mixed-component isolation diagnostics", md)
    @test occursin("Status: **failed_threshold**", md)
end
