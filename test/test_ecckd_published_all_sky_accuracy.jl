using JSON

include(joinpath(@__DIR__, "..", "validation", "ecckd_published_all_sky_accuracy.jl"))

@testset "published ecCKD all-sky accuracy" begin
    result = main()
    @test isfile(PUBLISHED_ALL_SKY_ACCURACY_JSON)
    @test isfile(PUBLISHED_ALL_SKY_ACCURACY_MD)
    artifact = JSON.parsefile(PUBLISHED_ALL_SKY_ACCURACY_JSON)
    @test artifact["case"] == "ecckd_published_all_sky_accuracy"
    @test artifact["status"] == "passed"
    @test artifact["model_count"] == 6
    @test length(artifact["models"]) == 6
    @test artifact["passed_count"] == artifact["model_count"]
    @test Set(model["ng_lw"] for model in artifact["models"]) == Set([32, 64])
    @test Set(model["ng_sw"] for model in artifact["models"]) == Set([32, 64, 96])
    @test all(model["hard_objective"] > 0 for model in artifact["models"])
    @test all(haskey(model, "toa_forcing_max_abs") for model in artifact["models"])
    @test all(haskey(model, "surface_forcing_max_abs") for model in artifact["models"])
    @test all(haskey(model, "component_boundary_errors") for model in artifact["models"])
    @test all(length(model["component_boundary_errors"]) == 4 for model in artifact["models"])
    @test all(Set(error["component"] for error in model["component_boundary_errors"]) ==
              Set(["lw", "sw"]) for model in artifact["models"])
    @test all(model["passed"] for model in artifact["models"])

    baseline = first(model for model in artifact["models"]
                     if model["case"] == "ecckd_32x32_all_sky_tropical_column")
    @test baseline["passed"]
    @test baseline["toa_forcing_max_abs"] < 0.3
    @test baseline["surface_forcing_max_abs"] < 0.3

    wide_rows = filter(model -> model["ng_lw"] == 64 || model["ng_sw"] > 32,
                       artifact["models"])
    @test !isempty(wide_rows)
    @test any(model -> model["limiting_metric"] in
              ("toa_forcing_abs_error", "surface_forcing_abs_error"),
              wide_rows)

    md = read(PUBLISHED_ALL_SKY_ACCURACY_MD, String)
    @test occursin("Published ecCKD All-Sky Accuracy", md)
    @test occursin("Tripleclouds/aerosol", md)
    @test occursin("Hard objective", md)
    @test occursin("LW TOA", md)
    @test occursin("SW surface", md)
end
