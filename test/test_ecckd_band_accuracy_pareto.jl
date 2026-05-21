using JSON

@testset "ecCKD band-count accuracy Pareto artifact" begin
    root = normpath(joinpath(@__DIR__, ".."))
    script = joinpath(root, "validation", "ecckd_band_accuracy_pareto.jl")
    test_project = joinpath(@__DIR__, "Project.toml")
    output = read(`$(Base.julia_cmd()) --project=$test_project $script`, String)

    @test occursin("ecCKD Band-Count Accuracy Pareto", output)
    @test occursin("Status: **passed**", output)
    @test occursin("Published ecCKD inventory entries: 6", output)

    md_path = joinpath(root, "validation", "results", "ecckd_band_accuracy_pareto.md")
    json_path = joinpath(root, "validation", "results", "ecckd_band_accuracy_pareto.json")
    csv_path = joinpath(root, "validation", "results", "ecckd_band_accuracy_pareto.csv")
    svg_path = joinpath(root, "validation", "results", "ecckd_band_accuracy_pareto.svg")
    @test isfile(md_path)
    @test isfile(json_path)
    @test isfile(csv_path)
    @test isfile(svg_path)

    result = JSON.parsefile(json_path)
    @test result["status"] == "passed"
    @test result["point_count"] >= 22
    @test result["published_inventory_count"] == 6
    @test result["passed_point_count"] >= 1
    @test any(row -> row["ng_lw"] == 32 && row["ng_sw"] == 32 && row["passed"],
              result["accuracy_points"])
    @test any(row -> row["ng_lw"] == 32 && row["ng_sw"] == 16,
              result["accuracy_points"])

    csv = read(csv_path, String)
    @test occursin("source,label,ng_lw,ng_sw,total_gpoints,passed", csv)
    @test length(split(chomp(csv), "\n")) == result["point_count"] + 1

    svg = read(svg_path, String)
    @test occursin("<svg", svg)
    @test occursin("ecCKD Accuracy vs Total G-points", svg)
    @test occursin("0.3 W m^-2 hard threshold", svg)
end
