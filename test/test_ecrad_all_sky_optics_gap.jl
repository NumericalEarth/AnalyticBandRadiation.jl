@testset "ecRad all-sky optics gap artifact" begin
    script = joinpath(@__DIR__, "..", "validation", "ecrad_all_sky_optics_gap.jl")
    include(script)
    result = markdown_all_sky_optics_gap((
        rows = [(
            candidate_kind = "smoke",
            variable = "od_sw_cloud",
            candidate_variable = "od_sw_cloud",
            units = "1",
            rmse = 1.0,
            max_abs = 1.0,
            mean_bias = 1.0,
            mean_abs = 1.0,
            reference_mean = 1.0,
            candidate_mean = 1.0,
        )],
        candidate_configuration = [(variable = "RH_AEROSOL_OPTICS", value = "false")],
    ))
    @test occursin("ecRad All-Sky Optics Gap", result)
    @test occursin("Candidate Configuration", result)

    json_path = joinpath(@__DIR__, "..", "validation", "results",
                         "ecrad_all_sky_optics_gap.json")
    md_path = joinpath(@__DIR__, "..", "validation", "results",
                       "ecrad_all_sky_optics_gap.md")
    @test isfile(json_path)
    @test isfile(md_path)

    json = read(json_path, String)
    @test occursin("\"case\": \"ecrad_all_sky_optics_gap\"", json)
    @test occursin("\"od_sw_cloud\"", json)
    @test occursin("\"od_lw_cloud\"", json)
end
