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
    @test occursin("Reference case", result)
    @test occursin("ecRad properties", result)
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

    json_32x64_path = joinpath(@__DIR__, "..", "validation", "results",
                               "ecrad_all_sky_optics_gap_32x64.json")
    md_32x64_path = joinpath(@__DIR__, "..", "validation", "results",
                             "ecrad_all_sky_optics_gap_32x64.md")
    @test isfile(json_32x64_path)
    @test isfile(md_32x64_path)

    json_32x64 = read(json_32x64_path, String)
    @test occursin("ecckd_32x64_all_sky_tropical_column.nc", json_32x64)
    @test occursin("ecckd-1.2_sw_climate_window-64b_ckd-definition.nc", json_32x64)

    gate_json_path = joinpath(@__DIR__, "..", "validation", "results",
                              "ecrad_all_sky_optics_gap_32x64_gate.json")
    gate_md_path = joinpath(@__DIR__, "..", "validation", "results",
                            "ecrad_all_sky_optics_gap_32x64_gate.md")
    @test isfile(gate_json_path)
    @test isfile(gate_md_path)

    gate_json = read(gate_json_path, String)
    @test occursin("\"RH_AEROSOL_OPTICS\"", gate_json)
    @test occursin("\"RH_IFS_AEROSOL_TABLE_OPTICS\"", gate_json)
    @test occursin("\"value\": \"true\"", gate_json)
    @test occursin("clear_region_current", gate_json)
    @test occursin("cloudy_region_delta_scaled", gate_json)
    @test occursin("od_sw+od_sw_cloud", gate_json)
    @test occursin("combined_ssa_sw", gate_json)
    @test occursin("boundary_materialized", gate_json)
    @test occursin("incoming_sw", gate_json)
    @test occursin("sw_albedo_direct", gate_json)
end
