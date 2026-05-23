@testset "ecRad reference-optics solver gap artifact" begin
    script = joinpath(@__DIR__, "..", "validation",
                      "ecrad_reference_optics_solver_gap.jl")
    include(script)

    default_json = joinpath(@__DIR__, "..", "validation", "results",
                            "ecrad_reference_optics_solver_gap.json")
    default_md = joinpath(@__DIR__, "..", "validation", "results",
                          "ecrad_reference_optics_solver_gap.md")
    diagnostic_json = joinpath(@__DIR__, "..", "validation", "results",
                               "ecrad_reference_optics_solver_gap_32x64.json")
    diagnostic_md = joinpath(@__DIR__, "..", "validation", "results",
                             "ecrad_reference_optics_solver_gap_32x64.md")

    @test isfile(default_json)
    @test isfile(default_md)
    @test isfile(diagnostic_json)
    @test isfile(diagnostic_md)

    default = read(default_json, String)
    diagnostic = read(diagnostic_json, String)
    @test occursin("\"case\": \"ecrad_reference_optics_solver_gap\"", default)
    @test occursin("ecrad_meridian_ecckd_tc_out_REFERENCE.nc", default)
    @test occursin("ecckd_32x64_all_sky_tropical_column.nc", diagnostic)
    @test occursin("ecrad_meridian_ecckd_32x64_all_sky_props_out.nc",
                   diagnostic)
    @test occursin("tripleclouds_alpha_p2", diagnostic)

    report = markdown_report((
        reference = "reference.nc",
        properties = "properties.nc",
        ecrad_output = "output.nc",
        modes = [(
            mode = "smoke",
            sw_up_rmse = 1.0,
            sw_down_rmse = 1.0,
            toa_net_abs_error = 1.0,
            surface_net_abs_error = 1.0,
            toa_net_mean_bias = 1.0,
            surface_net_mean_bias = 1.0,
            output_toa_net_abs_error = nothing,
            output_surface_net_abs_error = nothing,
            reference_output_toa_net_abs_error = nothing,
            reference_output_surface_net_abs_error = nothing,
            clear_direct_max_abs = nothing,
        )],
    ))
    @test occursin("Reference case", report)
    @test occursin("ecRad output", report)
    @test occursin("Output TOA net max abs", report)
    @test occursin("reference_output_toa_net_abs_error", diagnostic)
    @test occursin("5.097602388559608e-5", diagnostic)
    @test occursin("4.57763671875e-5", diagnostic)
end
