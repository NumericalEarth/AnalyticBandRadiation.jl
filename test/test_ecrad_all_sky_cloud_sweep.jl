@testset "ecRad all-sky cloud sweep artifact" begin
    script = joinpath(@__DIR__, "..", "validation", "ecrad_all_sky_cloud_sweep.jl")
    include(script)
    result = markdown_cloud_sweep_report((
        best_trial = "smoke",
        best_configuration = [(variable = "RH_AEROSOL_OPTICS", value = "false")],
        trials = [(
            name = "smoke",
            all_sky_worst_threshold_ratio = 1.0,
            toa_forcing_abs_error_w_m2 = 1.0,
            surface_forcing_abs_error_w_m2 = 1.0,
            toa_sw_cloud_effect_max_abs_w_m2 = 1.0,
            surface_sw_cloud_effect_max_abs_w_m2 = 1.0,
            profile_sw_cloud_effect_max_abs_w_m2 = 1.0,
            heating_cloud_effect_max_abs_k_day = 1.0,
        )],
    ))
    @test occursin("ecRad All-Sky Cloud Sweep", result)
    @test occursin("Best Configuration", result)

    json_path = joinpath(@__DIR__, "..", "validation", "results",
                         "ecrad_all_sky_cloud_sweep.json")
    md_path = joinpath(@__DIR__, "..", "validation", "results",
                       "ecrad_all_sky_cloud_sweep.md")
    @test isfile(json_path)
    @test isfile(md_path)

    json = read(json_path, String)
    @test occursin("\"case\": \"ecrad_all_sky_cloud_sweep\"", json)
    @test occursin("\"best_trial\":", json)
    @test occursin("\"trials\":", json)
end
