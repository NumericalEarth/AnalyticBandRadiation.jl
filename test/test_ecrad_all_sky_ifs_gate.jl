using Test

include(joinpath(@__DIR__, "..", "validation", "ecrad_all_sky_ifs_gate.jl"))

@testset "ecRad all-sky IFS gate artifact" begin
    result = all_sky_ifs_gate()
    @test result.case == "ecrad_all_sky_ifs_gate"
    @test result.status == "passed"
    @test result.accuracy_gate_ecckd_all_sky_passed
    @test result.cloud_scattering_tables_passed
    @test result.cloud_sweep.passed
    @test result.cloud_sweep.best_worst_threshold_ratio <= 1
    @test occursin("tripleclouds", result.cloud_sweep.best_trial)
    @test result.reference_optics_solver.passed
    @test result.reference_optics_solver.toa_net_abs_error <=
          result.reference_optics_solver.threshold
    @test result.reference_optics_solver.surface_net_abs_error <=
          result.reference_optics_solver.threshold

    main()
    @test isfile(ALL_SKY_IFS_GATE_JSON)
    @test isfile(ALL_SKY_IFS_GATE_MD)
    json = read(ALL_SKY_IFS_GATE_JSON, String)
    @test occursin("\"case\": \"ecrad_all_sky_ifs_gate\"", json)
    @test occursin("\"status\": \"passed\"", json)
    @test occursin("\"accuracy_gate_ecckd_all_sky_passed\": true", json)
    @test occursin("\"cloud_scattering_tables_passed\": true", json)
end
