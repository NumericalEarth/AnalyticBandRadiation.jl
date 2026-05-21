using Test

include(joinpath(@__DIR__, "..", "validation", "ecckd_32b_baseline_check.jl"))

@testset "ecCKD 32b baseline check artifact" begin
    result = ecckd_32b_baseline_check()
    @test result.case == "ecckd_32b_baseline_check"
    @test result.status == "passed"
    @test occursin("32b", result.longwave_definition)
    @test occursin("32b", result.shortwave_definition)
    @test result.lw_gpoints == 32
    @test result.sw_gpoints == 32
    @test result.weights_normalized
    @test result.full_32x32_reduced_accuracy_anchor_passed

    main()
    @test isfile(ECCKD_32B_JSON)
    @test isfile(ECCKD_32B_MD)
    json = read(ECCKD_32B_JSON, String)
    @test occursin("\"case\": \"ecckd_32b_baseline_check\"", json)
    @test occursin("\"status\": \"passed\"", json)
    @test occursin("\"lw_gpoints\": 32", json)
    @test occursin("\"sw_gpoints\": 32", json)
    @test occursin("\"full_32x32_reduced_accuracy_anchor_passed\": true", json)
end
