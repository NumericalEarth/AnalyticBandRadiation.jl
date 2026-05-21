using Test

include(joinpath(@__DIR__, "..", "validation", "official_ecckd_training.jl"))

@testset "official reduced ecCKD gas-optics training artifact" begin
    preflight = reduced_optimization_preflight()
    result = official_training_report(preflight)
    @test result.case == "official_reduced_ecckd_gas_optics_training"
    @test result.status == "passed"
    @test result.training_source == "official ecCKD tabulated gas-optics reduced hard-gate objective"
    @test result.parameter_count == 48
    @test result.ng_sw == 16
    @test result.final_objective < result.initial_objective
    @test result.objective_reduction > 0
    @test result.objective_ratio < 1
    @test !result.hard_accuracy_target_met
    @test result.final_objective_target_ratio > 1
    @test result.reactant_check.status == "passed"
    @test result.enzyme_check.status == "passed"
    @test result.topology_candidate_scan.status == "all_32x16_topologies_fail_forcing_gate"
    @test occursin("does not close the reduced-accuracy gate", result.notes)

    main()
    @test isfile(OFFICIAL_TRAINING_JSON)
    @test isfile(OFFICIAL_TRAINING_MD)
    json = read(OFFICIAL_TRAINING_JSON, String)
    @test occursin("\"case\": \"official_reduced_ecckd_gas_optics_training\"", json)
    @test occursin("\"status\": \"passed\"", json)
    @test occursin("\"hard_accuracy_target_met\": false", json)
    @test occursin("\"training_source\": \"official ecCKD tabulated gas-optics reduced hard-gate objective\"", json)
    @test occursin("\"topology_candidate_scan\"", json)
end
