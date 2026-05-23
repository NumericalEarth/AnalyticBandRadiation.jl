if !isdefined(@__MODULE__, :write_sw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateTableHumidityTribinWeightedDescentValidation
include(joinpath(@__DIR__, "..", "validation",
                 "ecckd_candidate_table_humidity_tribin_weighted_descent.jl"))
end

@testset "ecCKD candidate table humidity-tribin weighted descent" begin
    mktempdir() do root
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_sw_flux_fixture(sw_path)
        reference =
            EcckdCandidateTableHumidityTribinWeightedDescentValidation.official_ecckd_definition_path(
                :shortwave_32,
            )
        candidate = joinpath(root, "humidity_tribin_weighted_candidate.nc")

        result =
            EcckdCandidateTableHumidityTribinWeightedDescentValidation.run_table_humidity_tribin_weighted_descent(;
                root,
                reference_path = reference,
                candidate_path = candidate,
                scenarios = ("rel-415",),
                initial_parameters = [-0.1, -0.1, -0.1, 0.0, -0.05, 0.0],
                scenario_weight_sets = (Dict("rel-415" => 2.0),),
                deltas = (-0.01, 0.01),
                iterations = 2,
            )
        @test result.status in (
            "humidity_tribin_weighted_descent_improved",
            "humidity_tribin_weighted_descent_no_descent",
        )
        @test result.present_count == 1
        @test result.scenario_count == 1
        @test result.tested_weight_set_count == 1
        @test result.best_tested_move_count >= 6
        @test result.best_accepted_move_count >= 0
        @test result.best_weighted_loss_ratio >= 0
        @test result.best_worst_loss_ratio >= 0
        @test result.best_final_improved_count >= 0
        @test length(result.best_final_rows) == 1
        @test isfile(candidate)
        @test occursin("not published-model recovery", result.interpretation)

        markdown =
            EcckdCandidateTableHumidityTribinWeightedDescentValidation.markdown_table_humidity_tribin_weighted_descent(
                result,
            )
        @test occursin("ecCKD Candidate Table Humidity-Tribin Weighted Descent",
                       markdown)
        @test occursin("Weight-Set Sweep", markdown)
    end
end
