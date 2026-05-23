if !isdefined(@__MODULE__, :write_sw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateTableHumiditySplitDescentValidation
include(joinpath(@__DIR__, "..", "validation",
                 "ecckd_candidate_table_humidity_split_descent.jl"))
end

@testset "ecCKD candidate table humidity-split descent" begin
    mktempdir() do root
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_sw_flux_fixture(sw_path)
        reference =
            EcckdCandidateTableHumiditySplitDescentValidation.official_ecckd_definition_path(
                :shortwave_32,
            )
        candidate = joinpath(root, "humidity_split_descent_candidate.nc")

        result =
            EcckdCandidateTableHumiditySplitDescentValidation.run_table_humidity_split_descent(;
                root,
                reference_path = reference,
                candidate_path = candidate,
                scenarios = ("rel-415",),
                initial_parameters = [-0.1, -0.1, 0.0, -0.05, 0.0],
                deltas = (-0.01, 0.01),
                iterations = 1,
            )
        @test result.status in (
            "humidity_split_descent_improved",
            "humidity_split_descent_no_descent",
            "humidity_split_descent_failed",
        )
        @test result.present_count == 1
        @test result.scenario_count == 1
        @test result.tested_move_count >= 1
        @test result.accepted_move_count >= 0
        @test isfinite(result.initial_aggregate_loss)
        @test isfinite(result.final_aggregate_loss)
        @test result.aggregate_loss_reduction_factor >= 0
        @test result.worst_loss_ratio_reduction_factor >= 0
        @test result.final_improved_count >= 0
        @test length(result.final_rows) == 1
        @test isfile(candidate)
        @test occursin("not published-model recovery", result.interpretation)

        markdown =
            EcckdCandidateTableHumiditySplitDescentValidation.markdown_table_humidity_split_descent(
                result,
            )
        @test occursin("ecCKD Candidate Table Humidity-Split Descent", markdown)
        @test occursin("Descent Summary", markdown)
    end
end
