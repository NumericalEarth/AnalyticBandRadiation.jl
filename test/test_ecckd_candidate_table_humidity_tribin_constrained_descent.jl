if !isdefined(@__MODULE__, :write_sw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateTableHumidityTribinConstrainedDescentValidation
include(joinpath(@__DIR__, "..", "validation",
                 "ecckd_candidate_table_humidity_tribin_constrained_descent.jl"))
end

@testset "ecCKD candidate table humidity-tribin constrained descent" begin
    mktempdir() do root
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_sw_flux_fixture(sw_path)
        reference =
            EcckdCandidateTableHumidityTribinConstrainedDescentValidation.official_ecckd_definition_path(
                :shortwave_32,
            )
        candidate = joinpath(root, "humidity_tribin_constrained_candidate.nc")

        result =
            EcckdCandidateTableHumidityTribinConstrainedDescentValidation.run_table_humidity_tribin_constrained_descent(;
                root,
                reference_path = reference,
                candidate_path = candidate,
                scenarios = ("rel-415",),
                initial_parameters = [-0.1, -0.1, -0.1, 0.0, -0.05, 0.0],
                worst_ratio_tolerances = (0.0, 1.0e-4),
                deltas = (-0.01, 0.01),
                iterations = 2,
            )
        @test result.status in (
            "humidity_tribin_constrained_descent_improved",
            "humidity_tribin_constrained_descent_no_descent",
        )
        @test result.present_count == 1
        @test result.scenario_count == 1
        @test result.tested_tolerance_count == 2
        @test result.best_tested_move_count >= 6
        @test result.best_accepted_move_count >= 0
        @test result.best_aggregate_loss >= 0
        @test result.best_worst_loss_ratio >= 0
        @test result.best_final_improved_count >= 0
        @test length(result.best_final_rows) == 1
        @test isfile(candidate)
        @test occursin("not published-model recovery", result.interpretation)

        markdown =
            EcckdCandidateTableHumidityTribinConstrainedDescentValidation.markdown_table_humidity_tribin_constrained_descent(
                result,
            )
        @test occursin("ecCKD Candidate Table Humidity-Tribin Constrained Descent",
                       markdown)
        @test occursin("Tolerance Sweep", markdown)
    end
end
