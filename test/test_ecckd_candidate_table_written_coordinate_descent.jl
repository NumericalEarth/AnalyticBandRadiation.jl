if !isdefined(@__MODULE__, :write_sw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateTableWrittenCoordinateDescentValidation
include(joinpath(@__DIR__, "..", "validation",
                 "ecckd_candidate_table_written_coordinate_descent.jl"))
end

@testset "ecCKD candidate table written coordinate descent" begin
    mktempdir() do root
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_sw_flux_fixture(sw_path)
        reference =
            EcckdCandidateTableWrittenCoordinateDescentValidation.official_ecckd_definition_path(
                :shortwave_32,
            )
        candidate = joinpath(root, "written_descent_candidate.nc")

        result =
            EcckdCandidateTableWrittenCoordinateDescentValidation.run_table_written_coordinate_descent(;
                root,
                reference_path = reference,
                candidate_path = candidate,
                scenarios = ("rel-415",),
                initial_parameters = [-0.1, 0.0, -0.05, 0.0],
                deltas = (-0.01, 0.01),
                iterations = 1,
            )
        @test result.status in (
            "written_coordinate_descent_improved",
            "written_coordinate_descent_no_descent",
            "written_coordinate_descent_failed",
        )
        @test result.present_count == 1
        @test result.scenario_count == 1
        @test result.iterations == 1
        @test result.accepted_move_count >= 0
        @test isfinite(result.initial_aggregate_loss)
        @test isfinite(result.final_aggregate_loss)
        @test result.aggregate_loss_reduction_factor >= 0
        @test result.final_improved_count >= 0
        @test length(result.final_rows) == 1
        @test length(result.history) >= 1
        @test isfile(candidate)
        @test occursin("not published-model recovery", result.interpretation)

        markdown =
            EcckdCandidateTableWrittenCoordinateDescentValidation.markdown_table_written_coordinate_descent(
                result,
            )
        @test occursin("ecCKD Candidate Table Written Coordinate Descent", markdown)
        @test occursin("Written Coordinate Descent", markdown)
    end
end
