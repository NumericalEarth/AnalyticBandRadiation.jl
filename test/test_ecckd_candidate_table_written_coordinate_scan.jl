if !isdefined(@__MODULE__, :write_sw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateTableWrittenCoordinateScanValidation
include(joinpath(@__DIR__, "..", "validation",
                 "ecckd_candidate_table_written_coordinate_scan.jl"))
end

@testset "ecCKD candidate table written coordinate scan" begin
    mktempdir() do root
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_sw_flux_fixture(sw_path)
        reference =
            EcckdCandidateTableWrittenCoordinateScanValidation.official_ecckd_definition_path(
                :shortwave_32,
            )
        candidate = joinpath(root, "written_scan_candidate.nc")

        result =
            EcckdCandidateTableWrittenCoordinateScanValidation.run_table_written_coordinate_scan(;
                root,
                reference_path = reference,
                candidate_path = candidate,
                scenarios = ("rel-415",),
                initial_parameters = [-0.1, 0.0, -0.05, 0.0],
                deltas = (-0.01, 0.01),
            )
        @test result.status in (
            "written_coordinate_scan_improved",
            "written_coordinate_scan_no_descent",
            "written_coordinate_scan_failed",
        )
        @test result.present_count == 1
        @test result.scenario_count == 1
        @test result.tested_move_count >= 1
        @test isfinite(result.seed_aggregate_loss)
        @test isfinite(result.best_aggregate_loss)
        @test result.aggregate_loss_reduction_factor >= 0
        @test result.best_improved_count >= 0
        @test length(result.best_rows) == 1
        @test isfile(candidate)
        @test occursin("not used for acceptance", result.interpretation)

        markdown =
            EcckdCandidateTableWrittenCoordinateScanValidation.markdown_table_written_coordinate_scan(
                result,
            )
        @test occursin("ecCKD Candidate Table Written Coordinate Scan", markdown)
        @test occursin("Written Coordinate Scan", markdown)
    end
end
