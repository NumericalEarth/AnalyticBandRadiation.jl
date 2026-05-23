if !isdefined(@__MODULE__, :write_sw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateTableMultisampleOptimizerValidation
include(joinpath(@__DIR__, "..", "validation",
                 "ecckd_candidate_table_multisample_optimizer.jl"))
end

@testset "ecCKD candidate table multisample optimizer" begin
    mktempdir() do root
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_sw_flux_fixture(sw_path)
        reference =
            EcckdCandidateTableMultisampleOptimizerValidation.official_ecckd_definition_path(
                :shortwave_32,
            )
        candidate = joinpath(root, "multisample_candidate.nc")

        result =
            EcckdCandidateTableMultisampleOptimizerValidation.run_table_multisample_optimizer(;
                root,
                reference_path = reference,
                candidate_path = candidate,
                scenarios = ("rel-415",),
                iterations = 1,
            )
        @test result.status in (
            "table_multisample_optimizer_passed",
            "table_multisample_optimizer_seed_passed",
            "table_multisample_optimizer_partial",
            "table_multisample_optimizer_failed",
        )
        @test result.present_count == 1
        @test result.scenario_count == 1
        @test result.iterations == 1
        @test result.accepted_steps >= 0
        @test isfinite(result.initial_aggregate_loss)
        @test isfinite(result.final_aggregate_loss)
        @test result.aggregate_loss_reduction_factor >= 0
        @test result.best_written_checkpoint_index >= 0
        @test result.selected_written_aggregate_loss >= 0
        @test isfile(candidate)
        @test result.writeback_status in (
            "table_writeback_multisample_passed",
            "table_writeback_multisample_partial",
        )
        @test length(result.writeback_rows) == 1
        @test occursin("not published-model recovery", result.interpretation)

        markdown =
            EcckdCandidateTableMultisampleOptimizerValidation.markdown_table_multisample_optimizer(
                result,
            )
        @test occursin("ecCKD Candidate Table Multisample Optimizer", markdown)
        @test occursin("Aggregate Optimizer", markdown)
    end
end
