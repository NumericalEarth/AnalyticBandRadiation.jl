if !isdefined(@__MODULE__, :write_sw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateTableWritebackContinuationValidation
include(joinpath(@__DIR__, "..", "validation",
                 "ecckd_candidate_table_writeback_continuation.jl"))
end

@testset "ecCKD candidate table-writeback continuation" begin
    mktempdir() do root
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_sw_flux_fixture(sw_path)
        reference =
            EcckdCandidateTableWritebackContinuationValidation.official_ecckd_definition_path(
                :shortwave_32,
            )
        candidate = joinpath(root, "scaled_continuation_candidate.nc")

        result =
            EcckdCandidateTableWritebackContinuationValidation.run_table_writeback_continuation(;
                root,
                reference_path = reference,
                candidate_path = candidate,
                iterations = 1,
            )
        @test result.status in
              ("table_writeback_continuation_passed",
               "table_writeback_continuation_failed")
        @test isfile(candidate)
        @test result.iterations == 1
        @test result.accepted_steps >= 0
        @test isfinite(result.initial_in_memory_loss)
        @test isfinite(result.final_in_memory_loss)
        @test result.final_in_memory_loss <= result.initial_in_memory_loss
        @test isfinite(result.baseline_projected_loss)
        @test isfinite(result.written_projected_loss)
        @test result.writeback_loss_reduction_factor >= 0
        @test result.candidate_metrics_status in ("passed", "failed")
        @test occursin("not full published-model recovery", result.interpretation)

        markdown =
            EcckdCandidateTableWritebackContinuationValidation.markdown_table_writeback_continuation(
                result,
            )
        @test occursin("ecCKD Candidate Table-Writeback Continuation", markdown)
        @test occursin("Continuation written SW projected loss", markdown)
    end
end
