if !isdefined(@__MODULE__, :write_sw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateTableWritebackProbeValidation
include(joinpath(@__DIR__, "..", "validation",
                 "ecckd_candidate_table_writeback_probe.jl"))
end

@testset "ecCKD candidate table-writeback probe" begin
    mktempdir() do root
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_sw_flux_fixture(sw_path)
        reference =
            EcckdCandidateTableWritebackProbeValidation.official_ecckd_definition_path(
                :shortwave_32,
            )
        candidate = joinpath(root, "scaled_candidate.nc")
        log_parameters = [-0.1, 0.0, -0.05, 0.0]

        result =
            EcckdCandidateTableWritebackProbeValidation.run_table_writeback_probe(;
                root,
                reference_path = reference,
                candidate_path = candidate,
                log_parameters,
            )
        @test result.status in
              ("table_writeback_probe_passed", "table_writeback_probe_failed")
        @test isfile(candidate)
        @test isfinite(result.baseline_projected_loss)
        @test isfinite(result.written_projected_loss)
        @test result.loss_reduction_factor >= 0
        @test result.baseline_nonfinite_flux_count >= 0
        @test result.written_nonfinite_flux_count >= 0
        @test result.candidate_metrics_status in ("passed", "failed")
        @test occursin("not full published-model recovery", result.interpretation)

        markdown =
            EcckdCandidateTableWritebackProbeValidation.markdown_table_writeback_probe(
                result,
            )
        @test occursin("ecCKD Candidate Table-Writeback Probe", markdown)
        @test occursin("Written SW projected loss", markdown)
    end
end
