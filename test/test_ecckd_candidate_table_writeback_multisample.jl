if !isdefined(@__MODULE__, :write_sw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateTableWritebackMultisampleValidation
include(joinpath(@__DIR__, "..", "validation",
                 "ecckd_candidate_table_writeback_multisample.jl"))
end

@testset "ecCKD candidate table-writeback multisample" begin
    mktempdir() do root
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_sw_flux_fixture(sw_path)
        reference =
            EcckdCandidateTableWritebackMultisampleValidation.official_ecckd_definition_path(
                :shortwave_32,
            )
        candidate = joinpath(root, "scaled_candidate.nc")
        EcckdCandidateTableWritebackMultisampleValidation.write_scaled_shortwave_candidate(
            reference,
            candidate,
            [-0.1, 0.0, -0.05, 0.0],
        )

        result =
            EcckdCandidateTableWritebackMultisampleValidation.run_table_writeback_multisample(;
                root,
                reference_path = reference,
                candidate_path = candidate,
                scenarios = ("rel-415",),
            )
        @test result.status in
              ("table_writeback_multisample_passed",
               "table_writeback_multisample_partial")
        @test result.scenario_count == 1
        @test result.present_count == 1
        @test result.improved_count >= 0
        @test length(result.rows) == 1
        @test result.rows[1].present
        @test isfinite(result.rows[1].baseline_loss)
        @test isfinite(result.rows[1].candidate_loss)
        @test result.rows[1].loss_ratio >= 0
        @test occursin("not full published-model recovery", result.interpretation)

        markdown =
            EcckdCandidateTableWritebackMultisampleValidation.markdown_table_writeback_multisample(
                result,
            )
        @test occursin("ecCKD Candidate Table-Writeback Multisample", markdown)
        @test occursin("Scenario Scores", markdown)
    end
end
