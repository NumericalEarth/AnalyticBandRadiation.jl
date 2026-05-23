if !isdefined(@__MODULE__, :write_sw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateTableParameterProbeValidation
include(joinpath(@__DIR__, "..", "validation",
                 "ecckd_candidate_table_parameter_probe.jl"))
end

@testset "ecCKD candidate table-parameter probe" begin
    mktempdir() do root
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_sw_flux_fixture(sw_path)
        candidate =
            EcckdCandidateTableParameterProbeValidation.official_ecckd_definition_path(
                :shortwave_32,
            )

        result =
            EcckdCandidateTableParameterProbeValidation.run_table_parameter_probe(;
                root,
                candidate_path = candidate,
            )
        @test result.status in
              ("table_parameter_probe_passed", "table_parameter_probe_failed")
        @test result.parameter_count == 4
        @test length(result.parameter_names) == 4
        @test isfinite(result.initial_loss)
        @test isfinite(result.final_loss)
        @test result.final_loss <= result.initial_loss
        @test result.gradient_norm >= 0
        @test result.gradient_method ==
              "central_finite_difference_on_in_memory_ckd_table_scales"
        @test occursin("not a full CKD-definition writeback",
                       result.interpretation)

        markdown =
            EcckdCandidateTableParameterProbeValidation.markdown_table_parameter_probe(
                result,
            )
        @test occursin("ecCKD Candidate Table-Parameter Probe", markdown)
        @test occursin("Initial SW projected loss", markdown)
    end
end
