using JSON

if !isdefined(@__MODULE__, :write_lw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateTransferOptimizerProbeValidation
include(joinpath(@__DIR__, "..", "validation",
                 "ecckd_candidate_transfer_optimizer_probe.jl"))
end

@testset "ecCKD candidate transfer optimizer probe" begin
    mktempdir() do root
        lw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_lw_fluxes_rel-415.h5"))
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_lw_flux_fixture(lw_path)
        write_sw_flux_fixture(sw_path)
        candidate =
            EcckdCandidateTransferOptimizerProbeValidation.official_ecckd_definition_path(
                :shortwave_32,
            )

        result =
            EcckdCandidateTransferOptimizerProbeValidation.run_transfer_optimizer_probe(;
                root,
                candidate_path = candidate,
            )
        @test result.status in ("optimizer_probe_passed", "optimizer_probe_failed")
        @test result.parameter_count == 4
        @test length(result.parameter_names) == 4
        @test isfinite(result.initial_loss)
        @test isfinite(result.final_loss)
        @test result.final_loss <= result.initial_loss
        @test result.gradient_norm >= 0
        @test result.gradient_method ==
              "central_finite_difference_on_band_projected_transfer_loss"
        @test occursin("not a published CKD-definition recovery", result.interpretation)

        markdown =
            EcckdCandidateTransferOptimizerProbeValidation.markdown_transfer_optimizer_probe(
                result,
            )
        @test occursin("ecCKD Candidate Transfer Optimizer Probe", markdown)
        @test occursin("Optimizer Probe", markdown)
        @test occursin("Initial projected loss", markdown)
    end
end
