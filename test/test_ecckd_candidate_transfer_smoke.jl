using JSON

if !isdefined(@__MODULE__, :write_lw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateTransferSmokeValidation
include(joinpath(@__DIR__, "..", "validation", "ecckd_candidate_transfer_smoke.jl"))
end

@testset "ecCKD candidate transfer smoke" begin
    mktempdir() do root
        lw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_lw_fluxes_rel-415.h5"))
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_lw_flux_fixture(lw_path)
        write_sw_flux_fixture(sw_path)
        candidate =
            EcckdCandidateTransferSmokeValidation.official_ecckd_definition_path(
                :shortwave_32,
            )

        result =
            EcckdCandidateTransferSmokeValidation.run_ecckd_candidate_transfer_smoke(;
                root,
                candidate_path = candidate,
            )
        @test result.status == "candidate_transfer_smoke_ready"
        @test result.layer_count > 0
        @test result.longwave_gpoints == 32
        @test result.shortwave_gpoints == 32
        @test isfinite(result.longwave_optical_depth_max)
        @test isfinite(result.shortwave_optical_depth_max)
        @test isfinite(result.objective_scores.longwave.loss)
        @test isfinite(result.objective_scores.shortwave.loss)
        @test result.objective_scores.longwave.loss >= 0
        @test result.objective_scores.shortwave.loss >= 0
        @test result.longwave_projection_band_count > 0
        @test result.shortwave_projection_band_count > 0
        @test isfinite(result.spectral_projection_objective_scores.longwave.loss)
        @test isfinite(result.spectral_projection_objective_scores.shortwave.loss)
        @test result.spectral_projection_objective_scores.longwave.loss >= 0
        @test result.spectral_projection_objective_scores.shortwave.loss >= 0
        @test occursin("spectral projection score", result.remaining_gap)

        markdown = EcckdCandidateTransferSmokeValidation.markdown_transfer_smoke(result)
        @test occursin("ecCKD Candidate Transfer Smoke", markdown)
        @test occursin("Transfer Metrics", markdown)
        @test occursin("original-objective broadband loss", markdown)
        @test occursin("spectral-projection original-objective loss", markdown)
    end
end
