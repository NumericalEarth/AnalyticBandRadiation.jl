using JSON

if !isdefined(@__MODULE__, :write_lw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateOriginalObjectiveScoreValidation
include(joinpath(@__DIR__, "..", "validation", "ecckd_candidate_original_objective_score.jl"))
end

@testset "ecCKD candidate original-objective score" begin
    mktempdir() do root
        lw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_lw_fluxes_rel-415.h5"))
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_lw_flux_fixture(lw_path)
        write_sw_flux_fixture(sw_path)
        candidate =
            EcckdCandidateOriginalObjectiveScoreValidation.official_ecckd_definition_path(
                :shortwave_32,
            )

        result =
            EcckdCandidateOriginalObjectiveScoreValidation.run_ecckd_candidate_original_objective_score(;
                root,
                candidate_path = candidate,
            )
        @test result.status == "candidate_objective_score_ready"
        @test result.candidate_metrics.status == "passed"
        @test length(result.sample_scores) == 2
        @test all(row -> row.zero_forward_loss == 0.0, result.sample_scores)
        @test all(row -> row.perturbed_forward_loss > 0.0, result.sample_scores)
        @test all(row -> row.loss_increases_under_perturbation, result.sample_scores)
        @test occursin("differentiable transfer", result.remaining_gap)

        json = JSON.parse(EcckdCandidateOriginalObjectiveScoreValidation.json_string(result))
        @test json["case"] == "ecckd_candidate_original_objective_score"
        @test json["status"] == "candidate_objective_score_ready"

        markdown =
            EcckdCandidateOriginalObjectiveScoreValidation.markdown_candidate_objective_score(
                result,
            )
        @test occursin("ecCKD Candidate Original-Objective Score", markdown)
        @test occursin("Sample Objective Scores", markdown)
        @test occursin("Perturbed-forward loss", markdown)
    end
end
