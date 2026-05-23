if !isdefined(@__MODULE__, :write_sw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end

module EcckdCandidateTableHumidityTribinProbeValidation
include(joinpath(@__DIR__, "..", "validation",
                 "ecckd_candidate_table_humidity_tribin_probe.jl"))
end

@testset "ecCKD candidate table humidity-tribin probe" begin
    mktempdir() do root
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_sw_flux_fixture(sw_path)
        reference =
            EcckdCandidateTableHumidityTribinProbeValidation.official_ecckd_definition_path(
                :shortwave_32,
            )
        candidate = joinpath(root, "humidity_tribin_candidate.nc")

        result =
            EcckdCandidateTableHumidityTribinProbeValidation.run_table_humidity_tribin_probe(;
                root,
                reference_path = reference,
                candidate_path = candidate,
                scenarios = ("rel-415",),
                initial_parameters = [-0.1, -0.1, -0.1, 0.0, -0.05, 0.0],
                deltas = (-0.01, 0.01),
            )
        @test result.status in (
            "humidity_tribin_probe_improved",
            "humidity_tribin_probe_no_descent",
        )
        @test result.present_count == 1
        @test result.scenario_count == 1
        @test result.tested_move_count == 6
        @test isfinite(result.seed_aggregate_loss)
        @test isfinite(result.best_aggregate_loss)
        @test result.aggregate_loss_reduction_factor >= 0
        @test result.worst_loss_ratio_reduction_factor >= 0
        @test length(result.best_aggregate_rows) == 1
        @test length(result.best_minimax_rows) == 1
        @test isfile(candidate)
        @test occursin("not published-model recovery", result.interpretation)

        markdown =
            EcckdCandidateTableHumidityTribinProbeValidation.markdown_table_humidity_tribin_probe(
                result,
            )
        @test occursin("ecCKD Candidate Table Humidity-Tribin Probe", markdown)
        @test occursin("Probe Summary", markdown)
    end
end
