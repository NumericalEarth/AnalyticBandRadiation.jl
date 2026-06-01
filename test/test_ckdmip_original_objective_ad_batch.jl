using Test
using JSON

if !isdefined(@__MODULE__, :write_lw_flux_fixture)
    include(joinpath(@__DIR__, "test_ckdmip_original_objective_dataset.jl"))
end
include(joinpath(@__DIR__, "..", "validation", "ckdmip_original_objective_ad_batch.jl"))

@testset "CKDMIP original objective optimizer batch" begin
    mktempdir() do root
        lw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_lw_fluxes_rel-415.h5"))
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_lw_flux_fixture(lw_path)
        write_sw_flux_fixture(sw_path)

        lw_sample = compact_training_sample(read_ckdmip_training_sample(lw_path))
        sw_sample = compact_training_sample(read_ckdmip_training_sample(sw_path))
        parameters = deterministic_flux_perturbation(4prod(size(lw_sample.flux_dn_true)))
        f = p -> ckdmip_flux_correction_batch_loss(p, lw_sample, sw_sample)
        initial = f(parameters)
        gradient = finite_difference_gradient(f, parameters)
        step = loss_reducing_step(f, parameters, gradient)
        @test initial > 0
        @test norm(gradient) > 0
        @test step.accepted
        @test step.final_loss < initial

        result = run_ckdmip_original_objective_ad_batch(; root)
        @test result.status == "optimizer_batch_ready"
        @test result.parameter_count == length(parameters)
        @test result.accepted_step
        @test result.final_loss < result.initial_loss
        @test result.loss_reduction_factor > 1
        @test result.finite_difference_relative_error < 1.0e-5

        json = JSON.parse(json_string(result))
        @test json["case"] == "ckdmip_original_objective_ad_batch"
        @test json["status"] == "optimizer_batch_ready"
        md = markdown_ad_batch(result)
        @test occursin("CKDMIP Original Objective Optimizer Batch", md)
        @test occursin("optimizer-readiness", md)
    end
end
