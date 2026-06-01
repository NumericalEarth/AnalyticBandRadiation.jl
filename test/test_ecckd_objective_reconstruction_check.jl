using JSON

module EcckdObjectiveReconstructionValidation
include(joinpath(@__DIR__, "..", "validation", "ecckd_objective_reconstruction_check.jl"))
end

@testset "ecCKD original-objective reconstruction check" begin
    root = normpath(joinpath(@__DIR__, ".."))
    json_path = joinpath(root, "validation", "results", "ecckd_objective_reconstruction_check.json")
    md_path = joinpath(root, "validation", "results", "ecckd_objective_reconstruction_check.md")
    redirect_stdout(devnull) do
        EcckdObjectiveReconstructionValidation.ecckd_objective_reconstruction_check_main()
    end

    @test isfile(json_path)
    @test isfile(md_path)
    @test occursin("ecCKD Objective Reconstruction Check", read(md_path, String))

    result = JSON.parsefile(json_path)
    @test result["case"] == "ecckd_objective_reconstruction_check"
    @test result["status"] in (
        "blocked_missing_original_training_assets",
        "ready_to_reconstruct_original_objective",
        "passed",
    )
    @test endswith(result["ckdmip_training_data_preflight"],
                   joinpath("validation", "results", "ckdmip_training_data_preflight.json"))
    @test all(item -> item in (
              "original LBL training database",
              "derived ecCKD training flux products",
          ), result["missing_for_exact_original_recovery"])
    @test result["ecckd_source_root"] !== nothing

    checks = Dict(check["name"] => check for check in result["checks"])
    @test checks["published ecCKD CKD-definition files"]["present"] == true
    @test checks["CKDMIP evaluation profiles and reference fluxes"]["present"] == true
    @test haskey(checks, "original LBL training database")
    @test haskey(checks, "derived ecCKD training flux products")
    @test checks["official ecCKD generator and optimizer source"]["present"] == true
    @test checks["official ecCKD objective weights and training scripts"]["present"] == true

    if !isempty(result["blockers"])
        blocker_text = join(result["blockers"], "\n")
        @test occursin("RH_CKDMIP_DATA_PATH", blocker_text) ||
              occursin("derived ecCKD", blocker_text)
    end
end
