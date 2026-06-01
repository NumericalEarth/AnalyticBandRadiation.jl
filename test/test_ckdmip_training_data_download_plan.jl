using JSON

module CKDMIPTrainingDataDownloadPlanValidation
include(joinpath(@__DIR__, "..", "validation", "ckdmip_training_data_download_plan.jl"))
end

@testset "CKDMIP training data download plan" begin
    root = normpath(joinpath(@__DIR__, ".."))
    json_path = joinpath(root, "validation", "results", "ckdmip_training_data_download_plan.json")
    md_path = joinpath(root, "validation", "results", "ckdmip_training_data_download_plan.md")

    redirect_stdout(devnull) do
        CKDMIPTrainingDataDownloadPlanValidation.ckdmip_training_data_download_plan_main()
    end

    @test isfile(json_path)
    @test isfile(md_path)
    @test occursin("CKDMIP Training Data Download Plan", read(md_path, String))

    result = JSON.parsefile(json_path)
    @test result["case"] == "ckdmip_training_data_download_plan"
    @test result["status"] == "manual_or_batch_download_required"
    @test result["target_env"] == "RH_CKDMIP_DATA_PATH"
    @test result["task_count"] == length(result["tasks"])
    @test result["task_count"] > 20
    @test "evaluation1/lw_spectra" in result["expected_layout_roots"]
    @test "evaluation2/sw_fluxes" in result["expected_layout_roots"]

    commands = [task["command"] for task in result["tasks"]]
    @test any(command -> occursin("ckdmip_evaluation1_concentrations_present.nc", command), commands)
    @test any(command -> occursin("lw_spectra/evaluation1", command), commands)
    @test any(command -> occursin("sw_spectra/evaluation2", command), commands)
    @test any(command -> occursin("ln -sf", command) &&
                         occursin("mmm/sw_spectra_extras/ckdmip_ssi.h5", command), commands)
    @test all(command -> occursin("\$RH_CKDMIP_DATA_PATH", command), commands)
end
