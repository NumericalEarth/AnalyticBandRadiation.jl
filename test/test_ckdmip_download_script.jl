@testset "CKDMIP download helper dry run" begin
    root = normpath(joinpath(@__DIR__, ".."))
    script = joinpath(root, "validation", "download_ckdmip_training_data.sh")

    @test success(`bash -n $script`)

    output = read(setenv(`bash $script`,
                         "RH_CKDMIP_DATA_PATH" => "/tmp/ckdmip-dryrun",
                         "CKDMIP_DRY_RUN" => "true"), String)
    @test occursin("DRY-RUN:", output)
    @test occursin("ckdmip_mmm_concentrations.nc", output)
    @test occursin("lw_spectra/evaluation1", output)
    @test occursin("sw_fluxes/evaluation2", output)
    @test occursin("ckdmip_training_data_preflight.jl", output)

    preflight_output = read(setenv(`bash $script`,
                                   "RH_CKDMIP_DATA_PATH" => "/tmp/ckdmip-dryrun",
                                   "CKDMIP_DRY_RUN" => "true",
                                   "CKDMIP_RUN_PREFLIGHT" => "true"), String)
    @test occursin("DRY-RUN:", preflight_output)
    @test occursin("--project=$(joinpath(root, "test"))", preflight_output)
    @test occursin(joinpath(root, "validation", "ckdmip_training_data_preflight.jl"),
                   preflight_output)

    err = Pipe()
    proc = run(pipeline(setenv(`bash $script`, "CKDMIP_DRY_RUN" => "true"); stderr = err);
               wait = false)
    wait(proc)
    close(err.in)
    @test proc.exitcode == 2
    @test occursin("Set RH_CKDMIP_DATA_PATH", read(err, String))
end
