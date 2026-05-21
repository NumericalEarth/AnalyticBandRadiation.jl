using NCDatasets

module CKDMIPFluxChunkConcatenationValidation
include(joinpath(@__DIR__, "..", "validation", "concat_ckdmip_flux_chunks.jl"))
end

@testset "ecCKD derived flux launcher dry run" begin
    root = normpath(joinpath(@__DIR__, ".."))
    script = joinpath(root, "validation", "generate_ecckd_derived_fluxes.sh")
    install_script = joinpath(root, "validation", "install_ecckd_derived_fluxes.sh")

    @test success(`bash -n $script`)
    @test success(`bash -n $install_script`)

    mktempdir() do workdir
        output = read(setenv(`bash $script`,
                             "RH_CKDMIP_DATA_PATH" => "/tmp/ckdmip-dryrun",
                             "RH_ECCKD_LBL_WORKDIR" => workdir,
                             "RH_ECCKD_DERIVED_FLUX_DRY_RUN" => "true"), String)
        @test occursin("Dry run only", output)
        @test occursin("run_lw_lbl_evaluation.sh", output)
        @test occursin("run_sw_lbl_evaluation.sh", output)
        @test occursin("install generated ckdmip_evaluation1_*_fluxes_{5gas,rel}-*.h5 files", output)

        lw_script = joinpath(workdir, "ecckd", "test", "run_lw_lbl_evaluation.sh")
        sw_script = joinpath(workdir, "ecckd", "test", "run_sw_lbl_evaluation.sh")
        config = joinpath(workdir, "ecckd", "test", "config.h")
        @test isfile(lw_script)
        @test isfile(sw_script)
        @test isfile(config)
        @test occursin("SCENARIOS=\"5gas-180 5gas-280 5gas-415 5gas-560 5gas-1120 5gas-2240 rel-180 rel-280 rel-415 rel-560 rel-1120 rel-2240\"",
                       read(lw_script, String))
        @test occursin("SCENARIOS=\"rel-180 rel-280 rel-415 rel-560 rel-1120 rel-2240\"",
                       read(sw_script, String))
        @test occursin("concat_ckdmip_flux_chunks.jl", read(lw_script, String))
        @test occursin("concat_ckdmip_flux_chunks.jl", read(sw_script, String))
        @test occursin("*** REUSING \$OUTFILE ***", read(lw_script, String))
        @test occursin("*** REUSING \$OUTFILE ***", read(sw_script, String))
        @test occursin("CKDMIP_DATA_DIR=/tmp/ckdmip-dryrun", read(config, String))
        @test occursin("WORK_DIR=$(workdir)/work", read(config, String))

        launcher = read(script, String)
        @test occursin("install_ecckd_derived_fluxes.sh", launcher)

        installer = read(install_script, String)
        @test occursin("evaluation1/lw_fluxes", installer)
        @test occursin("evaluation1/sw_fluxes", installer)
        @test occursin("ckdmip_evaluation1_lw_fluxes_5gas-*.h5", installer)
        @test occursin("ckdmip_evaluation1_lw_fluxes_rel-*.h5", installer)
        @test occursin("ckdmip_evaluation1_sw_fluxes_rel-*.h5", installer)
    end

    err = Pipe()
    proc = run(pipeline(setenv(`bash $script`,
                               "RH_ECCKD_DERIVED_FLUX_DRY_RUN" => "true");
                        stderr = err);
               wait = false)
    wait(proc)
    close(err.in)
    @test proc.exitcode == 2
    @test occursin("Set RH_CKDMIP_DATA_PATH", read(err, String))

    mktempdir() do workdir
        ckdmip_root = joinpath(workdir, "ckdmip")
        lbl_root = joinpath(workdir, "lbl")
        lw_source = joinpath(lbl_root, "work", "lw_lbl_fluxes")
        sw_source = joinpath(lbl_root, "work", "sw_lbl_fluxes")
        mkpath(lw_source)
        mkpath(sw_source)
        write(joinpath(lw_source, "ckdmip_evaluation1_lw_fluxes_5gas-180.h5"), "lw 5gas\n")
        write(joinpath(lw_source, "ckdmip_evaluation1_lw_fluxes_rel-180.h5"), "lw rel\n")
        write(joinpath(sw_source, "ckdmip_evaluation1_sw_fluxes_rel-180.h5"), "sw rel\n")

        output = read(setenv(`bash $install_script`,
                             "RH_CKDMIP_DATA_PATH" => ckdmip_root,
                             "RH_ECCKD_LBL_WORKDIR" => lbl_root), String)
        @test occursin("Installed 3 derived ecCKD flux product", output)
        @test read(joinpath(ckdmip_root, "evaluation1", "lw_fluxes",
                            "ckdmip_evaluation1_lw_fluxes_5gas-180.h5"), String) == "lw 5gas\n"
        @test read(joinpath(ckdmip_root, "evaluation1", "lw_fluxes",
                            "ckdmip_evaluation1_lw_fluxes_rel-180.h5"), String) == "lw rel\n"
        @test read(joinpath(ckdmip_root, "evaluation1", "sw_fluxes",
                            "ckdmip_evaluation1_sw_fluxes_rel-180.h5"), String) == "sw rel\n"
    end
end

@testset "CKDMIP flux chunk concatenation" begin
    mktempdir() do dir
        inputs = String[]
        for i in 1:3
            path = joinpath(dir, "chunk$(i).h5")
            ds = NCDataset(path, "c")
            defDim(ds, "column", Inf)
            defDim(ds, "half_level", 2)
            defDim(ds, "gas", 1)
            ds.attrib["scenario"] = "fixture"
            v = defVar(ds, "flux_up_lw", Float32, ("half_level", "column"))
            v.attrib["units"] = "W m-2"
            v[:, 1:2] = fill(Float32(i), 2, 2)
            g = defVar(ds, "reference_surface_mole_fraction", Float32, ("gas",);
                       attrib = Dict("units" => "1"))
            g[:] = Float32[415e-6]
            close(ds)
            push!(inputs, path)
        end

        output = joinpath(dir, "joined.h5")
        CKDMIPFluxChunkConcatenationValidation.concat_ckdmip_flux_chunks(inputs, output)

        ds = NCDataset(output)
        try
            @test ds.dim["column"] == 6
            @test ds.dim["half_level"] == 2
            @test ds.attrib["scenario"] == "fixture"
            @test occursin("concat_ckdmip_flux_chunks", ds.attrib["history"])
            @test ds["flux_up_lw"].attrib["units"] == "W m-2"
            @test ds["flux_up_lw"][:, :] == Float32[
                1 1 2 2 3 3
                1 1 2 2 3 3
            ]
            @test ds["reference_surface_mole_fraction"][:] == Float32[415e-6]
        finally
            close(ds)
        end
    end
end
