using JSON

module CKDMIPTrainingDataPreflightValidation
include(joinpath(@__DIR__, "..", "validation", "ckdmip_training_data_preflight.jl"))
end

@testset "CKDMIP training data preflight" begin
    root = normpath(joinpath(@__DIR__, ".."))
    json_path = joinpath(root, "validation", "results", "ckdmip_training_data_preflight.json")
    md_path = joinpath(root, "validation", "results", "ckdmip_training_data_preflight.md")

    result = withenv("RH_CKDMIP_DATA_PATH" => "") do
        CKDMIPTrainingDataPreflightValidation.run_ckdmip_training_data_preflight()
    end
    @test result.case == "ckdmip_training_data_preflight"
    @test result.status == "missing_ckdmip_data_root"
    @test result.ckdmip_data_root === nothing
    @test result.expected_training_flux_file_count == 52
    @test result.upstream_training_flux_file_count == 34
    @test result.derived_training_flux_file_count == 18
    @test "ckdmip_evaluation1_lw_fluxes_present.h5" in result.expected_training_flux_files
    @test "ckdmip_evaluation1_lw_fluxes_present.h5" in result.upstream_training_flux_files
    @test "ckdmip_evaluation1_sw_fluxes_rel-415.h5" in result.expected_training_flux_files
    @test "ckdmip_evaluation1_sw_fluxes_rel-415.h5" in result.derived_training_flux_files
    @test !("ckdmip_evaluation2_lw_fluxes_rel-415.h5" in result.expected_training_flux_files)
    @test !("ckdmip_evaluation2_sw_fluxes_rel-415.h5" in result.expected_training_flux_files)
    @test any(blocker -> occursin("RH_CKDMIP_DATA_PATH", blocker), result.blockers)

    markdown = CKDMIPTrainingDataPreflightValidation.markdown_preflight(result)
    @test occursin("CKDMIP Training Data Preflight", markdown)

    mktempdir() do ckdmip_root
        for entry in CKDMIPTrainingDataPreflightValidation.expected_layout_roots()
            mkpath(joinpath(ckdmip_root, entry))
        end
        for entry in (
            "mmm/lw_spectra",
            "mmm/sw_spectra",
            "idealized/lw_spectra",
            "idealized/sw_spectra",
            "evaluation1/lw_spectra",
            "evaluation1/sw_spectra",
            "evaluation2/lw_spectra",
            "evaluation2/sw_spectra",
        )
            write(joinpath(ckdmip_root, entry, "marker.dat"), "fixture\n")
        end

        manifest = CKDMIPTrainingDataPreflightValidation.run_ecckd_published_training_manifest()
        upstream = [
            file for file in CKDMIPTrainingDataPreflightValidation.expected_training_flux_files(manifest)
            if !CKDMIPTrainingDataPreflightValidation.derived_training_flux_file(file)
        ]
        upstream_paths = vcat(
            CKDMIPTrainingDataPreflightValidation.key_files(),
            CKDMIPTrainingDataPreflightValidation.expected_flux_path.(upstream),
        )
        for path in upstream_paths
            full_path = joinpath(ckdmip_root, path)
            mkpath(dirname(full_path))
            write(full_path, "fixture\n")
        end

        derived_result = withenv("RH_CKDMIP_DATA_PATH" => ckdmip_root) do
            CKDMIPTrainingDataPreflightValidation.run_ckdmip_training_data_preflight()
        end
        @test derived_result.status == "ready_for_derived_flux_generation"
        @test isempty(derived_result.blockers)
        @test length(derived_result.derived_flux_generation_blockers) == 18
        @test all(blocker -> occursin("derived ecCKD training flux product", blocker),
                  derived_result.derived_flux_generation_blockers)

        derived = [
            file for file in CKDMIPTrainingDataPreflightValidation.expected_training_flux_files(manifest)
            if CKDMIPTrainingDataPreflightValidation.derived_training_flux_file(file)
        ]
        for path in CKDMIPTrainingDataPreflightValidation.expected_flux_path.(derived)
            full_path = joinpath(ckdmip_root, path)
            mkpath(dirname(full_path))
            write(full_path, "fixture\n")
        end

        ready_result = withenv("RH_CKDMIP_DATA_PATH" => ckdmip_root) do
            CKDMIPTrainingDataPreflightValidation.run_ckdmip_training_data_preflight()
        end
        @test ready_result.status == "ready_for_original_ecckd_objective"
        @test isempty(ready_result.blockers)
        @test isempty(ready_result.derived_flux_generation_blockers)
    end
end
