using Test
using JSON
using NCDatasets

include(joinpath(@__DIR__, "..", "validation", "ckdmip_original_objective_dataset.jl"))

function write_lw_flux_fixture(path)
    mkpath(dirname(path))
    dataset = NCDataset(path, "c")
    try
        defDim(dataset, "column", 1)
        defDim(dataset, "gas", 1)
        defDim(dataset, "level", 2)
        defDim(dataset, "half_level", 3)
        defDim(dataset, "band_lw", 2)
        dataset.attrib["scenario"] = "rel-415"
        dataset.attrib["constituent_id"] = "fixture"
        defVar(dataset, "pressure_hl", Float64, ("half_level", "column"))[:, :] =
            reshape([100.0, 500.0, 1000.0], 3, 1)
        defVar(dataset, "temperature_hl", Float64, ("half_level", "column"))[:, :] =
            reshape([220.0, 250.0, 280.0], 3, 1)
        defVar(dataset, "reference_surface_mole_fraction", Float64, ("gas",))[:] = [1.0]
        defVar(dataset, "mole_fraction_fl", Float64, ("level", "gas", "column"))[:, :, :] =
            reshape([1.0, 1.0], 2, 1, 1)
        defVar(dataset, "band_wavenumber1_lw", Float64, ("band_lw",))[:] = [10.0, 20.0]
        defVar(dataset, "band_wavenumber2_lw", Float64, ("band_lw",))[:] = [20.0, 30.0]
        flux_dn = zeros(Float64, 2, 3, 1)
        flux_up = zeros(Float64, 2, 3, 1)
        flux_dn[:, :, 1] = [1.0 2.0 4.0; 3.0 5.0 8.0]
        flux_up[:, :, 1] = [8.0 6.0 4.0; 7.0 6.0 5.0]
        defVar(dataset, "band_flux_dn_lw", Float64, ("band_lw", "half_level", "column"))[:, :, :] = flux_dn
        defVar(dataset, "band_flux_up_lw", Float64, ("band_lw", "half_level", "column"))[:, :, :] = flux_up
        defVar(dataset, "flux_dn_lw", Float64, ("half_level", "column"))[:, :] =
            reshape([4.0, 7.0, 12.0], 3, 1)
        defVar(dataset, "flux_up_lw", Float64, ("half_level", "column"))[:, :] =
            reshape([15.0, 12.0, 9.0], 3, 1)
    finally
        close(dataset)
    end
end

function write_sw_flux_fixture(path)
    mkpath(dirname(path))
    dataset = NCDataset(path, "c")
    try
        defDim(dataset, "column", 1)
        defDim(dataset, "mu0", 1)
        defDim(dataset, "gas", 1)
        defDim(dataset, "level", 2)
        defDim(dataset, "half_level", 3)
        defDim(dataset, "band_sw", 2)
        dataset.attrib["scenario"] = "rel-415"
        dataset.attrib["constituent_id"] = "fixture"
        defVar(dataset, "pressure_hl", Float64, ("half_level", "column"))[:, :] =
            reshape([100.0, 500.0, 1000.0], 3, 1)
        defVar(dataset, "temperature_hl", Float64, ("half_level", "column"))[:, :] =
            reshape([220.0, 250.0, 280.0], 3, 1)
        defVar(dataset, "mu0", Float64, ("mu0",))[:] = [0.5]
        defVar(dataset, "reference_surface_mole_fraction", Float64, ("gas",))[:] = [1.0]
        defVar(dataset, "mole_fraction_fl", Float64, ("level", "gas", "column"))[:, :, :] =
            reshape([1.0, 1.0], 2, 1, 1)
        defVar(dataset, "band_wavenumber1_sw", Float64, ("band_sw",))[:] = [1000.0, 2000.0]
        defVar(dataset, "band_wavenumber2_sw", Float64, ("band_sw",))[:] = [2000.0, 3000.0]
        flux_dn = zeros(Float64, 2, 3, 1, 1)
        flux_up = zeros(Float64, 2, 3, 1, 1)
        flux_dn[:, :, 1, 1] = [10.0 8.0 5.0; 20.0 17.0 11.0]
        flux_up[:, :, 1, 1] = [1.0 2.0 3.0; 3.0 4.0 5.0]
        defVar(dataset, "band_flux_dn_sw", Float64, ("band_sw", "half_level", "mu0", "column"))[:, :, :, :] = flux_dn
        defVar(dataset, "band_flux_up_sw", Float64, ("band_sw", "half_level", "mu0", "column"))[:, :, :, :] = flux_up
        defVar(dataset, "band_flux_dn_direct_sw", Float64, ("band_sw", "half_level", "mu0", "column"))[:, :, :, :] = flux_dn
        defVar(dataset, "flux_dn_sw", Float64, ("half_level", "mu0", "column"))[:, :, :] =
            reshape([30.0, 25.0, 16.0], 3, 1, 1)
        defVar(dataset, "flux_up_sw", Float64, ("half_level", "mu0", "column"))[:, :, :] =
            reshape([4.0, 6.0, 8.0], 3, 1, 1)
        defVar(dataset, "flux_dn_direct_sw", Float64, ("half_level", "mu0", "column"))[:, :, :] =
            reshape([30.0, 25.0, 16.0], 3, 1, 1)
        defDim(dataset, "wavenumber", 2)
        defVar(dataset, "wavenumber", Float64, ("wavenumber",))[:] = [1000.0, 2000.0]
        defVar(dataset, "spectral_flux_dn_surf_sw", Float64, ("wavenumber", "mu0", "column"))[:, :, :] =
            reshape([5.0, 11.0], 2, 1, 1)
    finally
        close(dataset)
    end
end

@testset "CKDMIP original objective dataset loader" begin
    pressure = [100.0, 500.0, 1000.0]
    weight = ckdmip_layer_weight(pressure)
    @test isapprox(sum(weight), 1.0; atol = 1e-14)
    @test weight ≈ (sqrt.(pressure[2:end]) .- sqrt.(pressure[1:2])) ./
                    sum(sqrt.(pressure[2:end]) .- sqrt.(pressure[1:2]))

    flux_dn = [10.0 20.0; 8.0 17.0; 5.0 11.0]
    flux_up = [1.0 3.0; 2.0 4.0; 3.0 5.0]
    sw_heating = ecckd_flux_heating_rate(pressure, flux_dn)
    lw_heating = ecckd_flux_heating_rate(pressure, flux_dn, flux_up)
    conversion = .-(ECCKD_GRAVITY / ECCKD_SPECIFIC_HEAT_AIR) ./ diff(pressure)
    @test sw_heating[:, 1] ≈ conversion .* [-2.0, -3.0]
    @test lw_heating[:, 1] ≈ conversion .* [-3.0, -4.0]

    mktempdir() do root
        lw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_lw_fluxes_rel-415.h5"))
        sw_path = joinpath(root, expected_flux_path("ckdmip_evaluation1_sw_fluxes_rel-415.h5"))
        write_lw_flux_fixture(lw_path)
        write_sw_flux_fixture(sw_path)

        lw = read_ckdmip_training_sample(lw_path)
        sw = read_ckdmip_training_sample(sw_path)
        @test lw.kind == "longwave"
        @test sw.kind == "shortwave"
        @test size(lw.flux_dn_true) == (3, 2)
        @test size(sw.flux_dn_true) == (3, 2)
        @test size(lw.heating_rate_true) == (2, 2)
        @test size(sw.heating_rate_true) == (2, 2)
        @test objective_self_loss(lw) == 0
        @test objective_self_loss(sw) == 0
        @test !lw.has_spectral_boundary
        @test sw.has_spectral_boundary

        result = run_ckdmip_original_objective_dataset(;
            root,
            training_flux_filenames = [
                "ckdmip_evaluation1_lw_fluxes_rel-415.h5",
                "ckdmip_evaluation1_sw_fluxes_rel-415.h5",
            ],
        )
        @test result.status == "dataset_samples_ready"
        @test result.training_flux_file_count == 2
        @test result.training_flux_schema_ok_count == 2
        @test all(schema -> schema.schema_ok, result.training_flux_schemas)
        @test length(result.samples) == 2
        @test all(sample -> sample.self_loss == 0, result.samples)

        md = markdown_objective_dataset(result)
        @test occursin("CKDMIP Original Objective Dataset", md)
        @test occursin("dataset_samples_ready", md)
        json = json_string(result)
        parsed = JSON.parse(json)
        @test parsed["case"] == "ckdmip_original_objective_dataset"
    end
end
