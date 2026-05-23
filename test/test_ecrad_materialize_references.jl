using NCDatasets

include(joinpath(@__DIR__, "..", "validation", "materialize_ecrad_references.jl"))

@testset "ecRad reference materializer conventions" begin
    p_interface = [100.0 200.0;
                   300.0 600.0;
                   900.0 1000.0]
    t_interface = [200.0 220.0;
                   260.0 280.0;
                   300.0 320.0]

    expected = [
        (200.0 * 100.0 + 260.0 * 300.0) / (100.0 + 300.0)  (220.0 * 200.0 + 280.0 * 600.0) / (200.0 + 600.0);
        (260.0 * 300.0 + 300.0 * 900.0) / (300.0 + 900.0)  (280.0 * 600.0 + 320.0 * 1000.0) / (600.0 + 1000.0)
    ]

    @test temperature_layer_from_interfaces(p_interface, t_interface) ≈ expected

    reference_path = joinpath(@__DIR__, "..", "validation", "reference", "ecrad",
                              "clear_sky_tropical_column.nc")
    if isfile(reference_path)
        NCDatasets.NCDataset(reference_path) do dataset
            @test haskey(dataset, "cos_solar_zenith_angle")
            @test haskey(dataset, "solar_irradiance")
            @test haskey(dataset, "o3")
            @test haskey(dataset, "ch4")
            @test haskey(dataset, "n2o")
            @test all(Array(dataset["cos_solar_zenith_angle"]) .>= 0)
            @test all(Array(dataset["solar_irradiance"]) .> 0)
        end
    end

    ecckd_reference_path = joinpath(@__DIR__, "..", "validation", "reference", "ecrad",
                                    "ecckd_clear_sky_tropical_column.nc")
    if isfile(ecckd_reference_path)
        NCDatasets.NCDataset(ecckd_reference_path) do dataset
            @test dataset.attrib["source_flux_scope"] == "total"
            @test haskey(dataset, "o3")
            @test haskey(dataset, "ch4")
            @test haskey(dataset, "n2o")
            @test haskey(dataset, "surface_longwave_up_spectral")
            @test haskey(dataset, "surface_albedo_spectral")
            @test haskey(dataset, "surface_albedo_direct_spectral")
        end
    end

    matched_64x64_reference_path =
        joinpath(@__DIR__, "..", "validation", "reference", "ecrad",
                 "ecckd_64x64_clear_sky_tropical_column.nc")
    if isfile(matched_64x64_reference_path)
        NCDatasets.NCDataset(matched_64x64_reference_path) do dataset
            @test haskey(dataset, "surface_longwave_up_spectral")
            @test haskey(dataset, "surface_albedo_spectral")
            @test size(dataset["surface_longwave_up_spectral"], 1) == 64
            @test size(dataset["surface_albedo_spectral"], 1) == 64
        end
    end

    matched_64x96_reference_path =
        joinpath(@__DIR__, "..", "validation", "reference", "ecrad",
                 "ecckd_64x96_clear_sky_tropical_column.nc")
    if isfile(matched_64x96_reference_path)
        NCDatasets.NCDataset(matched_64x96_reference_path) do dataset
            @test haskey(dataset, "surface_longwave_up_spectral")
            @test haskey(dataset, "surface_albedo_spectral")
            @test size(dataset["surface_longwave_up_spectral"], 1) == 64
            @test size(dataset["surface_albedo_spectral"], 1) == 96
        end
    end

    for (case_name, ng_lw, ng_sw) in (
        ("ecckd_32x64_clear_sky_tropical_column.nc", 32, 64),
        ("ecckd_32x96_clear_sky_tropical_column.nc", 32, 96),
        ("ecckd_64x32_clear_sky_tropical_column.nc", 64, 32),
    )
        matched_reference_path =
            joinpath(@__DIR__, "..", "validation", "reference", "ecrad", case_name)
        if isfile(matched_reference_path)
            NCDatasets.NCDataset(matched_reference_path) do dataset
                @test haskey(dataset, "surface_longwave_up_spectral")
                @test haskey(dataset, "surface_albedo_spectral")
                @test size(dataset["surface_longwave_up_spectral"], 1) == ng_lw
                @test size(dataset["surface_albedo_spectral"], 1) == ng_sw
            end
        end
    end

    all_sky_reference_path = joinpath(@__DIR__, "..", "validation", "reference", "ecrad",
                                      "all_sky_tropical_column.nc")
    if isfile(all_sky_reference_path)
        NCDatasets.NCDataset(all_sky_reference_path) do dataset
            @test haskey(dataset, "re_liquid")
            @test haskey(dataset, "re_ice")
            @test haskey(dataset, "inv_cloud_effective_size")
            @test haskey(dataset, "fractional_std")
            @test haskey(dataset, "overlap_param")
            @test haskey(dataset, "aerosol_mmr")
            @test haskey(dataset, "surface_albedo_direct_spectral")
            @test haskey(dataset, "surface_albedo_spectral")
            @test haskey(dataset, "surface_albedo_direct_gpoint")
            @test dimnames(dataset["overlap_param"]) == ("overlap_interface", "column")
            @test dimnames(dataset["aerosol_mmr"]) == ("layer", "aerosol_type", "column")
            @test dimnames(dataset["surface_albedo_direct_gpoint"]) ==
                  ("shortwave_gpoint", "column")
        end
    end

    ecckd_all_sky_reference_path = joinpath(@__DIR__, "..", "validation", "reference", "ecrad",
                                            "ecckd_all_sky_tropical_column.nc")
    if isfile(ecckd_all_sky_reference_path)
        NCDatasets.NCDataset(ecckd_all_sky_reference_path) do dataset
            @test dataset.attrib["source_flux_scope"] == "total"
            @test occursin("ecrad_meridian_ecckd_tc_out_REFERENCE.nc",
                           dataset.attrib["source_output"])
            @test haskey(dataset, "o3")
            @test haskey(dataset, "ch4")
            @test haskey(dataset, "n2o")
            @test haskey(dataset, "surface_longwave_up_spectral")
            @test haskey(dataset, "surface_albedo_direct_gpoint")
            @test haskey(dataset, "lw_up_clear")
            @test haskey(dataset, "sw_down_clear")
        end
    end

    for (case_name, ng_lw, ng_sw) in (
        ("ecckd_32x64_all_sky_tropical_column.nc", 32, 64),
        ("ecckd_32x96_all_sky_tropical_column.nc", 32, 96),
        ("ecckd_64x32_all_sky_tropical_column.nc", 64, 32),
        ("ecckd_64x64_all_sky_tropical_column.nc", 64, 64),
        ("ecckd_64x96_all_sky_tropical_column.nc", 64, 96),
    )
        matched_all_sky_path =
            joinpath(@__DIR__, "..", "validation", "reference", "ecrad", case_name)
        if isfile(matched_all_sky_path)
            NCDatasets.NCDataset(matched_all_sky_path) do dataset
                @test dataset.attrib["source_flux_scope"] == "total"
                @test haskey(dataset, "cloud_fraction")
                @test haskey(dataset, "aerosol_mmr")
                @test haskey(dataset, "lw_up_clear")
                @test haskey(dataset, "sw_down_clear")
                @test haskey(dataset, "surface_longwave_up_spectral")
                @test haskey(dataset, "surface_albedo_spectral")
                @test size(dataset["surface_longwave_up_spectral"], 1) == ng_lw
                @test size(dataset["surface_albedo_spectral"], 1) == ng_sw
                if case_name == "ecckd_32x64_all_sky_tropical_column.nc"
                    @test haskey(dataset, "surface_albedo_direct_gpoint")
                    @test haskey(dataset, "toa_shortwave_down_spectral")
                    @test size(dataset["surface_albedo_direct_gpoint"], 1) == ng_sw
                    @test occursin("radiative_properties_ecckd_32x64.nc",
                                   dataset.attrib["saved_shortwave_boundary_source"])
                    @test occursin("surface_albedo_direct_gpoint",
                                   dataset.attrib["saved_shortwave_boundaries"])
                end
            end
        end
    end
end
