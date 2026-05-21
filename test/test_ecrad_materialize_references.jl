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
end
