using NCDatasets

@testset "cloud scattering table" begin
    table = CloudScatteringTable(
        medium = "liquid-water",
        particle_type = "droplet",
        wavenumber = [100.0, 200.0],
        effective_radius = [1.0e-6, 3.0e-6],
        mass_extinction_coefficient = [10.0 30.0; 20.0 40.0],
        single_scattering_albedo = [0.8 1.0; 0.6 0.9],
        asymmetry_factor = [0.7 0.9; 0.5 0.8],
    )

    @test eltype(table) == Float64
    @test table.medium == "liquid-water"
    props = cloud_scattering_properties(table, 1, 2.0e-6)
    @test props.mass_extinction_coefficient ≈ 20.0
    @test props.single_scattering_albedo ≈ 0.9
    @test props.asymmetry_factor ≈ 0.8

    low = cloud_scattering_properties(table, 2, 0.1e-6)
    high = cloud_scattering_properties(table, 2, 10.0e-6)
    @test low.mass_extinction_coefficient ≈ 20.0
    @test high.mass_extinction_coefficient ≈ 40.0
end

@testset "cloud scattering g-point mapping" begin
    table = CloudScatteringTable(
        medium = "liquid-water",
        particle_type = "droplet",
        wavenumber = [100.0, 200.0, 300.0],
        effective_radius = [1.0e-6, 3.0e-6],
        mass_extinction_coefficient = [10.0 20.0; 30.0 40.0; 50.0 60.0],
        single_scattering_albedo = fill(0.8, 3, 2),
        asymmetry_factor = fill(0.5, 3, 2),
    )
    mapping = EcCKDSpectralMapping(
        wavenumber1 = [90.0, 190.0, 290.0],
        wavenumber2 = [110.0, 210.0, 310.0],
        gpoint_fraction = [1.0 0.0; 0.5 0.5; 0.0 1.0],
    )

    props = cloud_scattering_gpoint_properties(table, mapping, 1.0e-6)
    @test length(props.mass_extinction_coefficient) == 2
    @test props.mass_extinction_coefficient[1] ≈ (10.0 + 0.5 * 30.0) / 1.5
    @test props.mass_extinction_coefficient[2] ≈ (0.5 * 30.0 + 50.0) / 1.5
    @test all(props.single_scattering_albedo .≈ 0.8)
    @test all(props.asymmetry_factor .≈ 0.5)

    thick_table = CloudScatteringTable(
        medium = "liquid-water",
        particle_type = "droplet",
        wavenumber = [100.0, 200.0, 300.0],
        effective_radius = [1.0e-6, 3.0e-6],
        mass_extinction_coefficient = [10.0 20.0; 30.0 40.0; 50.0 60.0],
        single_scattering_albedo = [0.95 0.9; 0.7 0.8; 0.4 0.5],
        asymmetry_factor = [0.85 0.8; 0.6 0.65; 0.3 0.4],
    )
    thick = cloud_scattering_gpoint_properties(
        thick_table, mapping, 1.0e-6;
        mapping_method = :ecrad,
        delta_eddington_average = true,
        thick_averaging = true,
    )
    thin = cloud_scattering_gpoint_properties(
        thick_table, mapping, 1.0e-6;
        mapping_method = :ecrad,
        delta_eddington_average = true,
        thick_averaging = false,
    )
    @test length(thick.single_scattering_albedo) == 2
    @test all(0 .<= thick.single_scattering_albedo .<= 1)
    @test thick.single_scattering_albedo != thin.single_scattering_albedo
end

@testset "NCDatasets cloud scattering table reader" begin
    path = tempname() * ".nc"
    NCDataset(path, "c") do ds
        defDim(ds, "wavenumber", 2)
        defDim(ds, "effective_radius", 3)
        ds.attrib["medium"] = "ice"
        ds.attrib["particle_type"] = "cloud-ice"
        wavenumber = defVar(ds, "wavenumber", Float64, ("wavenumber",))
        radius = defVar(ds, "effective_radius", Float64, ("effective_radius",))
        mass_ext = defVar(ds, "mass_extinction_coefficient", Float64,
                          ("wavenumber", "effective_radius"))
        ssa = defVar(ds, "single_scattering_albedo", Float64,
                     ("wavenumber", "effective_radius"))
        asymmetry = defVar(ds, "asymmetry_factor", Float64,
                           ("wavenumber", "effective_radius"))
        wavenumber[:] = [100.0, 200.0]
        radius[:] = [1.0e-6, 2.0e-6, 3.0e-6]
        mass_ext[:, :] = [10.0 20.0 30.0; 40.0 50.0 60.0]
        ssa[:, :] = fill(0.9, 2, 3)
        asymmetry[:, :] = fill(0.7, 2, 3)
    end

    table = read_cloud_scattering_table(path)
    @test table isa CloudScatteringTable
    @test table.medium == "ice"
    @test table.particle_type == "cloud-ice"
    @test table.wavenumber == [100.0, 200.0]
    @test table.effective_radius == [1.0e-6, 2.0e-6, 3.0e-6]
    @test table.mass_extinction_coefficient[2, 3] == 60.0
end

@testset "NCDatasets ecCKD spectral mapping reader" begin
    path = tempname() * ".nc"
    NCDataset(path, "c") do ds
        defDim(ds, "wavenumber", 2)
        defDim(ds, "g_point", 3)
        w1 = defVar(ds, "wavenumber1", Float64, ("wavenumber",))
        w2 = defVar(ds, "wavenumber2", Float64, ("wavenumber",))
        frac = defVar(ds, "gpoint_fraction", Float64, ("wavenumber", "g_point"))
        w1[:] = [100.0, 200.0]
        w2[:] = [150.0, 250.0]
        frac[:, :] = [1.0 0.0 0.0; 0.0 0.25 0.75]
    end

    mapping = read_ecckd_spectral_mapping(path)
    @test mapping isa EcCKDSpectralMapping
    @test mapping.wavenumber1 == [100.0, 200.0]
    @test mapping.wavenumber2 == [150.0, 250.0]
    @test size(mapping.gpoint_fraction) == (2, 3)
    @test mapping.interval_weight == [1.0, 1.0]
end

@testset "official ecRad cloud scattering files" begin
    root = normpath(joinpath(@__DIR__, ".."))
    liquid_path = joinpath(root, "validation", "external", "ecrad", "data",
                           "mie_droplet_scattering.nc")
    ice_path = joinpath(root, "validation", "external", "ecrad", "data",
                        "baum-general-habit-mixture_ice_scattering.nc")

    if isfile(liquid_path) && isfile(ice_path)
        liquid = read_cloud_scattering_table(liquid_path)
        ice = read_cloud_scattering_table(ice_path)
        @test liquid.medium == "liquid-water"
        @test ice.medium == "ice"
        @test size(liquid.mass_extinction_coefficient) ==
              (length(liquid.wavenumber), length(liquid.effective_radius))
        @test size(ice.mass_extinction_coefficient) ==
              (length(ice.wavenumber), length(ice.effective_radius))
        @test minimum(liquid.mass_extinction_coefficient) >= 0
        @test minimum(ice.mass_extinction_coefficient) >= 0
        @test all(0 .<= liquid.single_scattering_albedo .<= 1)
        @test all(0 .<= ice.single_scattering_albedo .<= 1)
        @test all(-1 .<= liquid.asymmetry_factor .<= 1)
        @test all(-1 .<= ice.asymmetry_factor .<= 1)

        sw_mapping = read_ecckd_spectral_mapping(joinpath(root, "validation",
            "external", "ecrad", "data",
            "ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc"))
        liquid_gpoints = cloud_scattering_gpoint_properties(liquid, sw_mapping, 10.0e-6)
        ice_gpoints = cloud_scattering_gpoint_properties(ice, sw_mapping, 30.0e-6)
        liquid_ecrad_gpoints = cloud_scattering_gpoint_properties(
            liquid, sw_mapping, 10.0e-6;
            mapping_method = :ecrad,
            delta_eddington_average = true,
        )
        @test length(liquid_gpoints.mass_extinction_coefficient) == 32
        @test length(ice_gpoints.mass_extinction_coefficient) == 32
        @test length(liquid_ecrad_gpoints.mass_extinction_coefficient) == 32
        @test minimum(liquid_gpoints.mass_extinction_coefficient) >= 0
        @test minimum(ice_gpoints.mass_extinction_coefficient) >= 0
        @test minimum(liquid_ecrad_gpoints.mass_extinction_coefficient) >= 0
        @test all(0 .<= liquid_gpoints.single_scattering_albedo .<= 1)
        @test all(0 .<= ice_gpoints.single_scattering_albedo .<= 1)
        @test all(0 .<= liquid_ecrad_gpoints.single_scattering_albedo .<= 1)
    else
        @info "Skipping official cloud scattering table check; ecRad data files are not present" liquid_path ice_path
        @test_skip "official cloud scattering table files are not present"
    end
end
