@testset "Layer cloud optics" begin
    @testset "cloud optical properties" begin
        cloud = CloudOpticalProperties(zeros(3), zeros(3))
        model = LayerCloudOpticsModel(cloud_water_path = [0.0, 0.05, 0.10],
                                      longwave_mass_absorption = 2.0,
                                      shortwave_mass_extinction = 5.0)

        cloud_optical_properties!(cloud, model, (;))

        @test cloud.longwave_optical_depth == [0.0, 0.10, 0.20]
        @test cloud.shortwave_optical_depth == [0.0, 0.25, 0.50]
    end

    @testset "add cloud optical depths to gas optics" begin
        longwave = LongwaveOpticalProperties([0.1, 0.2], [50.0, 60.0])
        shortwave = ShortwaveOpticalProperties([0.01, 0.02])
        cloud = CloudOpticalProperties([0.3, 0.4], [0.5, 0.6])

        add_cloud_optical_depths!(longwave, shortwave, cloud)

        @test longwave.optical_depth ≈ [0.4, 0.6]
        @test shortwave.optical_depth ≈ [0.51, 0.62]
        @test shortwave.rayleigh_optical_depth ≈ [0.0, 0.0]
        @test shortwave.scattering_asymmetry ≈ [0.0, 0.0]
    end

    @testset "cloud shortwave scattering is composed separately" begin
        longwave = LongwaveOpticalProperties([0.1, 0.2], [50.0, 60.0])
        shortwave = ShortwaveOpticalProperties([0.01, 0.02])
        cloud = CloudOpticalProperties([0.3, 0.4], [0.5, 0.6];
                                       shortwave_scattering_optical_depth = [0.05, 0.06],
                                       shortwave_scattering_asymmetry = [0.7, 0.8])

        add_cloud_optical_depths!(longwave, shortwave, cloud)

        @test longwave.optical_depth ≈ [0.4, 0.6]
        @test shortwave.optical_depth ≈ [0.51, 0.62]
        @test shortwave.rayleigh_optical_depth ≈ [0.05, 0.06]
        @test shortwave.scattering_asymmetry ≈ [0.7, 0.8]
    end

    @testset "cloud scattering asymmetry is optical-depth weighted" begin
        longwave = LongwaveOpticalProperties([0.1], [50.0])
        shortwave = ShortwaveOpticalProperties([0.01];
                                               scattering_optical_depth = [0.2],
                                               scattering_asymmetry = [0.1])
        cloud = CloudOpticalProperties([0.0], [0.0];
                                       shortwave_scattering_optical_depth = [0.3],
                                       shortwave_scattering_asymmetry = [0.6])

        add_cloud_optical_depths!(longwave, shortwave, cloud)

        @test shortwave.rayleigh_optical_depth ≈ [0.5]
        @test shortwave.scattering_asymmetry ≈ [(0.2 * 0.1 + 0.3 * 0.6) / 0.5]
    end

    @testset "liquid/ice cloud optical properties with cloud fraction" begin
        cloud = CloudOpticalProperties(zeros(3), zeros(3))
        model = LayerLiquidIceCloudOpticsModel(
            liquid_water_path = [0.10, 0.20, 0.30],
            ice_water_path = [0.30, 0.20, 0.10],
            cloud_fraction = [1.0, 0.5, 0.0],
            liquid_longwave_mass_absorption = 2.0,
            ice_longwave_mass_absorption = 4.0,
            liquid_shortwave_mass_extinction = 5.0,
            ice_shortwave_mass_extinction = 7.0,
        )

        cloud_optical_properties!(cloud, model, (;))

        @test cloud.longwave_optical_depth ≈ [1.4, 0.6, 0.0]
        @test cloud.shortwave_optical_depth ≈ [2.6, 1.2, 0.0]
        @test cloud.shortwave_scattering_optical_depth ≈ [0.0, 0.0, 0.0]
    end

    @testset "liquid/ice cloud scattering partition" begin
        cloud = CloudOpticalProperties(zeros(2), zeros(2))
        model = LayerLiquidIceCloudOpticsModel(
            liquid_water_path = [0.10, 0.20],
            ice_water_path = [0.30, 0.20],
            cloud_fraction = [1.0, 0.5],
            liquid_longwave_mass_absorption = 2.0,
            ice_longwave_mass_absorption = 4.0,
            liquid_shortwave_mass_extinction = 5.0,
            ice_shortwave_mass_extinction = 7.0,
            liquid_shortwave_single_scattering_albedo = 0.8,
            ice_shortwave_single_scattering_albedo = 0.5,
            liquid_shortwave_scattering_asymmetry = 0.85,
            ice_shortwave_scattering_asymmetry = 0.75,
        )

        cloud_optical_properties!(cloud, model, (;))

        @test cloud.longwave_optical_depth ≈ [1.4, 0.6]
        @test cloud.shortwave_optical_depth ≈ [1.15, 0.45]
        @test cloud.shortwave_scattering_optical_depth ≈ [1.45, 0.75]
        @test cloud.shortwave_scattering_asymmetry ≈ [
            (0.85 * 0.4 + 0.75 * 1.05) / 1.45,
            (0.85 * 0.8 + 0.75 * 0.7) / 1.5,
        ]
    end

    @testset "liquid/ice cloud fraction exponent" begin
        cloud = CloudOpticalProperties(zeros(1), zeros(1))
        model = LayerLiquidIceCloudOpticsModel(
            liquid_water_path = [0.10],
            ice_water_path = [0.30],
            cloud_fraction = [0.25],
            liquid_longwave_mass_absorption = 2.0,
            ice_longwave_mass_absorption = 4.0,
            liquid_shortwave_mass_extinction = 5.0,
            ice_shortwave_mass_extinction = 7.0,
            cloud_fraction_exponent = 0.5,
        )

        cloud_optical_properties!(cloud, model, (;))

        @test cloud.longwave_optical_depth ≈ [0.5 * 1.4]
        @test cloud.shortwave_optical_depth ≈ [0.5 * 2.6]
    end

    @testset "cloudy-region liquid/ice optical properties keep cloud fraction separate" begin
        cloud = CloudyRegionCloudOpticalProperties(
            zeros(2),
            zeros(1),
            zeros(2),
            zeros(2),
            shortwave_scattering_optical_depth = zeros(2),
            shortwave_scattering_asymmetry = zeros(2),
        )
        model = LayerLiquidIceCloudOpticsModel(
            liquid_water_path = [0.10, 0.20],
            ice_water_path = [0.30, 0.20],
            cloud_fraction = [0.25, 0.50],
            liquid_longwave_mass_absorption = 2.0,
            ice_longwave_mass_absorption = 4.0,
            liquid_shortwave_mass_extinction = 5.0,
            ice_shortwave_mass_extinction = 7.0,
            liquid_shortwave_single_scattering_albedo = 0.8,
            ice_shortwave_single_scattering_albedo = 0.5,
            liquid_shortwave_scattering_asymmetry = 0.85,
            ice_shortwave_scattering_asymmetry = 0.75,
            cloud_fraction_exponent = 0.5,
        )
        atmosphere = (; overlap_parameter = [0.7])

        cloudy_region_optical_properties!(cloud, model, atmosphere)

        @test cloud.cloud_fraction ≈ [0.25, 0.50]
        @test cloud.overlap_parameter ≈ [0.7]
        @test cloud.longwave_optical_depth ≈ [1.4, 1.2]
        @test cloud.shortwave_optical_depth ≈ [1.15, 0.9]
        @test cloud.shortwave_scattering_optical_depth ≈ [1.45, 1.5]
        @test cloud.shortwave_scattering_asymmetry ≈ [
            (0.85 * 0.4 + 0.75 * 1.05) / 1.45,
            (0.85 * 0.8 + 0.75 * 0.7) / 1.5,
        ]
    end

    @testset "liquid/ice model reads atmosphere properties" begin
        cloud = CloudOpticalProperties(zeros(2), zeros(2))
        model = LayerLiquidIceCloudOpticsModel(
            liquid_water_path = 0.0,
            ice_water_path = 0.0,
            cloud_fraction = 0.0,
            liquid_longwave_mass_absorption = 1.0,
            ice_longwave_mass_absorption = 2.0,
            liquid_shortwave_mass_extinction = 3.0,
            ice_shortwave_mass_extinction = 4.0,
        )
        atmosphere = (;
            liquid_water_path = [0.1, 0.2],
            ice_water_path = [0.3, 0.4],
            cloud_fraction = [0.5, 1.0],
        )

        cloud_optical_properties!(cloud, model, atmosphere)

        @test cloud.longwave_optical_depth ≈ [0.35, 1.0]
        @test cloud.shortwave_optical_depth ≈ [0.75, 2.2]
    end

    @testset "absorptive all-sky shortwave reduces surface flux" begin
        nlayers = 2
        clear_fluxes = RadiativeFluxes(
            longwave_up = zeros(nlayers + 1),
            longwave_down = zeros(nlayers + 1),
            shortwave_up = zeros(nlayers + 1),
            shortwave_down = zeros(nlayers + 1),
        )
        cloudy_fluxes = RadiativeFluxes(
            longwave_up = zeros(nlayers + 1),
            longwave_down = zeros(nlayers + 1),
            shortwave_up = zeros(nlayers + 1),
            shortwave_down = zeros(nlayers + 1),
        )
        clear_shortwave = ShortwaveOpticalProperties(zeros(nlayers))
        cloudy_shortwave = ShortwaveOpticalProperties(zeros(nlayers))
        longwave = LongwaveOpticalProperties(zeros(nlayers), zeros(nlayers))
        cloud = CloudOpticalProperties(zeros(nlayers), [log(2.0), 0.0])
        boundary = ShortwaveBoundaryConditions(toa_shortwave_down = 400.0,
                                               surface_albedo = 0.0)

        add_cloud_optical_depths!(longwave, cloudy_shortwave, cloud)
        radiative_fluxes!(clear_fluxes, CloudlessShortwave(), clear_shortwave, nothing, boundary)
        radiative_fluxes!(cloudy_fluxes, CloudlessShortwave(), cloudy_shortwave, nothing, boundary)

        @test clear_fluxes.shortwave_down[end] == 400.0
        @test cloudy_fluxes.shortwave_down[end] ≈ 200.0
    end

    @testset "mapped cloud scattering composition" begin
        shortwave = ShortwaveOpticalProperties(zeros(2, 2);
                                               scattering_optical_depth = zeros(2, 2),
                                               scattering_asymmetry = zeros(2, 2))
        liquid = (
            mass_extinction_coefficient = [10.0, 20.0],
            single_scattering_albedo = [0.8, 0.5],
            asymmetry_factor = [0.7, 0.6],
        )
        ice = (
            mass_extinction_coefficient = [30.0, 40.0],
            single_scattering_albedo = [0.6, 0.25],
            asymmetry_factor = [0.5, 0.4],
        )

        add_mapped_cloud_scattering!(shortwave, liquid, ice,
                                     [0.1, 0.0], [0.0, 0.2], [1.0, 0.5])

        @test shortwave.optical_depth[:, 1] ≈ [0.2, 1.0]
        @test shortwave.rayleigh_optical_depth[:, 1] ≈ [0.8, 1.0]
        @test shortwave.scattering_asymmetry[:, 1] ≈ [0.7, 0.6]
        @test shortwave.optical_depth[:, 2] ≈ [1.2, 3.0]
        @test shortwave.rayleigh_optical_depth[:, 2] ≈ [1.8, 1.0]
        @test shortwave.scattering_asymmetry[:, 2] ≈ [0.5, 0.4]

        scaled = ShortwaveOpticalProperties(zeros(2, 1);
                                            scattering_optical_depth = zeros(2, 1),
                                            scattering_asymmetry = zeros(2, 1))
        add_mapped_cloud_scattering!(scaled, liquid, ice,
                                     [0.1], [0.0], [1.0];
                                     liquid_extinction_scale = 2.0)
        @test scaled.optical_depth[:, 1] ≈ [0.4, 2.0]
        @test scaled.rayleigh_optical_depth[:, 1] ≈ [1.6, 2.0]

        delta_scaled = ShortwaveOpticalProperties(zeros(1, 1);
                                                  scattering_optical_depth = zeros(1, 1),
                                                  scattering_asymmetry = zeros(1, 1))
        add_mapped_cloud_scattering!(
            delta_scaled,
            (mass_extinction_coefficient = [10.0],
             single_scattering_albedo = [0.8],
             asymmetry_factor = [0.7]),
            (mass_extinction_coefficient = [0.0],
             single_scattering_albedo = [0.0],
             asymmetry_factor = [0.0]),
            [0.1],
            [0.0],
            [1.0];
            delta_eddington_scale = true,
        )
        @test delta_scaled.optical_depth[1, 1] ≈ 0.2
        @test delta_scaled.rayleigh_optical_depth[1, 1] ≈ 0.8 * (1 - 0.7^2)
        @test delta_scaled.scattering_asymmetry[1, 1] ≈ 0.7 / 1.7
    end
end

@testset "Layer aerosol optics" begin
    @testset "aerosol optical properties" begin
        aerosol = AerosolOpticalProperties(zeros(3), zeros(3))
        model = LayerAerosolOpticsModel(aerosol_path = [0.0, 0.02, 0.04],
                                        longwave_mass_absorption = 1.5,
                                        shortwave_mass_extinction = 4.0)

        aerosol_optical_properties!(aerosol, model, (;))

        @test aerosol.longwave_optical_depth == [0.0, 0.03, 0.06]
        @test aerosol.shortwave_optical_depth == [0.0, 0.08, 0.16]
        @test aerosol.shortwave_scattering_optical_depth == [0.0, 0.0, 0.0]
        @test aerosol.shortwave_scattering_asymmetry == [0.0, 0.0, 0.0]
    end

    @testset "aerosol shortwave scattering partition" begin
        aerosol = AerosolOpticalProperties(zeros(2), zeros(2))
        model = LayerAerosolOpticsModel(aerosol_path = [0.02, 0.04],
                                        longwave_mass_absorption = 1.5,
                                        shortwave_mass_extinction = 4.0,
                                        shortwave_single_scattering_albedo = 0.75,
                                        shortwave_scattering_asymmetry = 0.65)

        aerosol_optical_properties!(aerosol, model, (;))

        @test aerosol.longwave_optical_depth ≈ [0.03, 0.06]
        @test aerosol.shortwave_optical_depth ≈ [0.02, 0.04]
        @test aerosol.shortwave_scattering_optical_depth ≈ [0.06, 0.12]
        @test aerosol.shortwave_scattering_asymmetry ≈ [0.65, 0.65]
    end

    @testset "add aerosol optical depths to gas optics" begin
        longwave = LongwaveOpticalProperties([0.1, 0.2], [50.0, 60.0])
        shortwave = ShortwaveOpticalProperties([0.01, 0.02])
        aerosol = AerosolOpticalProperties([0.03, 0.04], [0.05, 0.06];
                                           shortwave_scattering_optical_depth = [0.01, 0.02],
                                           shortwave_scattering_asymmetry = [0.3, 0.4])

        add_aerosol_optical_depths!(longwave, shortwave, aerosol)

        @test longwave.optical_depth ≈ [0.13, 0.24]
        @test shortwave.optical_depth ≈ [0.06, 0.08]
        @test shortwave.rayleigh_optical_depth ≈ [0.01, 0.02]
        @test shortwave.scattering_asymmetry ≈ [0.3, 0.4]
    end

    @testset "absorptive gas cloud aerosol composition" begin
        longwave = LongwaveOpticalProperties([0.1, 0.2], [50.0, 60.0])
        shortwave = ShortwaveOpticalProperties([0.01, 0.02])
        cloud = CloudOpticalProperties([0.3, 0.4], [0.5, 0.6])
        aerosol = AerosolOpticalProperties([0.03, 0.04], [0.05, 0.06])

        add_cloud_optical_depths!(longwave, shortwave, cloud)
        add_aerosol_optical_depths!(longwave, shortwave, aerosol)

        @test longwave.optical_depth ≈ [0.43, 0.64]
        @test shortwave.optical_depth ≈ [0.56, 0.68]
        @test shortwave.rayleigh_optical_depth ≈ [0.0, 0.0]
    end
end
