@testset "CloudlessLongwave precomputed-optics solver" begin
    nlayers = 3
    fluxes = RadiativeFluxes(
        longwave_up = zeros(nlayers + 1),
        longwave_down = zeros(nlayers + 1),
        shortwave_up = zeros(nlayers + 1),
        shortwave_down = zeros(nlayers + 1),
    )
    atmosphere = nothing

    @testset "no-atmosphere limit" begin
        optics = LongwaveOpticalProperties(zeros(nlayers), zeros(nlayers))
        boundary = LongwaveBoundaryConditions(surface_longwave_up = 300.0)

        radiative_fluxes!(fluxes, CloudlessLongwave(), optics, atmosphere, boundary)

        @test fluxes.longwave_up == fill(300.0, nlayers + 1)
        @test fluxes.longwave_down == zeros(nlayers + 1)
    end

    @testset "single-layer absorber" begin
        tau = [log(2.0)]
        source = [100.0]
        one_layer_fluxes = RadiativeFluxes(
            longwave_up = zeros(2),
            longwave_down = zeros(2),
            shortwave_up = zeros(2),
            shortwave_down = zeros(2),
        )
        optics = LongwaveOpticalProperties(tau, source)
        boundary = LongwaveBoundaryConditions(surface_longwave_up = 300.0)

        radiative_fluxes!(one_layer_fluxes, CloudlessLongwave(), optics, atmosphere, boundary)

        @test one_layer_fluxes.longwave_up[2] == 300.0
        @test one_layer_fluxes.longwave_up[1] ≈ 200.0
        @test one_layer_fluxes.longwave_down[1] == 0.0
        @test one_layer_fluxes.longwave_down[2] ≈ 50.0
    end

    @testset "two-layer pure absorber closed-form Schwarzschild solution" begin
        tau = [0.3, 1.1]
        source = [180.0, 260.0]
        surface_up = 340.0
        toa_down = 25.0
        two_layer_fluxes = RadiativeFluxes(
            longwave_up = zeros(3),
            longwave_down = zeros(3),
            shortwave_up = zeros(3),
            shortwave_down = zeros(3),
        )
        optics = LongwaveOpticalProperties(tau, source)
        boundary = LongwaveBoundaryConditions(surface_longwave_up = surface_up,
                                              toa_longwave_down = toa_down)

        radiative_fluxes!(two_layer_fluxes, CloudlessLongwave(), optics,
                          atmosphere, boundary)

        transmittance = exp.(-tau)
        expected_up = zeros(3)
        expected_down = zeros(3)
        expected_up[3] = surface_up
        expected_up[2] = surface_up * transmittance[2] +
            source[2] * (1 - transmittance[2])
        expected_up[1] = expected_up[2] * transmittance[1] +
            source[1] * (1 - transmittance[1])
        expected_down[1] = toa_down
        expected_down[2] = toa_down * transmittance[1] +
            source[1] * (1 - transmittance[1])
        expected_down[3] = expected_down[2] * transmittance[2] +
            source[2] * (1 - transmittance[2])

        @test two_layer_fluxes.longwave_up ≈ expected_up rtol = 1.0e-12 atol = 1.0e-12
        @test two_layer_fluxes.longwave_down ≈ expected_down rtol = 1.0e-12 atol = 1.0e-12
    end

    @testset "weighted spectral accumulation" begin
        tau = [0.0 0.0;
               log(2.0) log(2.0)]
        source = [0.0 0.0;
                  100.0 100.0]
        weights = [0.25, 0.75]
        two_layer_fluxes = RadiativeFluxes(
            longwave_up = zeros(3),
            longwave_down = zeros(3),
            shortwave_up = zeros(3),
            shortwave_down = zeros(3),
        )
        optics = LongwaveOpticalProperties(tau, source; weights)
        boundary = LongwaveBoundaryConditions(surface_longwave_up = 300.0)

        radiative_fluxes!(two_layer_fluxes, CloudlessLongwave(), optics, atmosphere, boundary)

        @test two_layer_fluxes.longwave_up[3] == 300.0
        @test two_layer_fluxes.longwave_up[1] ≈ 187.5
        @test two_layer_fluxes.longwave_down[1] == 0.0
        @test two_layer_fluxes.longwave_down[3] ≈ 56.25
    end

    @testset "spectral surface boundary" begin
        tau = [0.0;
               log(2.0)][:, :]
        source = zeros(2, 1)
        weights = [0.25, 0.75]
        one_layer_fluxes = RadiativeFluxes(
            longwave_up = zeros(2),
            longwave_down = zeros(2),
            shortwave_up = zeros(2),
            shortwave_down = zeros(2),
        )
        optics = LongwaveOpticalProperties(tau, source; weights)
        boundary = LongwaveBoundaryConditions(surface_longwave_up = [100.0, 500.0])

        radiative_fluxes!(one_layer_fluxes, CloudlessLongwave(), optics, atmosphere, boundary)

        @test one_layer_fluxes.longwave_up[2] ≈ 400.0
        @test one_layer_fluxes.longwave_up[1] ≈ 25.0 + 187.5
        @test one_layer_fluxes.longwave_down == zeros(2)
    end

    @testset "ecRad no-scattering interface sources" begin
        tau = [0.5]
        source = [125.0]
        source_top = [100.0]
        source_bottom = [200.0]
        one_layer_fluxes = RadiativeFluxes(
            longwave_up = zeros(2),
            longwave_down = zeros(2),
            shortwave_up = zeros(2),
            shortwave_down = zeros(2),
        )
        optics = LongwaveOpticalProperties(tau, source; source_top, source_bottom)
        boundary = LongwaveBoundaryConditions(surface_longwave_up = 300.0)

        radiative_fluxes!(one_layer_fluxes, CloudlessLongwave(), optics, atmosphere, boundary)

        coeff = 1.66 * tau[1]
        transmittance = exp(-coeff)
        gradient = (source_bottom[1] - source_top[1]) / coeff
        expected_source_up = gradient + source_top[1] -
            transmittance * (gradient + source_bottom[1])
        expected_source_down = -gradient + source_bottom[1] -
            transmittance * (-gradient + source_top[1])

        @test one_layer_fluxes.longwave_up[2] == 300.0
        @test one_layer_fluxes.longwave_up[1] ≈ 300.0 * transmittance + expected_source_up
        @test one_layer_fluxes.longwave_down[1] == 0.0
        @test one_layer_fluxes.longwave_down[2] ≈ expected_source_down
    end

    @testset "longwave scattering reduces to no scattering for zero ssa" begin
        tau = [0.5 0.2]
        source = [125.0 175.0]
        source_top = [100.0 150.0]
        source_bottom = [200.0 250.0]
        no_scattering_fluxes = RadiativeFluxes(
            longwave_up = zeros(3),
            longwave_down = zeros(3),
            shortwave_up = zeros(3),
            shortwave_down = zeros(3),
        )
        scattering_fluxes = RadiativeFluxes(
            longwave_up = zeros(3),
            longwave_down = zeros(3),
            shortwave_up = zeros(3),
            shortwave_down = zeros(3),
        )
        no_scattering = LongwaveOpticalProperties(
            tau, source; source_top, source_bottom)
        scattering = LongwaveOpticalProperties(
            tau, source;
            source_top,
            source_bottom,
            single_scattering_albedo = zeros(1, 2),
            scattering_asymmetry = zeros(1, 2),
        )
        boundary = LongwaveBoundaryConditions(surface_longwave_up = 300.0)

        radiative_fluxes!(no_scattering_fluxes, CloudlessLongwave(),
                          no_scattering, atmosphere, boundary)
        radiative_fluxes!(scattering_fluxes, CloudlessLongwave(),
                          scattering, atmosphere, boundary)

        @test scattering_fluxes.longwave_up ≈ no_scattering_fluxes.longwave_up
        @test scattering_fluxes.longwave_down ≈ no_scattering_fluxes.longwave_down
    end
end

@testset "CloudOverlapLongwave all-sky access point" begin
    atmosphere = nothing
    clear = LongwaveOpticalProperties(
        [0.1 0.1],
        zeros(1, 2);
        source_top = [100.0 120.0],
        source_bottom = [120.0 140.0],
    )
    cloudy = LongwaveOpticalProperties(
        [0.5 0.5],
        zeros(1, 2);
        source_top = [100.0 120.0],
        source_bottom = [120.0 140.0],
    )
    clear_fluxes = RadiativeFluxes(
        longwave_up = zeros(3),
        longwave_down = zeros(3),
        shortwave_up = zeros(3),
        shortwave_down = zeros(3),
    )
    cloudy_fluxes = RadiativeFluxes(
        longwave_up = zeros(3),
        longwave_down = zeros(3),
        shortwave_up = zeros(3),
        shortwave_down = zeros(3),
    )
    overlap_fluxes = RadiativeFluxes(
        longwave_up = zeros(3),
        longwave_down = zeros(3),
        shortwave_up = zeros(3),
        shortwave_down = zeros(3),
    )
    boundary = LongwaveBoundaryConditions(surface_longwave_up = 300.0)

    radiative_fluxes!(clear_fluxes, CloudlessLongwave(), clear, atmosphere, boundary)
    radiative_fluxes!(cloudy_fluxes, CloudlessLongwave(), cloudy, atmosphere, boundary)

    clear_overlap = LongwaveCloudOverlapOpticalProperties(clear, cloudy, [0.0, 0.0])
    radiative_fluxes!(overlap_fluxes, CloudOverlapLongwave(), clear_overlap,
                      atmosphere, boundary)
    @test overlap_fluxes.longwave_up ≈ clear_fluxes.longwave_up
    @test overlap_fluxes.longwave_down ≈ clear_fluxes.longwave_down

    full_cloud_overlap =
        LongwaveCloudOverlapOpticalProperties(clear, cloudy, [1.0, 1.0])
    radiative_fluxes!(overlap_fluxes, CloudOverlapLongwave(), full_cloud_overlap,
                      atmosphere, boundary)
    @test overlap_fluxes.longwave_up ≈ cloudy_fluxes.longwave_up
    @test overlap_fluxes.longwave_down ≈ cloudy_fluxes.longwave_down

    mixed_overlap = LongwaveCloudOverlapOpticalProperties(clear, cloudy, [0.5, 0.5])
    radiative_fluxes!(overlap_fluxes, CloudOverlapLongwave(), mixed_overlap,
                      atmosphere, boundary)
    @test all(isfinite, overlap_fluxes.longwave_up)
    @test all(isfinite, overlap_fluxes.longwave_down)

    tripleclouds_overlap = LongwaveCloudOverlapOpticalProperties(
        clear, cloudy, [0.2, 0.8];
        overlap_parameter = [0.5],
        fractional_std = [1.5, 2.0],
    )
    radiative_fluxes!(
        overlap_fluxes,
        CloudOverlapLongwave(overlap = :tripleclouds_alpha),
        tripleclouds_overlap,
        atmosphere,
        boundary,
    )
    @test all(isfinite, overlap_fluxes.longwave_up)
    @test all(isfinite, overlap_fluxes.longwave_down)

    clear_tripleclouds =
        LongwaveCloudOverlapOpticalProperties(clear, cloudy, [0.0, 0.0])
    radiative_fluxes!(
        overlap_fluxes,
        CloudOverlapLongwave(overlap = :tripleclouds_alpha),
        clear_tripleclouds,
        atmosphere,
        boundary,
    )
    @test overlap_fluxes.longwave_up ≈ clear_fluxes.longwave_up
    @test overlap_fluxes.longwave_down ≈ clear_fluxes.longwave_down

    @test_throws DimensionMismatch LongwaveCloudOverlapOpticalProperties(
        clear, cloudy, [0.5])
    @test_throws ArgumentError CloudOverlapLongwave(overlap = :invalid)
end
