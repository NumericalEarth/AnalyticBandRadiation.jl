@testset "CloudlessShortwave precomputed-optics solver" begin
    @testset "transparent atmosphere with black surface" begin
        nlayers = 3
        fluxes = RadiativeFluxes(
            longwave_up = zeros(nlayers + 1),
            longwave_down = zeros(nlayers + 1),
            shortwave_up = zeros(nlayers + 1),
            shortwave_down = zeros(nlayers + 1),
        )
        optics = ShortwaveOpticalProperties(zeros(nlayers))
        boundary = ShortwaveBoundaryConditions(toa_shortwave_down = 500.0,
                                               surface_albedo = 0.0)

        radiative_fluxes!(fluxes, CloudlessShortwave(), optics, nothing, boundary)

        @test fluxes.shortwave_down == fill(500.0, nlayers + 1)
        @test fluxes.shortwave_up == zeros(nlayers + 1)
    end

    @testset "no-atmosphere horizontal-flux limit with solar geometry" begin
        fluxes = RadiativeFluxes(
            longwave_up = zeros(1),
            longwave_down = zeros(1),
            shortwave_up = zeros(1),
            shortwave_down = zeros(1),
        )
        optics = ShortwaveOpticalProperties(Float64[])
        atmosphere = (; geometry = (; cos_zenith = 0.25))
        boundary = ShortwaveBoundaryConditions(toa_shortwave_down = 700.0,
                                               surface_albedo = 0.9,
                                               surface_albedo_direct = 0.3)

        radiative_fluxes!(fluxes, CloudlessShortwave(), optics, atmosphere, boundary)

        @test fluxes.shortwave_down == [700.0]
        @test fluxes.shortwave_up == [210.0]
        @test fluxes.shortwave_down[1] - fluxes.shortwave_up[1] == 490.0
    end

    @testset "single-layer absorber and reflecting surface" begin
        fluxes = RadiativeFluxes(
            longwave_up = zeros(2),
            longwave_down = zeros(2),
            shortwave_up = zeros(2),
            shortwave_down = zeros(2),
        )
        optics = ShortwaveOpticalProperties([log(2.0)])
        boundary = ShortwaveBoundaryConditions(toa_shortwave_down = 400.0,
                                               surface_albedo = 0.25)

        radiative_fluxes!(fluxes, CloudlessShortwave(), optics, nothing, boundary)

        @test fluxes.shortwave_down[1] == 400.0
        @test fluxes.shortwave_down[2] ≈ 200.0
        @test fluxes.shortwave_up[2] ≈ 50.0
        @test fluxes.shortwave_up[1] ≈ 25.0
    end

    @testset "weighted spectral accumulation" begin
        nlayers = 2
        fluxes = RadiativeFluxes(
            longwave_up = zeros(nlayers + 1),
            longwave_down = zeros(nlayers + 1),
            shortwave_up = zeros(nlayers + 1),
            shortwave_down = zeros(nlayers + 1),
        )
        tau = [0.0 0.0;
               log(2.0) log(2.0)]
        optics = ShortwaveOpticalProperties(tau; weights = [0.25, 0.75])
        boundary = ShortwaveBoundaryConditions(toa_shortwave_down = 400.0,
                                               surface_albedo = 0.25)

        radiative_fluxes!(fluxes, CloudlessShortwave(), optics, nothing, boundary)

        @test fluxes.shortwave_down[1] == 400.0
        @test fluxes.shortwave_down[3] ≈ 175.0
        @test fluxes.shortwave_up[3] ≈ 43.75
        @test fluxes.shortwave_up[1] ≈ 29.6875
    end

    @testset "per-g-point surface albedo" begin
        fluxes = RadiativeFluxes(
            longwave_up = zeros(2),
            longwave_down = zeros(2),
            shortwave_up = zeros(2),
            shortwave_down = zeros(2),
        )
        tau = reshape([0.0, 0.0], 2, 1)
        optics = ShortwaveOpticalProperties(tau; weights = [0.25, 0.75])
        boundary = ShortwaveBoundaryConditions(toa_shortwave_down = 400.0,
                                               surface_albedo = [0.2, 0.4])

        radiative_fluxes!(fluxes, CloudlessShortwave(), optics, nothing, boundary)

        @test fluxes.shortwave_down == fill(400.0, 2)
        @test fluxes.shortwave_up == fill(140.0, 2)
    end

    @testset "direct surface albedo can differ from diffuse albedo" begin
        fluxes = RadiativeFluxes(
            longwave_up = zeros(2),
            longwave_down = zeros(2),
            shortwave_up = zeros(2),
            shortwave_down = zeros(2),
        )
        optics = ShortwaveOpticalProperties([0.0])
        boundary = ShortwaveBoundaryConditions(toa_shortwave_down = 400.0,
                                               surface_albedo = 0.25,
                                               surface_albedo_direct = 0.1)

        radiative_fluxes!(fluxes, CloudlessShortwave(), optics, nothing, boundary)

        @test fluxes.shortwave_down == fill(400.0, 2)
        @test fluxes.shortwave_up == fill(40.0, 2)
    end

    @testset "cos-zenith scales the direct optical path" begin
        fluxes = RadiativeFluxes(
            longwave_up = zeros(2),
            longwave_down = zeros(2),
            shortwave_up = zeros(2),
            shortwave_down = zeros(2),
        )
        optics = ShortwaveOpticalProperties([log(2.0)])
        atmosphere = (; geometry = (; cos_zenith = 0.5))
        boundary = ShortwaveBoundaryConditions(toa_shortwave_down = 400.0,
                                               surface_albedo = 0.25)

        radiative_fluxes!(fluxes, CloudlessShortwave(), optics, atmosphere, boundary)

        @test fluxes.shortwave_down[1] == 400.0
        @test fluxes.shortwave_down[2] ≈ 100.0
        @test fluxes.shortwave_up[2] ≈ 25.0
        @test fluxes.shortwave_up[1] ≈ 6.25
    end

    @testset "Rayleigh scattering uses ecRad-style two-stream adding" begin
        fluxes = RadiativeFluxes(
            longwave_up = zeros(2),
            longwave_down = zeros(2),
            shortwave_up = zeros(2),
            shortwave_down = zeros(2),
        )
        optics = ShortwaveOpticalProperties([0.0];
                                            rayleigh_optical_depth = [log(2.0)])
        boundary = ShortwaveBoundaryConditions(toa_shortwave_down = 400.0,
                                               surface_albedo = 0.0)

        radiative_fluxes!(fluxes, CloudlessShortwave(), optics, nothing, boundary)

        @test fluxes.shortwave_down[1] == 400.0
        @test fluxes.shortwave_down[2] ≈ 296.07982702762007
        @test fluxes.shortwave_up[1] ≈ 103.92017297152903
        @test fluxes.shortwave_up[2] == 0.0
        @test fluxes.shortwave_down[2] + fluxes.shortwave_up[1] ≈ 400.0
    end

    @testset "single-layer conservative scattering closes energy at mu0 one" begin
        fluxes = RadiativeFluxes(
            longwave_up = zeros(2),
            longwave_down = zeros(2),
            shortwave_up = zeros(2),
            shortwave_down = zeros(2),
        )
        optics = ShortwaveOpticalProperties([0.0];
                                            scattering_optical_depth = [0.7],
                                            scattering_asymmetry = [0.0])
        atmosphere = (; geometry = (; cos_zenith = 1.0))
        boundary = ShortwaveBoundaryConditions(toa_shortwave_down = 400.0,
                                               surface_albedo = 0.0)

        radiative_fluxes!(fluxes, CloudlessShortwave(), optics, atmosphere, boundary)

        @test fluxes.shortwave_down[1] == 400.0
        @test fluxes.shortwave_up[2] == 0.0
        @test fluxes.shortwave_down[2] + fluxes.shortwave_up[1] ≈
              fluxes.shortwave_down[1] rtol = 1.0e-10
        @test all(>=(0), fluxes.shortwave_down)
        @test all(>=(0), fluxes.shortwave_up)
    end

    @testset "Rayleigh backscatter compatibility option preserves ecRad path" begin
        fluxes = RadiativeFluxes(
            longwave_up = zeros(2),
            longwave_down = zeros(2),
            shortwave_up = zeros(2),
            shortwave_down = zeros(2),
        )
        optics = ShortwaveOpticalProperties([0.0];
                                            rayleigh_optical_depth = [log(2.0)])
        boundary = ShortwaveBoundaryConditions(toa_shortwave_down = 400.0,
                                               surface_albedo = 0.0)

        radiative_fluxes!(fluxes,
                          CloudlessShortwave(rayleigh_backscatter_fraction = 0.5),
                          optics,
                          nothing,
                          boundary)

        @test fluxes.shortwave_down[2] ≈ 296.07982702762007
        @test fluxes.shortwave_up[1] ≈ 103.92017297152903
    end

    @testset "scattering asymmetry reaches two-stream coefficients" begin
        fluxes = RadiativeFluxes(
            longwave_up = zeros(2),
            longwave_down = zeros(2),
            shortwave_up = zeros(2),
            shortwave_down = zeros(2),
        )
        optics = ShortwaveOpticalProperties([0.0];
                                            scattering_optical_depth = [log(2.0)],
                                            scattering_asymmetry = [0.8])
        boundary = ShortwaveBoundaryConditions(toa_shortwave_down = 400.0,
                                               surface_albedo = 0.0)

        radiative_fluxes!(fluxes, CloudlessShortwave(), optics, nothing, boundary)

        @test fluxes.shortwave_down[1] == 400.0
        @test fluxes.shortwave_down[2] ≈ 407.619005357698
        @test fluxes.shortwave_up[1] == 0.0
        @test fluxes.shortwave_up[2] == 0.0
    end
end

@testset "CloudOverlapShortwave all-sky access point" begin
    nlayers = 2
    fluxes = RadiativeFluxes(
        longwave_up = zeros(nlayers + 1),
        longwave_down = zeros(nlayers + 1),
        shortwave_up = zeros(nlayers + 1),
        shortwave_down = zeros(nlayers + 1),
    )
    clear = ShortwaveOpticalProperties(zeros(nlayers))
    cloudy = ShortwaveOpticalProperties([log(2.0), log(2.0)])
    optics = ShortwaveCloudOverlapOpticalProperties(clear, cloudy, [0.0, 1.0])
    boundary = ShortwaveBoundaryConditions(toa_shortwave_down = 400.0,
                                           surface_albedo = 0.0)

    radiative_fluxes!(fluxes, CloudOverlapShortwave(overlap = :average),
                      optics, nothing, boundary)

    @test fluxes.shortwave_down[1] ≈ 400.0
    @test fluxes.shortwave_down[2] ≈ 0.5 * 400.0 + 0.5 * 200.0
    @test fluxes.shortwave_down[3] ≈ 100.0
    @test fluxes.shortwave_up == zeros(nlayers + 1)

    adding_fluxes = RadiativeFluxes(
        longwave_up = zeros(nlayers + 1),
        longwave_down = zeros(nlayers + 1),
        shortwave_up = zeros(nlayers + 1),
        shortwave_down = zeros(nlayers + 1),
    )
    radiative_fluxes!(adding_fluxes, CloudOverlapShortwave(overlap = :adding),
                      optics, nothing, boundary)
    @test adding_fluxes.shortwave_down[1] ≈ 400.0
    @test adding_fluxes.shortwave_down[2] ≈ 400.0
    @test adding_fluxes.shortwave_down[3] ≈ 200.0
    @test adding_fluxes.shortwave_up == zeros(nlayers + 1)

    matrix_fluxes = RadiativeFluxes(
        longwave_up = zeros(nlayers + 1),
        longwave_down = zeros(nlayers + 1),
        shortwave_up = zeros(nlayers + 1),
        shortwave_down = zeros(nlayers + 1),
    )
    radiative_fluxes!(matrix_fluxes, CloudOverlapShortwave(overlap = :matrix_maximum),
                      optics, nothing, boundary)
    @test matrix_fluxes.shortwave_down[1] ≈ 400.0
    @test 0.0 <= matrix_fluxes.shortwave_down[3] <= 400.0
    @test matrix_fluxes.shortwave_up == zeros(nlayers + 1)

    alpha_fluxes = RadiativeFluxes(
        longwave_up = zeros(nlayers + 1),
        longwave_down = zeros(nlayers + 1),
        shortwave_up = zeros(nlayers + 1),
        shortwave_down = zeros(nlayers + 1),
    )
    alpha_optics = ShortwaveCloudOverlapOpticalProperties(
        clear, cloudy, [0.0, 1.0]; overlap_parameter = [0.5])
    radiative_fluxes!(alpha_fluxes, CloudOverlapShortwave(overlap = :matrix_alpha),
                      alpha_optics, nothing, boundary)
    @test alpha_fluxes.shortwave_down[1] ≈ 400.0
    @test 0.0 <= alpha_fluxes.shortwave_down[3] <= 400.0
    @test alpha_fluxes.shortwave_up == zeros(nlayers + 1)

    tripleclouds_fluxes = RadiativeFluxes(
        longwave_up = zeros(nlayers + 1),
        longwave_down = zeros(nlayers + 1),
        shortwave_up = zeros(nlayers + 1),
        shortwave_down = zeros(nlayers + 1),
    )
    tripleclouds_optics = ShortwaveCloudOverlapOpticalProperties(
        clear,
        cloudy,
        [0.0, 1.0];
        overlap_parameter = [0.5],
        fractional_std = [1.0, 1.0],
    )
    radiative_fluxes!(tripleclouds_fluxes,
                      CloudOverlapShortwave(overlap = :tripleclouds_alpha),
                      tripleclouds_optics, nothing, boundary)
    @test tripleclouds_fluxes.shortwave_down[1] ≈ 400.0
    @test 0.0 <= tripleclouds_fluxes.shortwave_down[3] <= 400.0
    @test tripleclouds_fluxes.shortwave_up == zeros(nlayers + 1)

    @test_throws DimensionMismatch ShortwaveCloudOverlapOpticalProperties(
        clear, cloudy, [0.0, 1.0]; overlap_parameter = [0.5, 0.5])
    @test_throws DimensionMismatch ShortwaveCloudOverlapOpticalProperties(
        clear, cloudy, [0.0, 1.0]; fractional_std = [1.0])
    @test_throws ArgumentError CloudOverlapShortwave(overlap = :invalid)
end
