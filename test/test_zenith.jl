@testset "solar_declination: cardinal dates" begin
    # Spencer-series values at solstices and equinoxes, in degrees.
    γ_equinox   = 2π * (79 - 1) / 365.25          # roughly 20 March
    γ_summer    = 2π * (172 - 1) / 365.25         # roughly 21 June
    γ_winter    = 2π * (355 - 1) / 365.25         # roughly 21 December
    @test solar_declination(γ_equinox) * 180 / π ≈ 0 atol = 1
    @test solar_declination(γ_summer) * 180 / π ≈ 23.44 atol = 1
    @test solar_declination(γ_winter) * 180 / π ≈ -23.44 atol = 1
end

@testset "cosine_solar_zenith: nightside = 0" begin
    t = DateTime(2026, 6, 21, 12, 0, 0)
    cosθ_noon_equator = cosine_solar_zenith(0.0, 0.0, t)
    cosθ_midnight = cosine_solar_zenith(π, 0.0, t)
    @test cosθ_noon_equator > 0.8
    @test cosθ_midnight ≈ 0 atol = 1e-6
end
