@testset "planck: Stefan-Boltzmann recovery" begin
    # π ∫ B(T, ν̃) dν̃ from 10 to 2510 cm⁻¹ must be within ~1% of σ T⁴ at 300 K.
    T = 300.0f0
    lw = AnalyticBandLongwave(Float32)
    nν = lw.nwavenumber
    dν̃ = (lw.wavenumber_max - lw.wavenumber_min) / (nν - 1)
    pi_B = sum(Float32(π) * planck_wavenumber(T, lw.wavenumber_min + (iv-1)*dν̃) * dν̃
               for iv in 1:nν)
    σ_SB = 5.670374419f-8
    @test pi_B ≈ σ_SB * T^4 rtol = 0.01
end

@testset "planck: monotonic in temperature" begin
    ν̃ = 667.0
    @test planck_wavenumber(250.0, ν̃) < planck_wavenumber(290.0, ν̃)
    @test planck_wavenumber(290.0, ν̃) < planck_wavenumber(320.0, ν̃)
end
