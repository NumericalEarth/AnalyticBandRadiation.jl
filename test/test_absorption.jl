@testset "h2o_line_kappa_ref: band structure" begin
    lw = AnalyticBandLongwave(Float32)
    @test h2o_line_kappa_ref(100f0, lw) ≈ lw.κ_rot
    @test h2o_line_kappa_ref(200f0, lw) ≈ lw.κ_rot
    @test h2o_line_kappa_ref(600f0, lw) < lw.κ_rot
    @test h2o_line_kappa_ref(600f0, lw) ≈ lw.κ_rot * exp(-(600f0 - 200f0) / lw.l_rot) rtol = 1f-5
    @test h2o_line_kappa_ref(1450f0, lw) ≈ lw.κ_vr
    @test h2o_line_kappa_ref(1600f0, lw) ≈ lw.κ_vr
    @test h2o_line_kappa_ref(2000f0, lw) < lw.κ_vr
    @test h2o_line_kappa_ref(3000f0, lw) == 0f0
end

@testset "co2_kappa_ref: peak + e-folding" begin
    lw = AnalyticBandLongwave(Float32)
    @test co2_kappa_ref(667f0, lw) ≈ lw.κ_CO₂
    @test co2_kappa_ref(667f0 + lw.l_CO₂, lw) ≈ lw.κ_CO₂ / Float32(ℯ) rtol = 1f-4
    @test co2_kappa_ref(400f0, lw) == 0f0
    @test co2_kappa_ref(900f0, lw) == 0f0
end

@testset "h2o_cont_kappa_ref: two-band split" begin
    lw = AnalyticBandLongwave(Float32)
    @test h2o_cont_kappa_ref(1000f0, lw) == lw.κ_cnt1
    @test h2o_cont_kappa_ref(2000f0, lw) == lw.κ_cnt2
    # Paper convention: 1700 cm⁻¹ belongs to the upper (weaker) band.
    @test h2o_cont_kappa_ref(1700f0, lw) == lw.κ_cnt2
    @test h2o_cont_kappa_ref(1699f0, lw) == lw.κ_cnt1
    @test lw.κ_cnt1 > lw.κ_cnt2
end

@testset "williams_delta_tau: positivity + rotation-band dominance" begin
    NF = Float32
    nlayers = 4
    σ_half = collect(NF.(range(0, 1, length = nlayers + 1)))
    geom = ColumnGrid(σ_half)
    T = fill(NF(280), nlayers)
    q = fill(NF(0.005), nlayers)
    pₛ = NF(100_000)
    g  = NF(9.80665)

    lw = AnalyticBandLongwave(NF)

    dν̃ = (lw.wavenumber_max - lw.wavenumber_min) / (lw.nwavenumber - 1)

    for k in 1:nlayers
        Δτ_win      = NumericalRadiation.williams_delta_tau(k, NF(1000), NF(0),   T, q, pₛ, geom, lw, g)
        Δτ_rot      = NumericalRadiation.williams_delta_tau(k, NF(400),  NF(0),   T, q, pₛ, geom, lw, g)
        Δτ_no_co2   = NumericalRadiation.williams_delta_tau(k, NF(667),  NF(0),   T, q, pₛ, geom, lw, g)
        Δτ_with_co2 = NumericalRadiation.williams_delta_tau(k, NF(667),  NF(280), T, q, pₛ, geom, lw, g)

        @test Δτ_win > 0
        @test Δτ_rot > Δτ_win
        @test Δτ_with_co2 > Δτ_no_co2

        for iv in 1:lw.nwavenumber
            ν̃ = lw.wavenumber_min + (iv - 1) * dν̃
            @test NumericalRadiation.williams_delta_tau(k, NF(ν̃), NF(0), T, q, pₛ, geom, lw, g) >= 0
        end
    end
end
