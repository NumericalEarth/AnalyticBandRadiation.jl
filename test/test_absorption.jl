@testset "h2o_line_kappa_ref: band structure" begin
    lw = WilliamsLongwave(Float32)
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
    lw = WilliamsLongwave(Float32)
    @test co2_kappa_ref(667f0, lw) ≈ lw.κ_co2
    @test co2_kappa_ref(667f0 + lw.l_co2, lw) ≈ lw.κ_co2 / Float32(ℯ) rtol = 1f-4
    @test co2_kappa_ref(400f0, lw) == 0f0
    @test co2_kappa_ref(900f0, lw) == 0f0
end

@testset "h2o_cont_kappa_ref: two-band split" begin
    lw = WilliamsLongwave(Float32)
    @test h2o_cont_kappa_ref(1000f0, lw) == lw.κ_cnt1
    @test h2o_cont_kappa_ref(2000f0, lw) == lw.κ_cnt2
    @test h2o_cont_kappa_ref(1700f0, lw) == lw.κ_cnt1
    @test lw.κ_cnt1 > lw.κ_cnt2
end

@testset "williams_delta_tau: positivity + rotation-band dominance" begin
    NF = Float32
    nlayers = 4
    σ_half = collect(NF.(range(0, 1, length = nlayers + 1)))
    geom = ColumnGeometry(σ_half)
    T = fill(NF(280), nlayers)
    q = fill(NF(0.005), nlayers)
    pₛ = NF(100_000)
    g  = NF(9.80665)

    lw_h2o = WilliamsLongwave(NF; do_co2 = false)
    lw_co2 = WilliamsLongwave(NF; do_co2 = true)

    dν̃ = (lw_h2o.wavenumber_max - lw_h2o.wavenumber_min) / (lw_h2o.nwavenumber - 1)

    for k in 1:nlayers
        Δτ_win = AnalyticBandRadiation.williams_delta_tau(k, NF(1000), T, q, pₛ, geom, lw_h2o, g)
        Δτ_rot = AnalyticBandRadiation.williams_delta_tau(k, NF(400),  T, q, pₛ, geom, lw_h2o, g)
        Δτ_no_co2   = AnalyticBandRadiation.williams_delta_tau(k, NF(667), T, q, pₛ, geom, lw_h2o, g)
        Δτ_with_co2 = AnalyticBandRadiation.williams_delta_tau(k, NF(667), T, q, pₛ, geom, lw_co2, g)

        @test Δτ_win > 0
        @test Δτ_rot > Δτ_win
        @test Δτ_with_co2 > Δτ_no_co2

        for iv in 1:lw_h2o.nwavenumber
            ν̃ = lw_h2o.wavenumber_min + (iv - 1) * dν̃
            @test AnalyticBandRadiation.williams_delta_tau(k, NF(ν̃), T, q, pₛ, geom, lw_h2o, g) >= 0
        end
    end
end
