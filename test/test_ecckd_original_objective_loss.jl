include(joinpath(@__DIR__, "..", "validation", "ecckd_original_objective_loss.jl"))

function finite_difference_gradient(f, x; step = 1.0e-6)
    gradient = similar(x)
    plus = copy(x)
    minus = copy(x)
    for i in eachindex(x)
        plus .= x
        minus .= x
        plus[i] += step
        minus[i] -= step
        gradient[i] = (f(plus) - f(minus)) / (2step)
    end
    return gradient
end

function relative_gradient_error(a, b)
    norm(a .- b) / max(norm(b), eps(eltype(b)))
end

function enzyme_gradient_for_loss(f, x)
    enzyme = Base.require(Base.PkgId(Base.UUID("7da242da-08ed-463a-9acd-ee780be4f1d9"), "Enzyme"))
    gradient = zeros(length(x))
    duplicated = Base.invokelatest(enzyme.Duplicated, copy(x), gradient)
    const_f = Base.invokelatest(enzyme.Const, f)
    Base.invokelatest(enzyme.autodiff, enzyme.Reverse, const_f,
                      enzyme.Active, duplicated)
    return gradient
end

function reactant_compile_loss_probe(x, flux_up_true, flux_dn_true, layer_weight)
    reactant = Base.require(Base.PkgId(Base.UUID("3c362404-f566-11ee-1572-e11a4b42c853"), "Reactant"))
    Base.invokelatest(reactant.set_default_backend, "cpu")
    Core.eval(@__MODULE__, :(const Reactant = $reactant))
    Core.eval(@__MODULE__, quote
        function sw_objective_loss_probe(x, flux_up_true, flux_dn_true, layer_weight)
            flux_up_fwd = reshape(x, 4, 2)
            heating = zeros(eltype(x), 3, 2)
            ecckd_sw_ckd_loss(;
                heating_rate_fwd = heating,
                heating_rate_true = heating,
                flux_dn_fwd = flux_dn_true,
                flux_up_fwd = flux_up_fwd,
                flux_dn_true = flux_dn_true,
                flux_up_true = flux_up_true,
                layer_weight = layer_weight,
                flux_weight = 0.4,
                flux_profile_weight = 0.0,
                broadband_weight = 0.0,
            )
        end
        x_ra = Reactant.to_rarray($x)
        up_ra = Reactant.to_rarray($flux_up_true)
        dn_ra = Reactant.to_rarray($flux_dn_true)
        layer_ra = Reactant.to_rarray($layer_weight)
        Reactant.@compile raise = true raise_first = true sync = true sw_objective_loss_probe(x_ra, up_ra, dn_ra, layer_ra)
    end)
    return true
end

@testset "ecCKD original objective loss assembly" begin
    layer_weight = [0.2, 0.3, 0.5]
    heating_true = zeros(3, 2)
    heating_fwd = [1.0e-6 2.0e-6;
                   -1.0e-6 0.5e-6;
                   0.0 -0.5e-6]
    flux_dn_true = [10.0 20.0;
                    9.0 18.0;
                    8.0 16.0;
                    7.0 14.0]
    flux_up_true = [1.0 2.0;
                    1.2 2.2;
                    1.4 2.4;
                    1.6 2.6]
    flux_dn_fwd = flux_dn_true .+ [0.1 -0.2;
                                   0.05 -0.1;
                                   0.02 -0.03;
                                   0.3 -0.4]
    flux_up_fwd = flux_up_true .+ [0.2 -0.1;
                                   0.04 -0.02;
                                   -0.03 0.01;
                                   0.0 0.0]

    no_profile_lw = ecckd_lw_ckd_loss(;
        heating_rate_fwd = heating_fwd,
        heating_rate_true = heating_true,
        flux_dn_fwd,
        flux_up_fwd,
        flux_dn_true,
        flux_up_true,
        layer_weight,
        flux_weight = 0.2,
        flux_profile_weight = 0.0,
        broadband_weight = 0.0,
    )
    manual_lw = 0.0
    for iband in 1:2
        manual_lw += ECCKD_HR_SECONDS_PER_DAY^2 *
            sum(layer_weight .* abs2.(heating_fwd[:, iband] .- heating_true[:, iband]))
        manual_lw += 0.2 *
            ((flux_dn_fwd[end, iband] - flux_dn_true[end, iband])^2 +
             (flux_up_fwd[1, iband] - flux_up_true[1, iband])^2)
    end
    manual_lw /= 2
    @test no_profile_lw ≈ manual_lw

    spectral_boundary_lw = ecckd_lw_ckd_loss(;
        heating_rate_fwd = heating_true,
        heating_rate_true = heating_true,
        flux_dn_fwd = flux_dn_true,
        flux_up_fwd = flux_up_true,
        flux_dn_true,
        flux_up_true,
        layer_weight,
        flux_weight = 0.0,
        flux_profile_weight = 0.0,
        broadband_weight = 0.0,
        spectral_boundary_weight = 0.1,
        flux_dn_fwd_orig = flux_dn_fwd,
        flux_up_fwd_orig = flux_up_fwd,
        spectral_flux_dn_surf = flux_dn_true[end, :],
        spectral_flux_up_toa = flux_up_true[1, :],
    )
    @test spectral_boundary_lw ≈
          0.1 * (sum(abs2, flux_dn_fwd[end, :] .- flux_dn_true[end, :]) +
                 sum(abs2, flux_up_fwd[1, :] .- flux_up_true[1, :]))

    sw_toa = ecckd_sw_ckd_loss(;
        heating_rate_fwd = heating_true,
        heating_rate_true = heating_true,
        flux_dn_fwd = flux_dn_true,
        flux_up_fwd,
        flux_dn_true,
        flux_up_true,
        layer_weight,
        flux_weight = 0.4,
        flux_profile_weight = 0.0,
        broadband_weight = 0.0,
    )
    @test sw_toa ≈ 0.4 * 20.0 *
                    sum(abs2, flux_up_fwd[1, :] .- flux_up_true[1, :])

    sw_broadband_no_up = ecckd_sw_ckd_loss(;
        heating_rate_fwd = heating_true,
        heating_rate_true = heating_true,
        flux_dn_fwd,
        flux_up_fwd,
        flux_dn_true,
        flux_up_true,
        layer_weight,
        flux_weight = 0.4,
        flux_profile_weight = 0.0,
        broadband_weight = 0.5,
        all_albedo_positive = false,
    )
    sw_broadband_up = ecckd_sw_ckd_loss(;
        heating_rate_fwd = heating_true,
        heating_rate_true = heating_true,
        flux_dn_fwd,
        flux_up_fwd,
        flux_dn_true,
        flux_up_true,
        layer_weight,
        flux_weight = 0.4,
        flux_profile_weight = 0.0,
        broadband_weight = 0.5,
        all_albedo_positive = true,
    )
    @test sw_broadband_up > sw_broadband_no_up

    sw_boundary = ecckd_sw_ckd_loss(;
        heating_rate_fwd = heating_true,
        heating_rate_true = heating_true,
        flux_dn_fwd = flux_dn_true,
        flux_up_fwd = flux_up_true,
        flux_dn_true,
        flux_up_true,
        layer_weight,
        flux_weight = 0.0,
        flux_profile_weight = 0.0,
        broadband_weight = 0.0,
        spectral_boundary_weight = 0.3,
        flux_dn_fwd_orig = flux_dn_fwd,
        spectral_flux_dn_surf = flux_dn_true[end, :],
    )
    @test sw_boundary ≈
          0.3 * sum(abs2, flux_dn_fwd[end, :] .- flux_dn_true[end, :])

    f(x) = ecckd_sw_ckd_loss(;
        heating_rate_fwd = heating_true,
        heating_rate_true = heating_true,
        flux_dn_fwd = flux_dn_true,
        flux_up_fwd = reshape(x, 4, 2),
        flux_dn_true,
        flux_up_true,
        layer_weight,
        flux_weight = 0.4,
        flux_profile_weight = 0.0,
        broadband_weight = 0.0,
    )
    parameters = vec(copy(flux_up_fwd))
    enzyme_gradient = enzyme_gradient_for_loss(f, parameters)
    finite_difference = finite_difference_gradient(f, parameters)
    @test relative_gradient_error(enzyme_gradient, finite_difference) < 1.0e-6
    @test reactant_compile_loss_probe(parameters, flux_up_true, flux_dn_true, layer_weight)
end
