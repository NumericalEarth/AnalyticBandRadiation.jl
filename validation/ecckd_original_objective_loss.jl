using LinearAlgebra

const ECCKD_HR_SECONDS_PER_DAY = 86400.0

function ecckd_interface_weight(layer_weight, flux_profile_weight)
    length(layer_weight) >= 1 || throw(ArgumentError("layer_weight must be nonempty"))
    return flux_profile_weight .* 0.5 .* (layer_weight[1:(end - 1)] .+
                                          layer_weight[2:end])
end

function squared_sum(values)
    total = zero(eltype(values))
    for value in values
        total += value * value
    end
    return total
end

function checked_objective_shapes(heating_rate_fwd, heating_rate_true,
                                  flux_dn_fwd, flux_up_fwd,
                                  flux_dn_true, flux_up_true,
                                  layer_weight)
    size(heating_rate_fwd) == size(heating_rate_true) ||
        throw(DimensionMismatch("heating-rate arrays must have matching shape"))
    size(flux_dn_fwd) == size(flux_dn_true) ||
        throw(DimensionMismatch("downwelling flux arrays must have matching shape"))
    size(flux_up_fwd) == size(flux_up_true) ||
        throw(DimensionMismatch("upwelling flux arrays must have matching shape"))
    size(flux_dn_fwd) == size(flux_up_fwd) ||
        throw(DimensionMismatch("upwelling and downwelling flux arrays must have matching shape"))
    size(flux_dn_fwd, 1) == size(heating_rate_fwd, 1) + 1 ||
        throw(DimensionMismatch("flux arrays must have one more vertical level than heating-rate arrays"))
    size(flux_dn_fwd, 2) == size(heating_rate_fwd, 2) ||
        throw(DimensionMismatch("flux and heating-rate arrays must have matching band count"))
    length(layer_weight) == size(heating_rate_fwd, 1) ||
        throw(DimensionMismatch("layer_weight must match layer count"))
    return size(heating_rate_fwd, 1), size(heating_rate_fwd, 2)
end

function ecckd_lw_ckd_loss(; heating_rate_fwd,
                           heating_rate_true,
                           flux_dn_fwd,
                           flux_up_fwd,
                           flux_dn_true,
                           flux_up_true,
                           layer_weight,
                           flux_weight,
                           flux_profile_weight,
                           broadband_weight,
                           spectral_boundary_weight = 0.0,
                           flux_dn_fwd_orig = flux_dn_fwd,
                           flux_up_fwd_orig = flux_up_fwd,
                           spectral_flux_dn_surf = nothing,
                           spectral_flux_up_toa = nothing)
    nlay, nband = checked_objective_shapes(
        heating_rate_fwd, heating_rate_true, flux_dn_fwd, flux_up_fwd,
        flux_dn_true, flux_up_true, layer_weight,
    )
    cost = zero(promote_type(eltype(heating_rate_fwd), Float64))
    hr_weight2 = ECCKD_HR_SECONDS_PER_DAY^2

    for iband in 1:nband
        cost += hr_weight2 *
            sum(layer_weight .* abs2.(heating_rate_fwd[:, iband] .-
                                      heating_rate_true[:, iband]))
        cost += flux_weight *
            ((flux_dn_fwd[end, iband] - flux_dn_true[end, iband])^2 +
             (flux_up_fwd[1, iband] - flux_up_true[1, iband])^2)
        if flux_profile_weight > 0
            interface_weight = ecckd_interface_weight(layer_weight, flux_profile_weight)
            cost += sum(interface_weight .*
                        (abs2.(flux_dn_fwd[2:nlay, iband] .-
                               flux_dn_true[2:nlay, iband]) .+
                         abs2.(flux_up_fwd[2:nlay, iband] .-
                               flux_up_true[2:nlay, iband])))
        end
    end

    broadband_heating = sum(heating_rate_fwd .- heating_rate_true; dims = 2)
    cost = cost * (1 - broadband_weight) / nband
    cost += broadband_weight * hr_weight2 *
        sum(layer_weight .* vec(abs2.(broadband_heating)))
    cost += broadband_weight * flux_weight *
        ((sum(flux_dn_fwd[end, :] .- flux_dn_true[end, :]))^2 +
         (sum(flux_up_fwd[1, :] .- flux_up_true[1, :]))^2)

    if flux_profile_weight > 0
        interface_weight = ecckd_interface_weight(layer_weight, flux_profile_weight)
        flux_dn_error = vec(sum(flux_dn_fwd[2:nlay, :] .-
                                flux_dn_true[2:nlay, :]; dims = 2))
        flux_up_error = vec(sum(flux_up_fwd[2:nlay, :] .-
                                flux_up_true[2:nlay, :]; dims = 2))
        cost += broadband_weight *
            sum(interface_weight .* (abs2.(flux_dn_error) .+ abs2.(flux_up_error)))
    end

    if spectral_boundary_weight > 0 &&
       spectral_flux_dn_surf !== nothing &&
       spectral_flux_up_toa !== nothing
        cost += spectral_boundary_weight *
            (squared_sum(flux_dn_fwd_orig[end, :] .- spectral_flux_dn_surf) +
             squared_sum(flux_up_fwd_orig[1, :] .- spectral_flux_up_toa))
    end
    return cost
end

function ecckd_sw_ckd_loss(; heating_rate_fwd,
                           heating_rate_true,
                           flux_dn_fwd,
                           flux_up_fwd,
                           flux_dn_true,
                           flux_up_true,
                           layer_weight,
                           flux_weight,
                           flux_profile_weight,
                           broadband_weight,
                           all_albedo_positive = true,
                           spectral_boundary_weight = 0.0,
                           flux_dn_fwd_orig = flux_dn_fwd,
                           spectral_flux_dn_surf = nothing)
    nlay, nband = checked_objective_shapes(
        heating_rate_fwd, heating_rate_true, flux_dn_fwd, flux_up_fwd,
        flux_dn_true, flux_up_true, layer_weight,
    )
    cost = zero(promote_type(eltype(heating_rate_fwd), Float64))
    hr_weight2 = ECCKD_HR_SECONDS_PER_DAY^2

    if flux_profile_weight == 0 && broadband_weight == 0 &&
       spectral_boundary_weight == 0
        heating_error = heating_rate_fwd .- heating_rate_true
        weighted_heating_error =
            sum(layer_weight .* vec(sum(abs2.(heating_error); dims = 2)))
        surface_down_error = flux_dn_fwd[end, :] .- flux_dn_true[end, :]
        toa_up_error = flux_up_fwd[1, :] .- flux_up_true[1, :]
        return hr_weight2 * weighted_heating_error +
            flux_weight * (sum(abs2, surface_down_error) +
                           20.0 * sum(abs2, toa_up_error))
    end

    for iband in 1:nband
        cost += hr_weight2 *
            sum(layer_weight .* abs2.(heating_rate_fwd[:, iband] .-
                                      heating_rate_true[:, iband]))
        cost += flux_weight *
            ((flux_dn_fwd[end, iband] - flux_dn_true[end, iband])^2 +
             20.0 * (flux_up_fwd[1, iband] - flux_up_true[1, iband])^2)
        if flux_profile_weight > 0
            interface_weight = ecckd_interface_weight(layer_weight, flux_profile_weight)
            cost += sum(interface_weight .*
                        (abs2.(flux_dn_fwd[2:nlay, iband] .-
                               flux_dn_true[2:nlay, iband]) .+
                         abs2.(flux_up_fwd[2:nlay, iband] .-
                               flux_up_true[2:nlay, iband])))
        end
    end

    if broadband_weight > 0
        broadband_heating = sum(heating_rate_fwd .- heating_rate_true; dims = 2)
        cost = cost * (1 - broadband_weight) / nband
        cost += broadband_weight * hr_weight2 *
            sum(layer_weight .* vec(abs2.(broadband_heating)))
        cost += broadband_weight * flux_weight *
            (sum(flux_dn_fwd[end, :] .- flux_dn_true[end, :]))^2
        if all_albedo_positive
            cost += broadband_weight * flux_weight *
                (sum(flux_up_fwd[1, :] .- flux_up_true[1, :]))^2
        end
        if flux_profile_weight > 0
            interface_weight = ecckd_interface_weight(layer_weight, flux_profile_weight)
            flux_dn_error = vec(sum(flux_dn_fwd[2:nlay, :] .-
                                    flux_dn_true[2:nlay, :]; dims = 2))
            cost += broadband_weight * sum(interface_weight .* abs2.(flux_dn_error))
            if all_albedo_positive
                flux_up_error = vec(sum(flux_up_fwd[2:nlay, :] .-
                                        flux_up_true[2:nlay, :]; dims = 2))
                cost += broadband_weight * sum(interface_weight .* abs2.(flux_up_error))
            end
        end
    end

    if spectral_boundary_weight > 0 && spectral_flux_dn_surf !== nothing
        cost += spectral_boundary_weight *
            squared_sum(flux_dn_fwd_orig[end, :] .- spectral_flux_dn_surf)
    end
    return cost
end
