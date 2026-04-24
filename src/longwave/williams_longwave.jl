"""
$(TYPEDEF)

Williams (2026) "Simple Spectral Model" (SSM) for clear-sky longwave
radiative transfer.

The scheme solves Schwarzschild's two-stream equations

    dF↑/dτ = F↑ − πB(T)
    dF↓/dτ = πB(T) − F↓

at each of `nwavenumber` evenly spaced wavenumbers between `wavenumber_min`
and `wavenumber_max`, with analytic mass absorption coefficients for H₂O
line (rotation + vibration–rotation + combination bands, [`h2o_line_kappa_ref`](@ref)),
a two-band H₂O continuum ([`h2o_cont_kappa_ref`](@ref)) and a Lorentzian CO₂
15 μm bending mode ([`co2_kappa_ref`](@ref)). All reference constants are at
(T_ref, p_ref, RH_ref) = (260 K, 500 hPa, 100 %).

Fields and defaults follow Williams (2026), Table 1.

References:
- Williams (2026), *J. Adv. Model. Earth Syst.*, doi:10.1029/2025MS005405.
- Armstrong (1968), doi:10.1016/0022-4073(68)90052-6 (diffusivity factor D).
- Mlawer et al. (1997), doi:10.1029/97JD00237 (continuum temperature scaling).

Fields are

$(TYPEDFIELDS)
"""
Base.@kwdef struct WilliamsLongwave{NF} <: AbstractLongwaveScheme
    "Number of evenly spaced wavenumber quadrature points"
    nwavenumber::Int = 41

    "Minimum wavenumber of the spectral integration range [cm⁻¹]"
    wavenumber_min::NF = NF(10)

    "Maximum wavenumber of the spectral integration range [cm⁻¹]"
    wavenumber_max::NF = NF(2500)

    "Peak absorption of the pure-rotation band κ_rot [m² kg⁻¹]"
    κ_rot::NF = NF(37)
    "e-folding decay length of the rotation band l_rot [cm⁻¹]"
    l_rot::NF = NF(56)
    "Peak absorption of the vibration–rotation band κ_vr [m² kg⁻¹]"
    κ_vr::NF = NF(5)
    "e-folding length of vibration–rotation band (low-ν side) l_vr1 [cm⁻¹]"
    l_vr1::NF = NF(37)
    "e-folding length of vibration–rotation band (high-ν side) l_vr2 [cm⁻¹]"
    l_vr2::NF = NF(52)

    "Continuum absorption below 1700 cm⁻¹ κ_cnt1 [m² kg⁻¹]"
    κ_cnt1::NF = NF(0.004)
    "Continuum absorption above 1700 cm⁻¹ κ_cnt2 [m² kg⁻¹]"
    κ_cnt2::NF = NF(0.0002)

    "Include CO₂ longwave absorption"
    do_co2::Bool = true
    "CO₂ concentration [ppmv]"
    co2_ppmv::NF = NF(280)
    "Peak absorption of the CO₂ 15 μm band κ_co2 [m² kg⁻¹]"
    κ_co2::NF = NF(110)
    "e-folding half-width of CO₂ band l_co2 [cm⁻¹]"
    l_co2::NF = NF(12)
    "Centre wavenumber of CO₂ bending mode ν̃_co2 [cm⁻¹]"
    ν̃_co2::NF = NF(667)

    "Two-stream diffusivity factor D (Armstrong 1968)"
    diffusivity::NF = NF(1.5)
    "Reference pressure for pressure broadening [Pa]"
    p_ref::NF = NF(50000)
    "Reference temperature for absorption coefficient fits [K]"
    T_ref::NF = NF(260)
    "Reference saturation water-vapor pressure at T_ref [Pa]"
    pv_ref::NF = NF(224.92)
    "Temperature-scaling exponent for the continuum (Mlawer et al. 1997) [K⁻¹]"
    σ_cont::NF = NF(0.02)
end

Adapt.@adapt_structure WilliamsLongwave

"""$(TYPEDSIGNATURES)
Construct a [`WilliamsLongwave`](@ref) with the given floating-point type.
"""
WilliamsLongwave(::Type{NF}; kwargs...) where NF = WilliamsLongwave{NF}(; kwargs...)

Base.eltype(::WilliamsLongwave{NF}) where NF = NF
Base.eltype(::Type{<:WilliamsLongwave{NF}}) where NF = NF

"""$(TYPEDSIGNATURES)
Column longwave radiative transfer for the Williams (2026) Simple Spectral
Model. Tendencies are accumulated into `dTdt` with `+=`/`-=`; diagnostic
fluxes are written into `diag`.

Sign convention: temperature tendency has units [K s⁻¹]; positive OLR, positive
downward surface flux, positive upward surface flux.
"""
function solve_longwave!(dTdt::AbstractVector,
                         diag::LongwaveDiagnostics{NF},
                         scheme::WilliamsLongwave{NF},
                         profile::ColumnProfile,
                         geometry::ColumnGeometry,
                         surface::ColumnSurface,
                         constants::PhysicalConstants) where NF

    T  = profile.temperature
    q  = profile.humidity
    pₛ = NF(profile.surface_pressure)
    nlayers = length(T)

    σ_SB = NF(constants.stefan_boltzmann)
    cₚ   = NF(constants.heat_capacity)
    g    = NF(constants.gravity)

    ϵ_ocean = NF(surface.ocean_emissivity)
    ϵ_land  = NF(surface.land_emissivity)
    sst     = NF(surface.sea_surface_temperature)
    lst     = NF(surface.land_surface_temperature)
    land_fraction = NF(surface.land_fraction)

    # Broadband Stefan–Boltzmann surface upward flux (for diagnostics).
    U_sfc_ocean = ifelse(isfinite(sst), ϵ_ocean * σ_SB * sst^4, zero(NF))
    U_sfc_land  = ifelse(isfinite(lst), ϵ_land  * σ_SB * lst^4, zero(NF))
    U_sfc_bb    = (1 - land_fraction) * U_sfc_ocean + land_fraction * U_sfc_land

    # Wavenumber quadrature.
    nν = scheme.nwavenumber
    dν̃ = (scheme.wavenumber_max - scheme.wavenumber_min) / NF(nν - 1)

    olr_sum::NF    = zero(NF)
    D_surf_sum::NF = zero(NF)

    for iv in 1:nν
        ν̃ = scheme.wavenumber_min + NF(iv - 1) * dν̃

        B_sfc_ocean = ifelse(isfinite(sst), planck_wavenumber(sst, ν̃), zero(NF))
        B_sfc_land  = ifelse(isfinite(lst), planck_wavenumber(lst, ν̃), zero(NF))

        # Hemispherical (π·B) spectral surface flux [W m⁻²], land-sea weighted.
        U_spec::NF = dν̃ * NF(π) * (
            (1 - land_fraction) * ϵ_ocean * B_sfc_ocean +
             land_fraction      * ϵ_land  * B_sfc_land
        )

        # ---- Upward sweep: k = nlayers → 1 ------------------------------
        U::NF = U_spec
        # Surface upward flux enters the bottom of the lowest layer.
        dTdt[nlayers] += surface_flux_to_tendency(U / cₚ, profile, geometry, constants)

        for k in nlayers:-1:1
            Δτ_k  = williams_delta_tau(k, ν̃, T, q, pₛ, geometry, scheme, g)
            tr_k  = exp(-Δτ_k)
            B_k   = planck_wavenumber(T[k], ν̃)
            U_new::NF = U * tr_k + dν̃ * NF(π) * B_k * (1 - tr_k)

            if k > 1
                # U_new leaves layer k at the top and enters layer k-1 at the bottom.
                dTdt[k]     -= flux_to_tendency(U_new / cₚ, profile, geometry, constants, k)
                dTdt[k - 1] += flux_to_tendency(U_new / cₚ, profile, geometry, constants, k - 1)
            else
                # k == 1: U_new is OLR escaping to space.
                dTdt[1] -= flux_to_tendency(U_new / cₚ, profile, geometry, constants, 1)
                olr_sum += U_new
            end
            U = U_new
        end

        # ---- Downward sweep: k = 1 → nlayers ----------------------------
        # TOA boundary: no longwave incoming from space.
        D::NF = zero(NF)

        for k in 1:(nlayers - 1)
            Δτ_k  = williams_delta_tau(k, ν̃, T, q, pₛ, geometry, scheme, g)
            tr_k  = exp(-Δτ_k)
            B_k   = planck_wavenumber(T[k], ν̃)
            D_new::NF = D * tr_k + dν̃ * NF(π) * B_k * (1 - tr_k)

            dTdt[k]     -= flux_to_tendency(D_new / cₚ, profile, geometry, constants, k)
            dTdt[k + 1] += flux_to_tendency(D_new / cₚ, profile, geometry, constants, k + 1)
            D = D_new
        end

        # Surface-adjacent layer: the downward flux that reaches the surface.
        Δτ_nl = williams_delta_tau(nlayers, ν̃, T, q, pₛ, geometry, scheme, g)
        tr_nl = exp(-Δτ_nl)
        B_nl  = planck_wavenumber(T[nlayers], ν̃)
        D_surf::NF = D * tr_nl + dν̃ * NF(π) * B_nl * (1 - tr_nl)

        dTdt[nlayers] -= surface_flux_to_tendency(D_surf / cₚ, profile, geometry, constants)
        D_surf_sum    += D_surf
    end

    diag.outgoing_longwave        = olr_sum
    diag.surface_longwave_down    = D_surf_sum
    diag.ocean_surface_longwave_up = U_sfc_ocean
    diag.land_surface_longwave_up  = U_sfc_land
    diag.surface_longwave_up       = U_sfc_bb

    return nothing
end
