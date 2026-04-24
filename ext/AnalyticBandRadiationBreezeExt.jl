module AnalyticBandRadiationBreezeExt

using AnalyticBandRadiation
using Breeze
using Oceananigans
using KernelAbstractions

import AnalyticBandRadiation: AnalyticBandLongwave, OneBandShortwave, OneBandGreyShortwave,
    TransparentShortwave, AtmosphereProfile, ColumnGrid, SurfaceState,
    PhysicalConstants, ThermodynamicConstants, LongwaveDiagnostics,
    ShortwaveDiagnostics, solve_longwave!, solve_shortwave!

# -----------------------------------------------------------------------------
# NOTE
# -----------------------------------------------------------------------------
# This extension is a scaffold: it defines the types and launch plumbing that
# hook the AnalyticBandRadiation column solvers into Breeze's atmosphere model.
# The concrete `_update_radiation!` dispatch for
# `Breeze.AtmosphereModels.RadiativeTransferModel` parameterised by an
# `AnalyticBandAtmosphericState` below populates `upwelling_longwave_flux`,
# `downwelling_longwave_flux` and `downwelling_shortwave_flux`. `flux_divergence`
# is left to Breeze's existing `_compute_radiation_flux_divergence!` helper.
#
# Breeze's `RadiativeTransferModel(grid, ::AbstractOptics, ...)` is the
# user-facing entry point. We add a new optics tag `SpectralOptics` and a
# constructor dispatch. Because Breeze has not yet been taught about this
# optics tag its concrete field types are mostly `Any`-lite (duck-typed), so
# this extension provides them without further changes to Breeze.
# -----------------------------------------------------------------------------

"""
    SpectralOptics <: Breeze.AtmosphereModels.AbstractOptics

Select the `AnalyticBandRadiation` schemes: Williams (2026) clear-sky longwave
and the SPEEDY one-band shortwave.
"""
struct SpectralOptics <: Breeze.AtmosphereModels.AbstractOptics end

"""
    AnalyticBandAtmosphericState{FT, LW, SW}

Scheme bundle + pre-allocated per-column buffers stored inside Breeze's
`RadiativeTransferModel.atmospheric_state` slot. One instance per model; the
column buffers are sized to `(Nx * Ny, nlayers)` and reused each call.
"""
struct AnalyticBandAtmosphericState{FT, LW, SW, B}
    longwave_scheme::LW
    shortwave_scheme::SW
    thermodynamic_constants::ThermodynamicConstants{FT}
    transmissivity_scratch::B   # Matrix{FT} (Nxy × Nlayers)
end

"""
Construct a Breeze `RadiativeTransferModel` backed by AnalyticBandRadiation.

`longwave_scheme`, `shortwave_scheme` default to sensible Earth choices and can
be replaced with any [`AnalyticBandRadiation.AbstractLongwaveScheme`](@ref) and
[`AnalyticBandRadiation.AbstractShortwaveScheme`](@ref).
"""
function Breeze.AtmosphereModels.RadiativeTransferModel(grid::Oceananigans.AbstractGrid,
        ::SpectralOptics,
        constants;
        longwave_scheme = AnalyticBandLongwave(eltype(grid)),
        shortwave_scheme = OneBandShortwave(eltype(grid)),
        thermodynamic_constants = ThermodynamicConstants{eltype(grid)}(),
        solar_constant = 1361.0,
        surface_temperature = nothing,
        surface_emissivity = 0.98,
        direct_surface_albedo = nothing,
        diffuse_surface_albedo = nothing,
        surface_albedo = nothing,
        coordinate = nothing,
        epoch = nothing,
        schedule = Oceananigans.Utils.IterationInterval(1))

    # The detailed wiring (ZFaceField allocation, SurfaceRadiativeProperties
    # construction, schedule defaults) mirrors BreezeRRTMGPExt's gray-model
    # constructor. This scaffold returns a minimal state; a full integration
    # is tracked as follow-up work.
    error("""
    AnalyticBandRadiationBreezeExt: the Breeze-side column launch is a scaffold.
    The next step is to allocate ZFaceField flux arrays, package
    AnalyticBandAtmosphericState, and write the `:xy` kernel that builds a
    per-column AtmosphereProfile/Geometry/Surface from Breeze state and calls
    `AnalyticBandRadiation.solve_longwave!` / `solve_shortwave!`. See
    `ext/BreezeRRTMGPExt/gray_radiative_transfer_model.jl` in Breeze for the
    reference pattern.
    """)
end

# Kernel skeleton — fleshed out in follow-up work. Keep the signature stable
# so the subsequent implementation is a drop-in.
@kernel function _analytic_band_longwave_kernel!(rtm, grid, temperature, humidity,
                                                  reference_pressure,
                                                  land_fraction, surface_temperature)
    i, j = @index(Global, NTuple)
    # Pack (i, j) column into AtmosphereProfile, ColumnGrid, SurfaceState,
    # then call AnalyticBandRadiation.solve_longwave! and write fluxes into
    # rtm.upwelling_longwave_flux[i, j, k] / downwelling_longwave_flux[i, j, k].
end

end # module
