"""
$(TYPEDEF)

Zero-atmosphere shortwave: TOA insolation reaches the surface unattenuated
and is reflected by the surface albedo. Temperature tendencies are zero.
Useful as a baseline and for tests of the surface energy budget.
"""
struct TransparentShortwave <: AbstractShortwaveScheme end

Adapt.@adapt_structure TransparentShortwave

"""$(TYPEDSIGNATURES)
Column transparent-atmosphere shortwave.
"""
function solve_shortwave!(dTdt::AbstractVector,
                          diag::ShortwaveDiagnostics{NF},
                          ::TransparentShortwave,
                          profile::ColumnProfile,
                          geometry::ColumnGeometry,
                          surface::ColumnSurface,
                          constants::PhysicalConstants,
                          thermodynamic::ThermodynamicConstants;
                          cloud_top_convective::Integer = length(profile.temperature) + 1) where NF
    S₀ = NF(constants.solar_constant)
    cos_zenith = NF(surface.cos_zenith)
    D = S₀ * cos_zenith

    diag.surface_shortwave_down       = D
    diag.ocean_surface_shortwave_down = D
    diag.land_surface_shortwave_down  = D

    ocean_up = NF(surface.ocean_albedo) * D
    land_up  = NF(surface.land_albedo)  * D
    albedo   = (1 - NF(surface.land_fraction)) * NF(surface.ocean_albedo) +
               NF(surface.land_fraction) * NF(surface.land_albedo)

    diag.ocean_surface_shortwave_up = ocean_up
    diag.land_surface_shortwave_up  = land_up
    diag.surface_shortwave_up       = albedo * D
    diag.albedo                     = albedo
    diag.outgoing_shortwave         = diag.surface_shortwave_up

    diag.cloud_cover        = zero(NF)
    diag.stratocumulus_cover = zero(NF)
    diag.cloud_top          = length(profile.temperature) + 1
    return nothing
end
