"""
$(TYPEDEF)

Trivial cloud model that reports zero cloud cover. Returned shape matches
[`DiagnosticClouds`](@ref) so downstream code can be agnostic.
"""
struct NoClouds <: AbstractShortwaveClouds end

Adapt.@adapt_structure NoClouds

"""
$(TYPEDEF)

Diagnostic clouds after the Fortran SPEEDY scheme (Kucharski, Molteni &
Bracco, 2006). Cloud cover is a combination of a relative-humidity term and
a precipitation term; the highest layer exceeding the RH threshold sets the
cloud top. An independent stratocumulus term is diagnosed at the surface
from dry-static-energy stability.

Fields are

$(TYPEDFIELDS)
"""
Base.@kwdef struct DiagnosticClouds{NF} <: AbstractShortwaveClouds
    "Relative humidity threshold for cloud cover = 0 [1]"
    relative_humidity_threshold_min::NF = NF(0.3)
    "Relative humidity threshold for cloud cover = 1 [1]"
    relative_humidity_threshold_max::NF = NF(1)
    "Specific humidity threshold for cloud cover [kg/kg]"
    specific_humidity_threshold_min::NF = NF(0.0002)
    "Weight for the √precipitation term [1]"
    precipitation_weight::NF = NF(0.2)
    "Cap on precipitation contributing to cloud cover [mm/day]"
    precipitation_max::NF = NF(10)
    "Cloud albedo at CLC = 1 [1]"
    cloud_albedo::NF = NF(0.6)
    "Stratocumulus cloud albedo [1]"
    stratocumulus_albedo::NF = NF(0.5)
    "Static-stability lower threshold for stratocumulus (GSES0) [J/kg]"
    stratocumulus_stability_min::NF = NF(0.25)
    "Static-stability upper threshold for stratocumulus (GSES1) [J/kg]"
    stratocumulus_stability_max::NF = NF(0.4)
    "Maximum stratocumulus cloud cover (CLSMAX) [1]"
    stratocumulus_cover_max::NF = NF(0.6)
    "Enable the stratocumulus parameterization"
    use_stratocumulus::Bool = true
    "Stratocumulus cloud factor (SPEEDY clfact) [1]"
    stratocumulus_cloud_factor::NF = NF(1.2)
end

Adapt.@adapt_structure DiagnosticClouds

DiagnosticClouds(::Type{NF}; kwargs...) where NF = DiagnosticClouds{NF}(; kwargs...)

"""$(TYPEDSIGNATURES)
Returned tuple from cloud diagnosis: `(cloud_cover, cloud_top, cloud_albedo,
stratocumulus_cover, stratocumulus_albedo)`.

For `NoClouds`, `cloud_top = nlayers + 1` (below the surface) so downstream
shortwave code skips the cloud-reflection branch.
"""
@inline function diagnose_clouds(::NoClouds, profile::ColumnProfile,
                                  geometry::ColumnGeometry,
                                  surface::ColumnSurface,
                                  constants::PhysicalConstants,
                                  thermodynamic::ThermodynamicConstants,
                                  cloud_top_convective::Integer)
    NF = eltype(profile.temperature)
    nlayers = length(profile.temperature)
    return (
        cloud_cover = zero(NF),
        cloud_top = nlayers + 1,
        cloud_albedo = zero(NF),
        stratocumulus_cover = zero(NF),
        stratocumulus_albedo = zero(NF),
    )
end

@inline function diagnose_clouds(clouds::DiagnosticClouds{NF},
                                  profile::ColumnProfile,
                                  geometry::ColumnGeometry,
                                  surface::ColumnSurface,
                                  constants::PhysicalConstants,
                                  thermodynamic::ThermodynamicConstants,
                                  cloud_top_convective::Integer) where NF

    T = profile.temperature
    q = profile.humidity
    Φ = profile.geopotential
    nlayers = length(T)
    pₛ = profile.surface_pressure
    σ_full = geometry.σ_full
    cₚ = constants.heat_capacity
    land_fraction = surface.land_fraction

    rh_min = clouds.relative_humidity_threshold_min
    rh_max = clouds.relative_humidity_threshold_max
    q_min  = clouds.specific_humidity_threshold_min
    precip_weight = clouds.precipitation_weight
    precip_max    = clouds.precipitation_max

    # Precipitation term — rain_rate in m/s; convert to mm/day.
    precip_term = min(precip_max, (NF(86400) * profile.rain_rate) / NF(1000))
    P = precip_weight * sqrt(max(zero(NF), precip_term))

    humidity_term::NF = zero(NF)
    cloud_top_humidity = nlayers + 1

    for k in 1:(nlayers - 1)
        q_k = q[k]
        qsat = saturation_humidity(T[k], σ_full[k] * pₛ, thermodynamic)
        if q_k > q_min && qsat > 0
            rh_k = q_k / qsat
            if rh_k >= rh_min
                rh_norm = max(zero(NF), (rh_k - rh_min) / (rh_max - rh_min))
                humidity_term = min(one(NF), rh_norm)^2
                cloud_top_humidity = min(k, cloud_top_humidity)
            end
        end
    end

    cloud_cover = min(one(NF), P + humidity_term)
    cloud_top   = min(cloud_top_humidity, cloud_top_convective)

    stratocumulus_cover::NF = zero(NF)
    if clouds.use_stratocumulus
        surface_k = nlayers
        above_k   = max(1, nlayers - 1)
        G = (cₚ * T[surface_k] + Φ[surface_k]) - (cₚ * T[above_k] + Φ[above_k])
        stab_min = clouds.stratocumulus_stability_min
        stab_max = clouds.stratocumulus_stability_max
        static_stability = clamp((G - stab_min) / (stab_max - stab_min), zero(NF), one(NF))
        cover_max = clouds.stratocumulus_cover_max
        cloud_factor = clouds.stratocumulus_cloud_factor
        strc_ocean = static_stability * max(cover_max - cloud_factor * cloud_cover, zero(NF))
        qsat_surf = saturation_humidity(T[surface_k], σ_full[surface_k] * pₛ, thermodynamic)
        rh_surf   = q[surface_k] / qsat_surf
        strc_land = strc_ocean * rh_surf
        stratocumulus_cover = (1 - land_fraction) * strc_ocean + land_fraction * strc_land
    end

    return (
        cloud_cover = cloud_cover,
        cloud_top = cloud_top,
        cloud_albedo = clouds.cloud_albedo,
        stratocumulus_cover = stratocumulus_cover,
        stratocumulus_albedo = clouds.stratocumulus_albedo,
    )
end
