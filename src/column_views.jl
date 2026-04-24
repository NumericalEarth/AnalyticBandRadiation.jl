"""
$(TYPEDEF)

Vertical geometry for a single column expressed in sigma-pressure coordinates.

`σ_full` has length `nlayers` and gives the midpoint of each layer.
`σ_half` has length `nlayers + 1` and gives the layer interfaces.
`σ_thick = diff(σ_half)` has length `nlayers`.
"""
struct ColumnGeometry{NF, V<:AbstractVector{NF}}
    σ_full::V
    σ_half::V
    σ_thick::V
end

function ColumnGeometry(σ_half::AbstractVector{NF}) where NF
    nlayers = length(σ_half) - 1
    σ_full  = @views (σ_half[1:nlayers] .+ σ_half[2:end]) ./ 2
    σ_thick = diff(σ_half)
    return ColumnGeometry{NF, typeof(σ_half)}(σ_full, σ_half, σ_thick)
end

"""
$(TYPEDEF)

Column thermodynamic profile and lower boundary quantities read by the
radiation solvers.

The arrays are indexed top-down: `k = 1` is the top of the atmosphere,
`k = nlayers` is the bottom (surface-adjacent) layer.
"""
struct ColumnProfile{NF, V<:AbstractVector{NF}}
    temperature::V
    humidity::V
    geopotential::V
    surface_pressure::NF
    rain_rate::NF
end

ColumnProfile(; temperature, humidity, geopotential = similar(temperature, 0),
              surface_pressure, rain_rate = zero(eltype(temperature))) =
    ColumnProfile{eltype(temperature), typeof(temperature)}(
        temperature, humidity, geopotential, surface_pressure, rain_rate)

"""
$(TYPEDEF)

Lower-boundary properties needed by the column radiation solvers.
`NaN` in `sea_surface_temperature` or `land_surface_temperature` means the
column has no ocean / no land respectively (the corresponding contribution
is skipped).
"""
Base.@kwdef struct ColumnSurface{NF}
    sea_surface_temperature::NF
    land_surface_temperature::NF
    land_fraction::NF
    ocean_albedo::NF = zero(NF)
    land_albedo::NF = zero(NF)
    cos_zenith::NF = zero(NF)
    ocean_emissivity::NF = one(NF)
    land_emissivity::NF = one(NF)
end

"""
$(TYPEDEF)

Physical constants consumed by the column radiation solvers.
"""
Base.@kwdef struct PhysicalConstants{NF}
    gravity::NF = NF(9.80665)
    heat_capacity::NF = NF(1004.64)
    stefan_boltzmann::NF = NF(5.670374419e-8)
    solar_constant::NF = NF(1361)
end

"""
$(TYPEDEF)

Thermodynamic constants needed for saturation-humidity calculations used by
the diagnostic cloud scheme.
"""
Base.@kwdef struct ThermodynamicConstants{NF}
    saturation_vapor_pressure_reference::NF = NF(610.78)
    latent_heat_condensation::NF = NF(2.501e6)
    gas_constant_vapor::NF = NF(461.50)
    freezing_temperature::NF = NF(273.15)
    molar_mass_ratio::NF = NF(0.622)
end

"""$(TYPEDSIGNATURES)
Sensible Earth defaults for the full set of physical constants needed by the
shortwave solver (constants + thermodynamic constants).
"""
default_earth_constants(::Type{NF}) where NF =
    (physical = PhysicalConstants{NF}(), thermodynamic = ThermodynamicConstants{NF}())

default_earth_constants() = default_earth_constants(Float64)

"""
$(TYPEDEF)

Clausius–Clapeyron saturation specific humidity at `(T, p)` given
`ThermodynamicConstants`. Returns `NaN` if the partial pressure is
unresolvable (e.g. zero total pressure).
"""
@inline function saturation_humidity(T, p, tc::ThermodynamicConstants)
    (; saturation_vapor_pressure_reference, latent_heat_condensation,
       gas_constant_vapor, freezing_temperature, molar_mass_ratio) = tc
    e_sat = saturation_vapor_pressure_reference *
            exp(latent_heat_condensation / gas_constant_vapor *
                (inv(freezing_temperature) - inv(T)))
    return molar_mass_ratio * e_sat / p
end

"""
$(TYPEDEF)

Column longwave diagnostic outputs written by [`solve_longwave!`](@ref).
"""
mutable struct LongwaveDiagnostics{NF}
    outgoing_longwave::NF
    surface_longwave_down::NF
    surface_longwave_up::NF
    ocean_surface_longwave_up::NF
    land_surface_longwave_up::NF
end

LongwaveDiagnostics{NF}() where NF = LongwaveDiagnostics{NF}(
    zero(NF), zero(NF), zero(NF), zero(NF), zero(NF))

"""
$(TYPEDEF)

Column shortwave diagnostic outputs written by [`solve_shortwave!`](@ref).
"""
mutable struct ShortwaveDiagnostics{NF}
    outgoing_shortwave::NF
    surface_shortwave_down::NF
    surface_shortwave_up::NF
    ocean_surface_shortwave_down::NF
    land_surface_shortwave_down::NF
    ocean_surface_shortwave_up::NF
    land_surface_shortwave_up::NF
    albedo::NF
    cloud_cover::NF
    cloud_top::Int
    stratocumulus_cover::NF
end

ShortwaveDiagnostics{NF}(nlayers::Integer = 1) where NF = ShortwaveDiagnostics{NF}(
    zero(NF), zero(NF), zero(NF), zero(NF), zero(NF), zero(NF), zero(NF),
    zero(NF), zero(NF), nlayers + 1, zero(NF))
