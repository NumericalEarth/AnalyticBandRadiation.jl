module AnalyticBandRadiationSpeedyWeatherExt

using AnalyticBandRadiation
using SpeedyWeather
using Adapt

import AnalyticBandRadiation: ColumnProfile, ColumnGeometry, ColumnSurface,
    PhysicalConstants, ThermodynamicConstants, LongwaveDiagnostics,
    ShortwaveDiagnostics, solve_longwave!, solve_shortwave!,
    WilliamsLongwave, TransparentShortwave, OneBandShortwave, OneBandGreyShortwave

# -----------------------------------------------------------------------------
# Longwave adapter
# -----------------------------------------------------------------------------

"""
    SpeedyWilliamsLongwave{NF} <: SpeedyWeather.AbstractLongwave

SpeedyWeather-native wrapper around [`AnalyticBandRadiation.WilliamsLongwave`](@ref).
Exposes the same per-column dispatch as the PR #1057 prototype so that a
`PrimitiveWetModel` can be constructed with

    model = PrimitiveWetModel(spectral_grid;
                              longwave_radiation = SpeedyWilliamsLongwave(spectral_grid))
"""
struct SpeedyWilliamsLongwave{NF} <: SpeedyWeather.AbstractLongwave
    scheme::WilliamsLongwave{NF}
end

Adapt.@adapt_structure SpeedyWilliamsLongwave

function SpeedyWilliamsLongwave(SG::SpeedyWeather.SpectralGrid; kwargs...)
    return SpeedyWilliamsLongwave(WilliamsLongwave{SG.NF}(; kwargs...))
end

SpeedyWeather.initialize!(::SpeedyWilliamsLongwave, ::SpeedyWeather.PrimitiveEquation) = nothing

# Re-export under the PR's original name for drop-in compatibility.
const SimpleSpectralLongwave = SpeedyWilliamsLongwave

@inline function _speedy_physical_constants(model)
    NF = typeof(model.planet.gravity)
    return PhysicalConstants{NF}(
        gravity          = model.planet.gravity,
        heat_capacity    = model.atmosphere.heat_capacity,
        stefan_boltzmann = model.atmosphere.stefan_boltzmann,
        solar_constant   = model.planet.solar_constant,
    )
end

@inline function _speedy_column_geometry(model)
    geom = model.geometry
    return ColumnGeometry(geom.σ_levels_full, geom.σ_levels_half, geom.σ_levels_thick)
end

function SpeedyWeather.parameterization!(ij::Integer, vars,
                                         rad::SpeedyWilliamsLongwave{NF},
                                         model) where NF
    nlayers = size(vars.grid.temperature_prev, 2)

    T  = @view vars.grid.temperature_prev[ij, :]
    q  = @view vars.grid.humidity_prev[ij, :]
    Φ  = vars.grid.geopotential isa Nothing ?
         zeros(eltype(T), nlayers) :
         @view vars.grid.geopotential[ij, :]
    pₛ = vars.grid.pressure_prev[ij]

    profile  = ColumnProfile(temperature = T, humidity = q,
                             geopotential = Φ, surface_pressure = pₛ)
    geometry = _speedy_column_geometry(model)
    surface  = ColumnSurface{NF}(
        sea_surface_temperature  = vars.prognostic.ocean.sea_surface_temperature[ij],
        land_surface_temperature = vars.prognostic.land.soil_temperature[ij, 1],
        land_fraction            = model.land_sea_mask.mask[ij],
    )
    constants = _speedy_physical_constants(model)
    diag = LongwaveDiagnostics{NF}()
    dTdt = @view vars.tendencies.grid.temperature[ij, :]

    solve_longwave!(dTdt, diag, rad.scheme, profile, geometry, surface, constants)

    vars.parameterizations.outgoing_longwave[ij]        = diag.outgoing_longwave
    vars.parameterizations.surface_longwave_down[ij]    = diag.surface_longwave_down
    vars.parameterizations.surface_longwave_up[ij]      = diag.surface_longwave_up
    vars.parameterizations.ocean.surface_longwave_up[ij] = diag.ocean_surface_longwave_up
    vars.parameterizations.land.surface_longwave_up[ij]  = diag.land_surface_longwave_up

    return nothing
end

end # module
