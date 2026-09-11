module NumericalRadiationSpeedyWeatherExt

using NumericalRadiation
using SpeedyWeather
using Adapt

import NumericalRadiation: AtmosphereProfile, ColumnGrid, SurfaceState,
    PhysicalConstants, LongwaveDiagnostics, solve_longwave!, AnalyticBandLongwave

# -----------------------------------------------------------------------------
# Longwave adapter
# -----------------------------------------------------------------------------

"""
    SpeedyAnalyticBandLongwave{NF} <: SpeedyWeather.AbstractLongwave

SpeedyWeather wrapper around [`NumericalRadiation.AnalyticBandLongwave`](@ref).
Allows for setting the default CO₂ concentration [ppmv], used when the model has
no `greenhouse_gases` component with a `co2` entry.

Usage:

```julia
spectral_grid = SpectralGrid()
longwave = SpeedyAnalyticBandLongwave(spectral_grid; CO₂ = 280)
model = PrimitiveWetModel(spectral_grid; radiation = Radiation(spectral_grid; longwave))
```
"""
struct SpeedyAnalyticBandLongwave{NF} <: SpeedyWeather.AbstractLongwave
    scheme::AnalyticBandLongwave{NF}
    default_CO₂::NF
end

Adapt.@adapt_structure SpeedyAnalyticBandLongwave

function SpeedyAnalyticBandLongwave(SG::SpeedyWeather.SpectralGrid; CO₂ = 280, kwargs...)
    return SpeedyAnalyticBandLongwave(AnalyticBandLongwave{SG.NF}(; kwargs...), SG.NF(CO₂))
end

SpeedyWeather.initialize!(::SpeedyAnalyticBandLongwave, ::SpeedyWeather.PrimitiveEquation) = nothing

# Re-export under the PR's original name for drop-in compatibility.
const SimpleSpectralLongwave = SpeedyAnalyticBandLongwave

@inline function speedy_physical_constants(model)
    NF = typeof(model.planet.gravity)
    return PhysicalConstants{NF}(
        gravity          = model.planet.gravity,
        heat_capacity    = model.atmosphere.heat_capacity,
        stefan_boltzmann = model.atmosphere.stefan_boltzmann,
        solar_constant   = model.planet.solar_constant,
    )
end

@inline function speedy_column_geometry(model)
    geom = model.geometry
    return ColumnGrid(geom.σ_levels_full, geom.σ_levels_half, geom.σ_levels_thick)
end

Base.@propagate_inbounds function SpeedyWeather.parameterization!(ij, vars,
                                                                   rad::SpeedyAnalyticBandLongwave{NF},
                                                                   model) where NF
    time_stepping = model.time_stepping
    T_all    = SpeedyWeather.get_prognostic_step(vars.grid.temperature, time_stepping, rad)
    q_all    = SpeedyWeather.get_prognostic_step(vars.grid.humidity, time_stepping, rad)
    dTdt_all = SpeedyWeather.get_tendency_step(vars.tendencies.grid.temperature, time_stepping, rad)
    sst_all  = SpeedyWeather.get_prognostic_step(vars.prognostic.ocean.sea_surface_temperature,
                                                 time_stepping, rad)

    T    = @view T_all[ij, :]
    q    = @view q_all[ij, :]
    Φ    = @view vars.dynamics.geopotential[ij, :]
    dTdt = @view dTdt_all[ij, :]
    pₛ   = vars.parameterizations.surface_pressure[ij]            # [Pa]

    CO₂ = let prog = vars.prognostic
        if hasproperty(prog, :greenhouse_gases) && haskey(prog.greenhouse_gases, :co2)
            NF(prog.greenhouse_gases.co2[])
        else
            rad.default_CO₂
        end
    end

    profile  = AtmosphereProfile(temperature = T, humidity = q,
                                 geopotential = Φ, surface_pressure = pₛ,
                                 CO₂ = CO₂)

    geometry = speedy_column_geometry(model)
    surface  = SurfaceState{NF}(
        sea_surface_temperature  = sst_all[ij],
        land_surface_temperature = vars.prognostic.land.soil_temperature[ij, 1],
        land_fraction            = model.land_sea_mask.land_fraction[ij],
    )
    constants = speedy_physical_constants(model)
    diag = LongwaveDiagnostics{NF}()

    solve_longwave!(dTdt, diag, rad.scheme, profile, geometry, surface, constants)

    vars.parameterizations.outgoing_longwave[ij]         = diag.outgoing_longwave
    vars.parameterizations.surface_longwave_down[ij]     = diag.surface_longwave_down
    vars.parameterizations.surface_longwave_up[ij]       = diag.surface_longwave_up
    vars.parameterizations.ocean.surface_longwave_up[ij] = diag.ocean_surface_longwave_up
    vars.parameterizations.land.surface_longwave_up[ij]  = diag.land_surface_longwave_up

    return nothing
end

end # module
