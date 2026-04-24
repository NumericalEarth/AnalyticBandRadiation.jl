"""
    AbstractRadiationScheme

Root type for column radiation schemes defined by `AnalyticBandRadiation`.
"""
abstract type AbstractRadiationScheme end

"""
    AbstractLongwaveScheme <: AbstractRadiationScheme

A column longwave radiative transfer scheme. Concrete subtypes implement
[`solve_longwave!`](@ref).
"""
abstract type AbstractLongwaveScheme <: AbstractRadiationScheme end

"""
    AbstractShortwaveScheme <: AbstractRadiationScheme

A column shortwave radiative transfer scheme. Concrete subtypes implement
[`solve_shortwave!`](@ref).
"""
abstract type AbstractShortwaveScheme <: AbstractRadiationScheme end

"""
    AbstractShortwaveTransmissivity <: AbstractRadiationScheme

A layer-transmissivity model used by shortwave radiative-transfer solvers.
"""
abstract type AbstractShortwaveTransmissivity <: AbstractRadiationScheme end

"""
    AbstractShortwaveClouds <: AbstractRadiationScheme

A diagnostic cloud model for the shortwave.
"""
abstract type AbstractShortwaveClouds <: AbstractRadiationScheme end
