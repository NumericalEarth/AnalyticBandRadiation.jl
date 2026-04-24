module AnalyticBandRadiation

using Adapt
using Dates
using DocStringExtensions

export AbstractRadiationScheme, AbstractLongwaveScheme, AbstractShortwaveScheme
export AbstractShortwaveTransmissivity, AbstractShortwaveClouds

export ColumnProfile, ColumnGeometry, ColumnSurface
export PhysicalConstants, ThermodynamicConstants, default_earth_constants

export LongwaveDiagnostics, ShortwaveDiagnostics

export planck_wavenumber

export WilliamsLongwave
export h2o_line_kappa_ref, h2o_cont_kappa_ref, co2_kappa_ref
export williams_delta_tau

export NoClouds, DiagnosticClouds
export ConstantShortwaveTransmissivity, BackgroundShortwaveTransmissivity
export TransparentShortwave, OneBandShortwave, OneBandGreyShortwave
export OneBandShortwaveRadiativeTransfer
export saturation_humidity

export solve_longwave!, solve_shortwave!

export solar_declination, equation_of_time, cosine_solar_zenith

include("abstract_types.jl")
include("column_views.jl")
include("flux_to_tendency.jl")
include("planck.jl")
include("zenith.jl")
include("longwave/williams_absorption.jl")
include("longwave/williams_longwave.jl")
include("shortwave/clouds.jl")
include("shortwave/transmissivity.jl")
include("shortwave/transparent_shortwave.jl")
include("shortwave/one_band_shortwave.jl")

end # module
