# API reference

## Column inputs

```@docs
ColumnProfile
ColumnGeometry
ColumnSurface
PhysicalConstants
ThermodynamicConstants
default_earth_constants
LongwaveDiagnostics
ShortwaveDiagnostics
```

## Longwave

```@docs
WilliamsLongwave
solve_longwave!
planck_wavenumber
h2o_line_kappa_ref
h2o_cont_kappa_ref
co2_kappa_ref
AnalyticBandRadiation.williams_delta_tau
```

## Shortwave

```@docs
TransparentShortwave
OneBandShortwave
OneBandGreyShortwave
OneBandShortwaveRadiativeTransfer
AnalyticBandRadiation.AbstractShortwaveTransmissivity
ConstantShortwaveTransmissivity
BackgroundShortwaveTransmissivity
AnalyticBandRadiation.AbstractShortwaveClouds
NoClouds
DiagnosticClouds
solve_shortwave!
```

## Solar geometry

```@docs
solar_declination
equation_of_time
cosine_solar_zenith
```
