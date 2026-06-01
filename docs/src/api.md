# API reference

## Umbrella

```@docs
RadiativeTransferColumn
reset!
radiative_heating!
radiation_workspace
```

## Staged runtime interface

```@docs
AbstractAtmosphericState
AbstractGasOpticsModel
AbstractCloudOpticsModel
AbstractAerosolOpticsModel
AbstractRadiativeTransferSolver
AbstractRadiationBackend
ColumnAtmosphere
RadiativeFluxes
LongwaveOpticalProperties
CloudlessLongwave
LongwaveBoundaryConditions
ShortwaveOpticalProperties
CloudlessShortwave
ShortwaveBoundaryConditions
CloudOpticalProperties
LayerCloudOpticsModel
add_cloud_optical_depths!
AerosolOpticalProperties
LayerAerosolOpticsModel
add_aerosol_optical_depths!
optical_properties!
cloud_optical_properties!
aerosol_optical_properties!
radiative_fluxes!
heating_rates!
```

## ecCKD Schema

```@docs
EcCKDDefinition
EcCKDSchemaSummary
EcCKDGasOpticsModel
EcCKDTabulatedGasOpticsModel
read_ecckd_definition
summarize_ecckd_definition
validate_ecckd_definition
```

## Validation metrics

```@docs
RadiationErrorMetrics
RadiationThresholds
radiation_error_metrics
radiative_flux_error_metrics
passes_thresholds
```

## Column inputs

```@docs
AtmosphereProfile
ColumnGrid
SurfaceState
PhysicalConstants
ThermodynamicConstants
default_earth_constants
LongwaveDiagnostics
ShortwaveDiagnostics
```

## Longwave

```@docs
AnalyticBandLongwave
solve_longwave!
planck_wavenumber
h2o_line_kappa_ref
h2o_cont_kappa_ref
co2_kappa_ref
NumericalRadiation.williams_delta_tau
```

## Shortwave

```@docs
NumericalRadiation.TransparentShortwave
NumericalRadiation.OneBandShortwave
NumericalRadiation.OneBandGreyShortwave
OneBandShortwaveRadiativeTransfer
NumericalRadiation.AbstractShortwaveTransmissivity
ConstantShortwaveTransmissivity
BackgroundShortwaveTransmissivity
NumericalRadiation.AbstractShortwaveClouds
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
