# Host-Model Access Points Check

Status: **passed**

| Symbol | Exported | Defined |
|---|---:|---:|
| `AbstractAtmosphericState` | true | true |
| `AbstractGasOpticsModel` | true | true |
| `AbstractCloudOpticsModel` | true | true |
| `AbstractAerosolOpticsModel` | true | true |
| `AbstractRadiativeTransferSolver` | true | true |
| `AbstractRadiationBackend` | true | true |
| `ColumnAtmosphere` | true | true |
| `RadiativeFluxes` | true | true |
| `LongwaveOpticalProperties` | true | true |
| `ShortwaveOpticalProperties` | true | true |
| `ShortwaveCloudOverlapOpticalProperties` | true | true |
| `CloudOpticalProperties` | true | true |
| `CloudyRegionCloudOpticalProperties` | true | true |
| `AerosolOpticalProperties` | true | true |
| `CloudScatteringTable` | true | true |
| `EcCKDSpectralMapping` | true | true |
| `EcCKDGasOpticsModel` | true | true |
| `LayerCloudOpticsModel` | true | true |
| `LayerLiquidIceCloudOpticsModel` | true | true |
| `LayerAerosolOpticsModel` | true | true |
| `CloudlessLongwave` | true | true |
| `CloudlessShortwave` | true | true |
| `CloudOverlapShortwave` | true | true |
| `optical_properties!` | true | true |
| `cloud_optical_properties!` | true | true |
| `cloudy_region_optical_properties!` | true | true |
| `aerosol_optical_properties!` | true | true |
| `read_cloud_scattering_table` | true | true |
| `read_ecckd_spectral_mapping` | true | true |
| `cloud_scattering_properties` | true | true |
| `cloud_scattering_gpoint_properties` | true | true |
| `radiative_fluxes!` | true | true |
| `heating_rates!` | true | true |
| `radiative_heating!` | true | true |
| `radiation_workspace` | true | true |
| `add_mapped_cloud_scattering!` | true | true |

## Component Smoke

| Capability | Result |
|---|---:|
| Gas optical properties callable | true |
| Cloud optical properties callable | true |
| Cloudy-region cloud optical properties callable | true |
| Aerosol optical properties callable | true |
| Solvers callable separately | true |
| Heating rates callable separately | true |
| Host can stop after gas optics | true |
| Host can replace solver or vertical integral | true |
