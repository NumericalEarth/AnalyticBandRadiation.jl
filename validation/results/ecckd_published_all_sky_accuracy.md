# Published ecCKD All-Sky Accuracy

Status: **passed**

Each row rewrites package candidate variables into the matched all-sky ecRad reference using the same Tripleclouds/aerosol configuration as the current all-sky IFS gate, but with model-specific ecCKD gas-optics and cloud-scattering mapping files.

| Model | LW | SW | Passed | TOA forcing | Surface forcing | LW TOA | LW surface | SW TOA | SW surface | Hard objective | Limiting metric |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| official ecCKD 1.0 32-LW x 32-SW all-sky climate model | 32 | 32 | true | 0.116188256566 W m^-2 | 0.0128536244544 W m^-2 | 0.111466467535 | 0.00646439036603 | 0.00687483939146 | 0.00764450172153 | 0.387294188553 | `toa_forcing_abs_error` |
| official ecCKD 1.0/1.2 32-LW x 64-SW all-sky climate/window model | 32 | 64 | true | 0.231391178591 W m^-2 | 0.17169319049 W m^-2 | 0.111466467535 | 0.00646439036603 | 0.23098837544 | 0.171384584585 | 0.771303928637 | `toa_forcing_abs_error` |
| official ecCKD 1.0/1.4 32-LW x 96-SW all-sky climate/vfine model | 32 | 96 | true | 0.116158872897 W m^-2 | 0.0129957468097 W m^-2 | 0.111466467535 | 0.00646439036603 | 0.00681068525694 | 0.00774535229959 | 0.387196242991 | `toa_forcing_abs_error` |
| official ecCKD 1.2/1.4 64-LW x 32-SW all-sky narrow/rgb model | 64 | 32 | true | 0.11171179402 W m^-2 | 0.0126801196332 W m^-2 | 0.106990004989 | 0.00603254364034 | 0.00687483939146 | 0.00764450172153 | 0.372372646733 | `toa_forcing_abs_error` |
| official ecCKD 1.2 64-LW x 64-SW all-sky climate model | 64 | 64 | true | 0.231170230722 W m^-2 | 0.171755120473 W m^-2 | 0.106990004989 | 0.00603254364034 | 0.23098837544 | 0.171384584585 | 0.770567435739 | `toa_forcing_abs_error` |
| official ecCKD 1.2/1.4 64-LW x 96-SW all-sky climate/vfine model | 64 | 96 | true | 0.111682410351 W m^-2 | 0.0127809702112 W m^-2 | 0.106990004989 | 0.00603254364034 | 0.00681068525694 | 0.00774535229959 | 0.372274701172 | `toa_forcing_abs_error` |

## Candidate Configuration

| Environment variable | Value |
|---|---|
| `RH_AEROSOL_OPTICS` | `true` |
| `RH_CANDIDATE_GAS_OPTICS` | `official_ecckd` |
| `RH_CLOUD_FRACTION_EXPONENT` | `1.0` |
| `RH_CLOUD_INHOM_OVERLAP_EXPONENT` | `2.0` |
| `RH_CLOUD_OVERLAP_LONGWAVE` | `true` |
| `RH_CLOUD_OVERLAP_LONGWAVE_RULE` | `tripleclouds_alpha` |
| `RH_CLOUD_OVERLAP_RULE` | `tripleclouds_alpha` |
| `RH_CLOUD_OVERLAP_SHORTWAVE` | `true` |
| `RH_CLOUD_OVERLAP_USE_CLOUDY_REGION_OPTICS` | `true` |
| `RH_CLOUD_SCATTERING_TABLE_OPTICS` | `true` |
| `RH_IFS_AEROSOL_TABLE_OPTICS` | `true` |
| `RH_LW_CLOUD_SCATTERING` | `true` |
