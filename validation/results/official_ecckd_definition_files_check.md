# Official ecCKD Definition Files Check

Status: **passed**

This check verifies that the package reader recognizes the official ecCKD definition files shipped in the ecRad checkout. It does not yet convert those definitions into production optical-depth kernels.

| Kind | Status | G-points | Bands | Required gases present | Source tables | Rayleigh tables | Path |
|---|---|---:|---:|---:|---:|---:|---|
| longwave | passed | 64 | 13 | true | true | false | `ecckd-1.2_lw_climate_narrow-64b_ckd-definition.nc` |
| shortwave | passed | 32 | 5 | true | false | true | `ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc` |

## Runtime LUT Ingestion

Status: **passed**

| Field | Value |
|---|---|
| Gas names | `composite, h2o, o3, co2, ch4, n2o, cfc11, cfc12` |
| H2O mole fraction sample | 0.005 |
| Pressure grid size | 53 |
| Temperature grid shape | `[53, 6]` |
| Longwave absorption shape | `[64, 8, 53, 6]` |
| Shortwave absorption shape | `[32, 8, 53, 6]` |
| Shortwave Rayleigh coefficients | 32 |
| Shortwave Rayleigh coefficients positive | true |
| Shortwave Rayleigh optical-depth max | 1.009856558801884 |
| Longwave source table present | true |
| Longwave source table shape | `[64, 231]` |
| Longwave source table finite | true |
| Finite optical depths | true |
| Nonnegative optical depths | true |
| Weights normalized | true |
