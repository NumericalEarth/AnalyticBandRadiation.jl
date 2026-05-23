# ecCKD Published Recovery Vector

Status: **passed**

Reference: `/shared/home/greg/.julia/artifacts/49ce668ce0861f9d5e8299d68af7138485eb5f19/ecrad-131ac980517719b7a859e3ccc117919a1d888a20/data/ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc`

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_vector_roundtrip_sw32_candidate.nc`

| Metric | Value |
|---|---:|
| Arrays | 9 |
| Parameters | 204896 |
| Round-trip max abs error | 0.0 |
| Round-trip L1 relative error | 0.0 |
| Recovery metrics status | passed |
| Worst log-coefficient RMSE | 0.0 |
| G-point weight max abs error | 0.0 |

## Arrays

| Name | Elements | Shape |
|---|---:|---|
| `ch4_molar_absorption_coeff` | 10176 | 32x53x6 |
| `co2_molar_absorption_coeff` | 10176 | 32x53x6 |
| `composite_molar_absorption_coeff` | 10176 | 32x53x6 |
| `gpoint_fraction` | 31840 | 995x32 |
| `h2o_molar_absorption_coeff` | 122112 | 32x53x6x12 |
| `n2o_molar_absorption_coeff` | 10176 | 32x53x6 |
| `o3_molar_absorption_coeff` | 10176 | 32x53x6 |
| `rayleigh_molar_scattering_coeff` | 32 | 32 |
| `solar_irradiance` | 32 | 32 |

This is the optimizer handoff vector for published-model recovery: the optimizer can update this parameter vector, write a CKD-definition candidate, and score it with the published recovery metrics.
