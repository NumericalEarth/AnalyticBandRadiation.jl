# Reduced ecCKD Retained TOA-Probe Expansion

Status: **toa_probe_expansion_rejected**

| Field | Value |
|---|---:|
| Configurations | 2 |
| Best label | toa_probe_32_medium_tiny |
| Best exact objective | 7.13146967612 |
| Best objective reduction | 0.0164856640157 |
| Best TOA forcing | 2.13944090284 W m^-2 |
| Best surface forcing | 2.07939172665 W m^-2 |
| Any accepted | false |

This diagnostic uses the current retained base with a TOA-only
linearized residual and accepts only strict Pareto-safe full-objective
updates that do not regress either TOA or surface forcing.

## Configurations

| Label | Candidates | Probe step | Max log scale | Base objective | Base TOA | Base surface | Best objective | Reduction | TOA | Surface | Pareto safe | Accepted |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| toa_probe_32_medium_tiny | 32 | 0.000244140625 | 0.001953125 | 7.14795534014 | 2.14438660204 | 2.05760969562 | 7.13146967612 | 0.0164856640157 | 2.13944090284 | 2.07939172665 | false | false |
| toa_probe_32_small_tiny | 32 | 0.0001220703125 | 0.0009765625 | 7.14795534014 | 2.14438660204 | 2.05760969562 | 7.13969180337 | 0.0082635367653 | 2.14190754101 | 2.06848449947 | false | false |
