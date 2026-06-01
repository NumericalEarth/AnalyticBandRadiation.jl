# Reduced ecCKD Slot-Blend Refinement

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Radius | 31 |
| Iterations requested | 1 |
| Iterations completed | 1 |
| Candidate count | 256 |
| Evaluated trial count | 1024 |
| Active move count | 417 |
| Base objective | 7.47836712724 |
| Candidate objective | 9.18814316568 |
| Final objective | 7.47836712724 |
| Objective reduction | 0 |
| Accepted | false |
| Accepted blend count | 0 |
| Total blend count | 7 |
| Best local index | 9 |
| Best move | g21 -> g18 |
| Best alpha | 0.25 |
| Best worst case | ecckd_clear_sky_tropical_column |
| Best worst metric | surface_forcing_max_abs |

This diagnostic convexly blends one optimized reduced shortwave slot
toward a nearby full official ecCKD g-point table. It is a bounded
quadrature/table deformation and is accepted only if the nonlinear
hard-gate objective improves.
