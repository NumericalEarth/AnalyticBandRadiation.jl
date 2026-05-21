# Reduced ecCKD Pair Slot-Blend Refinement

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Radius | 31 |
| Base blend count | 7 |
| Ranked single count | 12 |
| Evaluated pair count | 48 |
| Base objective | 7.47836712724 |
| Best single objective | 8.39717935543 |
| Best pair objective | 9.10566044494 |
| Final objective | 7.47836712724 |
| Objective reduction | 0 |
| Accepted | false |
| Accepted blend count | 0 |
| Total blend count | 7 |
| Best first move | g12 -> g18 alpha 0.001953125 |
| Best second move | g21 -> g20 alpha 0.0078125 |
| Best worst case | ecckd_clear_sky_tropical_column |
| Best worst metric | surface_forcing_max_abs |

This diagnostic evaluates bounded pairs from the best single slot-blend
candidates and accepts the pair only if the nonlinear hard-gate objective
improves.
