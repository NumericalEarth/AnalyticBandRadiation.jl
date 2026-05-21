# Reduced ecCKD Global Entry Refinement

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Target case | ecckd_clear_sky_tropical_column |
| Target metric | toa_forcing_max_abs |
| Target value | 2.56760882716 |
| Target threshold | 0.3 |
| Pressure move count | 4 |
| Starting active move count | 14 |
| Ranked candidate count | 22000 |
| Candidate limit | 256 |
| Candidate offset | 3328 |
| Iterations requested | 1 |
| Sweep windows | 8 |
| Sweep stride | 256 |
| Sweep max no-improve | 2 |
| Accepted move count | 7 |
| Initial objective | 8.55869609052 |
| Final objective | 8.5582182932 |
| Objective reduction | 0.000477797321992 |
| Improved | true |

This artifact broadens the active-entry search beyond entries inside the
accepted pressure-band moves. It ranks all reduced shortwave slots by the
current worst case's interpolation stencil and gas amounts, then applies
the same full hard-objective acceptance rule as the main refinement.
