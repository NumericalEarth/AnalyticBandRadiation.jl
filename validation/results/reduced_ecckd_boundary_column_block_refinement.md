# Reduced ecCKD Boundary-Column Block Refinement

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Target case | ecckd_rcemip_style_column_subset |
| Target metric | toa_forcing_max_abs |
| Target boundary | toa |
| Target column | 6 |
| Target residual | -2.49642847799 W m^-2 |
| Starting active move count | 413 |
| Ranked candidate count | 18528 |
| Block count | 48 |
| Evaluated trial count | 288 |
| Base objective | 8.32142825996 |
| Base boundary objective | 8.32142825996 |
| Candidate objective | 8.32133115454 |
| Objective reduction | 9.71054268941e-05 |
| Accepted | true |
| Accepted move count | 4 |
| Best local g-point index | 13 |
| Best component | static_absorption |
| Best block size | 4 |
| Best direction | positive |
| Best step | 0.001953125 |

This diagnostic groups entries used by the current worst boundary column
and applies one coherent nonnegative scaling to each block. It accepts
only candidates that reduce the full normalized hard-gate objective.
