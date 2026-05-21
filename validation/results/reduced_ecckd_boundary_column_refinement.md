# Reduced ecCKD Boundary-Column Refinement

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Target case | ecckd_rcemip_style_column_subset |
| Target metric | toa_forcing_max_abs |
| Target boundary | toa |
| Target column | 6 |
| Target residual | -2.49645664034 W m^-2 |
| Pressure move count | 4 |
| Starting active move count | 408 |
| Ranked candidate count | 18528 |
| Evaluated candidate count | 128 |
| Evaluated trial count | 1024 |
| Base objective | 8.32152213447 |
| Base boundary objective | 8.32152213447 |
| Candidate objective | 8.3214840192 |
| Objective reduction | 3.81152726447e-05 |
| Accepted | true |
| Accepted move count | 1 |
| Best direction | positive |
| Best step | 0.001953125 |

This diagnostic scans active coefficient-table entries that are actually
used by the column producing the current worst boundary forcing error.
It accepts only a move that reduces the full normalized hard-gate
objective across the reduced validation cases.
