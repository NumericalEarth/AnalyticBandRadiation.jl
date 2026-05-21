# Reduced ecCKD Hard-Gate Sparse Subset Search

Status: **failed_threshold**

Method: `overlapping sparse nonnegative 16-term shortwave flux-basis hard-gate subset search`

| Field | Value |
|---|---:|
| Selected shortwave g-points | `1, 2, 9, 10, 12, 13, 14, 16, 21, 22, 25, 27, 28, 30, 31, 32` |
| Exact objective | `89.1338250789` |
| Approximate objective | `59.3838504434` |
| Worst case | `ecckd_clear_sky_tropical_column` |
| Worst TOA forcing error | `7.39985344849` W m^-2 |
| Worst surface forcing error | `16.536096435` W m^-2 |
| Passed hard thresholds | `false` |

This is a flux-basis feasibility diagnostic. It uses nonnegative sparse weights over official shortwave g-point responses, then exact-checks the selected support through the normal tabulated model path. It does not replace the physical reduced coefficient-table artifact unless it passes and is promoted into `reduced_ecckd_accuracy.jl`.

## Search History

| Pass | Approximate objective | Exact objective | Indices |
|---:|---:|---:|---|
| 0 | 253.1282683 | 214.213692302 | `1, 4, 9, 10, 12, 13, 14, 16, 21, 22, 25, 27, 28, 30, 31, 32` |
| 1 | 59.3838504434 | 89.1338250789 | `1, 2, 9, 10, 12, 13, 14, 16, 21, 22, 25, 27, 28, 30, 31, 32` |
