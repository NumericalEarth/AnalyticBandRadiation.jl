# Reduced ecCKD Gas Pressure-Band Refinement

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Target case | ecckd_rcemip_style_column_subset |
| Target metric | toa_forcing_max_abs |
| Target boundary | toa |
| Target column | 6 |
| Target residual | -2.49623974222 W m^-2 |
| Pressure-band count | 8 |
| Ranked candidate count | 1136 |
| Evaluated candidate count | 96 |
| Evaluated trial count | 768 |
| Base objective | 8.32079914074 |
| Base boundary objective | 8.32079914074 |
| Candidate objective | 8.32070420146 |
| Objective reduction | 9.49392793359e-05 |
| Accepted | true |
| Accepted move count | 1 |
| Best component | static_absorption |
| Best local g-point index | 2 |
| Best official g-point | 4 |
| Best gas index | 1 |
| Best pressure range | 47:53 |
| Best direction | positive |
| Best step | 0.001953125 |

This diagnostic applies one nonnegative gas-specific pressure-band scale
to the current reduced shortwave table and accepts only if the nonlinear
hard-gate objective improves.
