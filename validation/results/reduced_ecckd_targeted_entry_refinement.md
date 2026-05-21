# Reduced ecCKD Targeted Entry Refinement

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Target case | ecckd_rcemip_style_column_subset |
| Target metric | toa_forcing_max_abs |
| Target value | 2.58446425485 |
| Target threshold | 0.3 |
| Pressure move count | 4 |
| Targeted candidate count | 1155 |
| Candidate limit | 256 |
| Iterations requested | 2 |
| Accepted move count | 2 |
| Pressure-refined objective | 8.61316394314 |
| Final objective | 8.61303768346 |
| Objective reduction | 0.000126259689353 |
| Improved | true |

This artifact is intentionally narrower than the main optimization preflight. It
only evaluates table entries ranked by the current worst case's interpolation
stencil and gas amounts, then applies the same full hard-objective acceptance
rule as the main active table-entry refinement.
