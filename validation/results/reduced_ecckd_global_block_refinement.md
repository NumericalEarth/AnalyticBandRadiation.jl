# Reduced ecCKD Global Block Refinement

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Target case | ecckd_rcemip_style_column_subset |
| Target metric | toa_forcing_max_abs |
| Target value | 2.56687278214 |
| Target threshold | 0.3 |
| Pressure move count | 4 |
| Starting active move count | 23 |
| Group limit | 16 |
| Block sizes | 2, 4, 8 |
| Iterations requested | 2 |
| Accepted move count | 10 |
| Initial objective | 8.55624260713 |
| Final objective | 8.55304192199 |
| Objective reduction | 0.00320068514242 |
| Improved | true |

This artifact tests coherent multi-entry table moves rather than one
scalar entry at a time. Candidate entries are grouped by component,
local shortwave slot, official g-point, and gas identity, then accepted
only when the nonlinear full hard-gate objective improves.
