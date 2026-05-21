# Reduced ecCKD Boundary-Table Topology Replacement

Status: **topology_replacement_rejected**

| Field | Value |
|---|---:|
| Radius | `2` |
| Evaluated candidates | `33 / 33` |
| Base objective | `7.13985029534` |
| Best objective | `8.91143243322` |
| Best improvement | `-1.77158213788` |
| Best replacement | `g16 -> g18` |
| Worst case | `ecckd_clear_sky_tropical_column` |
| Worst metric | `toa_or_profile_max` |
| Worst value | `2.14556161196` |

This diagnostic replaces one shortwave slot in the current boundary-aware table-continuation model with the raw official ecCKD table for a neighboring g-point. It is intentionally conservative: accepted evidence would still need promotion through `reduced_ecckd_accuracy.jl`.
