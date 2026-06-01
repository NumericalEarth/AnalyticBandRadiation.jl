# Reduced ecCKD Boundary-Table Topology Weight Refit

Status: **topology_weight_refit_rejected**

| Field | Value |
|---|---:|
| Radius | `2` |
| Weight iterations | `800` |
| Evaluated candidates | `12 / 33` |
| Base objective | `8.7693298333` |
| Best objective | `43.5953802998` |
| Best improvement | `-34.8260504665` |
| Best replacement | `g4 -> g2` |
| Worst case | `ecckd_rcemip_style_column_subset` |
| Worst metric | `surface_or_profile_max` |
| Worst value | `13.0786140899` |

This diagnostic replaces one shortwave slot in the retained boundary-aware
table-continuation model, then refits the reduced shortwave weights with
the hard-gate max-norm weight optimizer before exact objective evaluation.
