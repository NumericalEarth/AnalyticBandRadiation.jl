# Reduced ecCKD Coefficient Continuation

Status: **coefficient_continuation_above_target**

| Field | Value |
|---|---:|
| Initial objective | 214.264529894 |
| Objective target | 1 |
| Final objective | 8.6288431477 |
| Greedy checkpoint used | true |
| Saved states evaluated | 5 |
| Best start label | preflight_post_coefficient_weight_refinement |
| Best start objective | 8.6288431477 |
| Greedy iterations | 18 |
| Greedy final objective | 9.01318970914 |
| Joint direction best objective | 8.94544944983 |
| Post-joint iterations | 2 |
| Post-joint final objective | 8.6288431477 |
| Weight refinement iterations | 1 |
| Weight refinement final objective | 8.6288431477 |
| Worst case | ecckd_rcemip_style_column_subset |
| Worst metric | toa_forcing_max_abs |
| Worst value | 2.58865294431 |
| Worst threshold | 0.3 |

This artifact intentionally stops before topology, pressure-band table,
and active-entry table scans. It captures only the cheap-to-resume
48-parameter coefficient/weight continuation so useful progress is not
lost when the full preflight becomes too expensive.
