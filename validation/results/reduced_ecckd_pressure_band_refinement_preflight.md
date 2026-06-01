# Reduced ecCKD Pressure-Band Refinement Preflight

Status: **pressure_band_move_improved**

| Field | Value |
|---|---:|
| Current full objective | 8.61488084949 |
| Target case | ecckd_rcemip_style_column_subset |
| Target metric | toa_forcing_max_abs |
| Target metric objective | 8.61488084949 |
| Pressure bands | 4 |
| Candidates | 16 |
| Best target-metric band | 1 |
| Best target-metric scale | 1.13315 |
| Best target-metric objective | 8.614963167 |
| Best full-objective band | 1 |
| Best full-objective scale | 0.882497 |
| Best full objective | 9.43100405638 |
| Global pressure-band accepted | false |
| Per-g pressure-band candidates | 64 |
| Best per-g target g-point | 4 |
| Best per-g target band | 3 |
| Best per-g target objective | 8.61442114049 |
| Best per-g full g-point | 4 |
| Best per-g full band | 3 |
| Best per-g full objective | 8.61442114049 |
| Per-g pressure-band accepted | true |
| Component per-g candidates | 32 |
| Best component per-g full component | static_absorption |
| Best component per-g full g-point | 4 |
| Best component per-g full band | 3 |
| Best component per-g full objective | 8.61442089472 |
| Component per-g pressure-band accepted | true |
| Iterative component per-g iterations completed | 2 |
| Iterative component per-g accepted moves | 2 |
| Iterative component per-g final objective | 8.61439354943 |
| Active table-entry candidates | 32 |
| Active table-entry component | static_absorption |
| Active table-entry g-point | 4 |
| Best active table-entry full objective | 8.61486190339 |
| Active table-entry accepted | true |
| Iterative active table-entry iterations completed | 2 |
| Iterative active table-entry accepted moves | 2 |
| Iterative active table-entry final objective | 8.61483757625 |
| Pairwise per-g candidates | 13 |
| Pairwise selected single candidates | 6 |
| Best pairwise full objective | 8.61442158026 |
| Pairwise pressure-band accepted | true |
| Iterative per-g iterations requested | 2 |
| Iterative per-g iterations completed | 2 |
| Iterative per-g accepted moves | 2 |
| Iterative per-g final objective | 8.61439378607 |
| Iterative per-g objective reduction | 0.000487063424638 |

Next required work: The accepted iterative component pressure-band moves are already promoted into the main reduced optimizer and reduced-accuracy artifact; move beyond bounded pressure-band table scales to a constrained multi-parameter coefficient-table or quadrature-bin optimizer against flux and heating residuals.
