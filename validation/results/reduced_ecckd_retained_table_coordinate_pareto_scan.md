# Reduced ecCKD Retained Table Coordinate Pareto Scan

Status: **table_coordinate_pareto_scan_rejected**

| Field | Value |
|---|---:|
| Base objective | 7.19352418454 |
| Base TOA forcing | 2.15805725536 W m^-2 |
| Base surface forcing | 2.02792974588 W m^-2 |
| Candidates | 64 |
| Trials | 640 |
| Best objective | 7.19352418454 |
| Best objective reduction | 0 |
| Best TOA forcing | 2.15805725536 W m^-2 |
| Best surface forcing | 2.02792974527 W m^-2 |
| Best component | static_absorption |
| Best g-point | 30 |
| Best move | negative 0.0001220703125 |
| Best TOA regressed | false |
| Best surface regressed | false |
| Any Pareto-safe | false |
| Best Pareto-safe objective | 7.19352418454 |
| Best Pareto-safe TOA forcing | 2.15805725536 W m^-2 |
| Best Pareto-safe surface forcing | 2.02792974588 W m^-2 |

This diagnostic evaluates exact one-coordinate nonlinear table moves on
the retained boundary-table-continuation base. A candidate is
Pareto-safe only if it lowers the full hard objective and does not worsen
either worst TOA or worst surface forcing error.
