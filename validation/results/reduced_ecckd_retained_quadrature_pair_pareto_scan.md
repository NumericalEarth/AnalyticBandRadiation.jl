# Reduced ecCKD Retained Quadrature Pair Pareto Scan

Status: **quadrature_pair_pareto_scan_rejected**

| Field | Value |
|---|---:|
| Base objective | 7.19352418454 |
| Base TOA forcing | 2.15805725536 W m^-2 |
| Base surface forcing | 2.02792974588 W m^-2 |
| Pairs | 120 |
| Candidates | 480 |
| Best objective | 7.79378260905 |
| Best objective reduction | -0.600258424503 |
| Best TOA forcing | 2.33699895857 W m^-2 |
| Best surface forcing | 2.33813478271 W m^-2 |
| Best move | g14 up, g25 down, 0.03125 |
| Best TOA regressed | true |
| Best surface regressed | true |
| Any Pareto-safe | false |
| Best Pareto-safe objective | 7.19352418454 |
| Best Pareto-safe TOA forcing | 2.15805725536 W m^-2 |
| Best Pareto-safe surface forcing | 2.02792974588 W m^-2 |

This diagnostic redistributes reduced shortwave quadrature weight between
two retained g-points at a time. A candidate is Pareto-safe only if it
lowers the full hard objective and does not worsen either worst TOA or
worst surface forcing error.
