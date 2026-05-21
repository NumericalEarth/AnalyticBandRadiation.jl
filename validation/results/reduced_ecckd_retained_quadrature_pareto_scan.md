# Reduced ecCKD Retained Quadrature Pareto Scan

Status: **quadrature_pareto_scan_rejected**

| Field | Value |
|---|---:|
| Base objective | 7.19352418454 |
| Base TOA forcing | 2.15805725536 W m^-2 |
| Base surface forcing | 2.02792974588 W m^-2 |
| Candidates | 128 |
| Best objective | 7.84146264025 |
| Best objective reduction | -0.647938455706 |
| Best TOA forcing | 2.34126285103 W m^-2 |
| Best surface forcing | 2.35243879207 W m^-2 |
| Best move | g14 positive 0.03125 |
| Best TOA regressed | true |
| Best surface regressed | true |
| Any Pareto-safe | false |
| Best Pareto-safe objective | 7.19352418454 |
| Best Pareto-safe TOA forcing | 2.15805725536 W m^-2 |
| Best Pareto-safe surface forcing | 2.02792974588 W m^-2 |

This diagnostic perturbs one reduced shortwave quadrature logit at a time
on the fully retained boundary-table-continuation base. A candidate is
Pareto-safe only if it lowers the full hard objective and does not worsen
either worst TOA or worst surface forcing error.
