# Reduced ecCKD Retained Capped Table Optimizer

Status: **retained_capped_table_optimizer_improved**

| Field | Value |
|---|---:|
| Base mode | retained_topology |
| Candidate scope | all_global_residual_probe |
| Residual mode | surface |
| Surface cap | 2.03 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Candidate count | 48 |
| Probe step | 0.001953125 |
| Max log scale | 0.0078125 |
| Base objective | 7.05452306522 |
| Base TOA forcing | 2.11635691956 W m^-2 |
| Base surface forcing | 2.0297810592 W m^-2 |
| Cap-safe row present | true |
| Best cap-safe objective | 7.00134059822 |
| Best cap-safe TOA forcing | 2.10040217947 W m^-2 |
| Best cap-safe surface forcing | 2.01397583862 W m^-2 |
| Best overall objective | 7.00134059822 |
| Best overall TOA forcing | 2.10040217947 W m^-2 |
| Best overall surface forcing | 2.01397583862 W m^-2 |
| Best unsafe objective | 7.15865669544 |
| Best unsafe TOA forcing | 2.13014303985 W m^-2 |
| Best unsafe surface forcing | 2.14759700863 W m^-2 |
| Accepted | true |
| Accepted moves | 48 |

This diagnostic reruns the constrained table-entry optimizer on the
current retained base, but composes only rows that improve the exact
objective, stay within the configured TOA tolerance, and remain below
the absolute surface forcing cap.
