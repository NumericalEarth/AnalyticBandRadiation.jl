# Reduced ecCKD Retained Mixed Component Pareto Scan

Status: **retained_mixed_component_scan_improved**

| Field | Value |
|---|---:|
| Base mode | retained_topology |
| Acceptance rule | bounded_frontier_surface_cap |
| Forcing-error tolerance | 0.01 W m^-2 |
| Surface cap | 2.03 W m^-2 |
| Candidates | 12 |
| Trials | 144 |
| Iterations | 20 / 20 |
| Base objective | 7.19352418454 |
| Final objective | 7.05452306522 |
| Objective reduction | 0.139001119327 |
| Base TOA forcing | 2.15805725536 W m^-2 |
| Base surface forcing | 2.02792974588 W m^-2 |
| Final TOA forcing | 2.11635691956 W m^-2 |
| Final surface forcing | 2.0297810592 W m^-2 |
| Best exact objective | 7.03955734015 |
| Best TOA forcing | 2.11186720204 W m^-2 |
| Best surface forcing | 2.03748904283 W m^-2 |
| Any strict Pareto-safe | true |
| Any surface-cap-safe | true |
| Any tolerance Pareto-safe | true |
| Accepted | true |
| Accepted moves | 20 |

This diagnostic reuses the mixed static/H2O pressure-temperature component
blocks on the fully retained current base. With `strict_pareto`, every
accepted move must reduce the hard objective without increasing TOA or
surface forcing beyond roundoff. With `bounded_frontier`, accepted moves
must still reduce the hard objective while keeping each forcing error
within the configured W m^-2 tolerance of the current iteration base.
If a positive surface cap is configured, accepted moves must also keep
the absolute worst-case surface forcing error at or below that cap.
