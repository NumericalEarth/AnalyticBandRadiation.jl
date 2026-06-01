# Reduced ecCKD Boundary-Table Triple Coordinate Scan

Status: **triple_coordinate_scan_improved**

| Field | Value |
|---|---:|
| Base mode | boundary_weight_refit_plus_retained_boundary_table_and_coordinate_descent_moves |
| Include Rayleigh candidates | true |
| Single trials | 384 |
| Selected singles | 16 |
| Triple trials | 560 |
| Base objective | 7.29585929551 |
| Best objective | 7.29585928897 |
| Best objective reduction | 6.53993037503e-09 |
| Base TOA forcing error | 2.18875778865 W m^-2 |
| Best TOA forcing error | 2.18875778669 W m^-2 |
| Base surface forcing error | 2.18224163998 W m^-2 |
| Best surface forcing error | 2.18224163998 W m^-2 |
| Accepted | true |

This exact triple scan starts from the fully retained non-recursive
boundary-table state, ranks exact single-coordinate moves, and evaluates
all triples from the retained single-move set against the full hard
objective.
