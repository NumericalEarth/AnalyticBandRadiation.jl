# Reduced ecCKD Boundary-Table Coordinate Descent

Status: **coordinate_descent_improved**

| Field | Value |
|---|---:|
| Base mode | boundary_weight_refit_plus_retained_boundary_table_moves |
| Include Rayleigh candidates | true |
| Candidates | 48 |
| Iteration limit | 6 |
| Completed iterations | 6 |
| Baseline objective | 7.29587581499 |
| Final objective | 7.29585929551 |
| Final objective reduction | 1.65194763708e-05 |
| Baseline TOA forcing error | 2.1887627445 W m^-2 |
| Final TOA forcing error | 2.18875778865 W m^-2 |
| Baseline surface forcing error | 2.18229856753 W m^-2 |
| Final surface forcing error | 2.18224163998 W m^-2 |
| Accepted | true |

This exact coordinate descent starts from the retained non-recursive
boundary-table state, including accepted single and pair Rayleigh moves,
then accumulates further single active-entry moves only when the full
hard objective improves.
