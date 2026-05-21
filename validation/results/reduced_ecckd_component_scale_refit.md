# Reduced ecCKD Component-Scale Refit

Status: **component_scale_refit_improved**

| Field | Value |
|---|---:|
| Base mode | boundary_table_post_descent_plus_continuation |
| Parameter count | 48 |
| Probe step | 0.00390625 |
| Max log scale | 0.125 |
| Best ridge lambda | 10000 |
| Coordinate iterations | 12 / 12 |
| Base objective | 7.29585913012 |
| Best linearized exact objective | 8.452155696 |
| Coordinate final objective | 7.13988657025 |
| Final objective | 7.13988657025 |
| Objective reduction | 0.155972559866 |
| Base TOA forcing error | 2.18875773904 W m^-2 |
| Final TOA forcing error | 2.14196597108 W m^-2 |
| Base surface forcing error | 2.18224163998 W m^-2 |
| Final surface forcing error | 2.14149643678 W m^-2 |
| Selected refit | coordinate |
| Accepted moves | 12 |
| Accepted | true |

This diagnostic applies whole-component log-scale factors to the retained
boundary-table 32x16 model. Static absorption, H2O absorption, and
Rayleigh scattering are scaled independently per retained shortwave
g-point, preserving nonnegative tabulated optics.
