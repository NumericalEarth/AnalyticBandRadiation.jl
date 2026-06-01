# Reduced ecCKD Slot-Blend Linearized Refit

Status: **slot_blend_linearized_rejected**

| Field | Value |
|---|---:|
| Base blend count | 7 |
| Candidate count | 32 |
| Single trial count | 1024 |
| Probe alpha | 0.0009765625 |
| Max alpha | 0.0078125 |
| Base objective | 7.47836712724 |
| Best ridge lambda | 100000000 |
| Best exact objective | 10.8272641841 |
| Best objective reduction | -3.3488970569 |
| Best positive delta count | 18 |
| Best TOA forcing error | 2.78645506239 W m^-2 |
| Best surface forcing error | 3.24817925524 W m^-2 |
| Accepted | false |
| Accepted blend count | 0 |
| Total blend count | 7 |

This diagnostic finite-differences ranked slot-blend directions, solves a
bounded nonnegative ridge update over them, then accepts only if the exact
nonlinear hard-gate objective improves.
