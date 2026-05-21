# Reduced ecCKD Grouped Quadrature Search

Status: **failed_threshold**

| Field | Value |
|---|---:|
| Candidates | 10 |
| ng_lw | 32 |
| ng_sw | 16 |
| Best label | weight_descending_flux0.1_heat1.0_boundary10.0 |
| Best objective | 16.0106691529 |
| Worst TOA forcing | 2.43386896118 W m^-2 |
| Worst surface forcing | 4.80320074587 W m^-2 |
| Passed hard thresholds | false |

Best groups: `1+3, 2+9, 4+12, 5+26, 6+8, 7+10, 11+29, 13+14, 15+17, 16+18, 19+31, 20+21, 22+30, 23+32, 24+25, 27+28`

This diagnostic searches deterministic 16-bin pairings generated from
different flux, heating-rate, and boundary-flux feature scalings and
anchor orderings. Each candidate is exact-evaluated as a grouped
spectral-weighted coefficient table against the hard clean ecCKD cases.

## Candidate Summary

| Label | Objective | TOA forcing | Surface forcing | Passed |
|---|---:|---:|---:|---:|
| weight_descending_flux0.1_heat1.0_boundary10.0 | 16.0106691529 | 2.43386896118 | 4.80320074587 | false |
| weight_descending_flux1.0_heat1.0_boundary1.0 | 19.8565998546 | 3.53906574306 | 5.95697995639 | false |
| spectral_flux0.25_heat1.0_boundary4.0 | 42.5755033263 | 12.7726509979 | 7.51742967466 | false |
| spectral_flux1.0_heat1.0_boundary1.0 | 42.576877057 | 12.7730631171 | 7.57564996682 | false |
| weight_ascending_flux0.25_heat1.0_boundary4.0 | 42.9269746681 | 1.70943877339 | 10.2784037487 | false |
| weight_ascending_flux1.0_heat1.0_boundary1.0 | 42.9685216609 | 1.70755769394 | 10.2784037487 | false |
| weight_ascending_flux0.1_heat1.0_boundary10.0 | 43.1007654187 | 1.70660518539 | 10.2362750979 | false |
| reverse_spectral_flux1.0_heat1.0_boundary1.0 | 74.4029217309 | 5.81073867668 | 22.3208765193 | false |
| weighted_greedy_first_flux1.0_heat1.0_boundary1.0 | 193.83098222 | 42.7357534419 | 58.1492946661 | false |
| adjacent_spectral_pairs | 549.944116676 | 66.0484885039 | 164.983235003 | false |
