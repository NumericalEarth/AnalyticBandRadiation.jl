# Reduced ecCKD Grouped Quadrature Weight Refit

Status: **weight_refit_improved**

| Field | Value |
|---|---:|
| Candidates | 10 |
| Iterations | 3000 |
| Best label | weight_descending_flux0.1_heat1.0_boundary10.0 |
| Best base objective | 16.0106691529 |
| Best refit objective | 15.5327466353 |
| Best objective reduction | 0.477922517621 |
| Best approximate objective | 15.5153673809 |
| Worst TOA forcing | 2.44473581761 W m^-2 |
| Worst surface forcing | 4.65982399058 W m^-2 |
| Max absolute weight delta | 3.20299444989e-05 |
| Passed hard thresholds | false |

Best groups: `1+3, 2+9, 4+12, 5+26, 6+8, 7+10, 11+29, 13+14, 15+17, 16+18, 19+31, 20+21, 22+30, 23+32, 24+25, 27+28`

This diagnostic keeps each grouped quadrature coefficient table fixed but
refits its 16 nonnegative shortwave weights on the hard-gate-scaled
flux, heating-rate, and boundary residual basis before exact evaluation.

## Candidate Summary

| Label | Base objective | Refit objective | TOA forcing | Surface forcing | Passed |
|---|---:|---:|---:|---:|---:|
| weight_descending_flux0.1_heat1.0_boundary10.0 | 16.0106691529 | 15.5327466353 | 2.44473581761 | 4.65982399058 | false |
| weight_descending_flux1.0_heat1.0_boundary1.0 | 19.8565998546 | 19.3667724914 | 3.46405736437 | 5.81003174742 | false |
| spectral_flux0.25_heat1.0_boundary4.0 | 42.5755033263 | 42.3686044554 | 12.7105813366 | 7.41704210201 | false |
| spectral_flux1.0_heat1.0_boundary1.0 | 42.576877057 | 42.3699451984 | 12.7109835595 | 7.47493195972 | false |
| weight_ascending_flux0.25_heat1.0_boundary4.0 | 42.9269746681 | 42.5850528742 | 1.71930889153 | 10.1545268649 | false |
| weight_ascending_flux1.0_heat1.0_boundary1.0 | 42.9685216609 | 42.6268931469 | 1.71743331301 | 10.1544954498 | false |
| weight_ascending_flux0.1_heat1.0_boundary10.0 | 43.1007654187 | 42.7609717409 | 1.71640194253 | 10.1134447362 | false |
| reverse_spectral_flux1.0_heat1.0_boundary1.0 | 74.4029217309 | 73.9210152729 | 5.77346185813 | 22.1763045819 | false |
| weighted_greedy_first_flux1.0_heat1.0_boundary1.0 | 193.83098222 | 193.270702763 | 42.6400845766 | 57.9812108289 | false |
| adjacent_spectral_pairs | 549.944116676 | 549.612107199 | 66.0275557097 | 164.88363216 | false |
