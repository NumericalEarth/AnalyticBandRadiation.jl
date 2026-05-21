# Reduced ecCKD Retained Post-Weight Bounded Weight Refit

Status: **retained_post_weight_bounded_weight_refit_improved**

| Field | Value |
|---|---:|
| Iterations | 50000 |
| Surface cap | 2.03 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Base objective | 6.97926430984 |
| Refit objective | 6.97212754254 |
| Objective reduction | 0.00713676729187 |
| Approximate hard-gate objective | 7.67476528352 |
| Base TOA forcing | 2.09206312139 W m^-2 |
| Base surface forcing | 2.00649036013 W m^-2 |
| Refit TOA forcing | 2.09163826276 W m^-2 |
| Refit surface forcing | 2.00694307171 W m^-2 |
| Accepted | true |
| Max absolute weight delta | 3.42765050393e-06 |

This diagnostic refits shortwave quadrature weights after the post-weight
surface-table update. It accepts objective and TOA improvements while
enforcing the configured absolute surface-forcing cap.
