# Reduced ecCKD Retained Post-Capped Weight Refit

Status: **retained_post_capped_weight_refit_improved**

| Field | Value |
|---|---:|
| Iterations | 25000 |
| Base objective | 6.98655990975 |
| Refit objective | 6.97926430984 |
| Objective reduction | 0.00729559991206 |
| Approximate hard-gate objective | 7.67854251783 |
| Base TOA forcing | 2.09245154437 W m^-2 |
| Base surface forcing | 2.01344999846 W m^-2 |
| Refit TOA forcing | 2.09206312139 W m^-2 |
| Refit surface forcing | 2.00649036013 W m^-2 |
| Accepted | true |
| Max absolute weight delta | 2.81506440697e-06 |

This diagnostic refits only the retained 16 shortwave quadrature weights
after composing the retained capped-table optimizer and continuation
moves. It is accepted only when the exact hard objective falls and
neither worst TOA nor worst surface forcing error regresses.
