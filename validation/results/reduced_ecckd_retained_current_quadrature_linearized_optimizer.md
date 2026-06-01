# Reduced ecCKD Retained Current Quadrature Linearized Optimizer

Status: **current_quadrature_linearized_optimizer_rejected**

| Field | Value |
|---|---:|
| Base mode | current_capped_post_weight |
| Residual mode | all_shortwave |
| Basis count | 15 |
| Probe step | 0.00390625 |
| Max logit delta | 0.03125 |
| Surface cap | 2.03 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Base objective | 6.97212754254 |
| Base TOA forcing | 2.09163826276 W m^-2 |
| Base surface forcing | 2.00694307171 W m^-2 |
| Best ridge lambda | 10000 |
| Best exact objective | 9.74757110026 |
| Best objective reduction | -2.77544355771 |
| Best TOA forcing | 2.2812210358 W m^-2 |
| Best surface forcing | 2.92427133008 W m^-2 |
| Accepted | false |
| Accepted ridge lambda | n/a |
| Accepted objective | 6.97212754254 |
| Accepted TOA forcing | 2.09163826276 W m^-2 |
| Accepted surface forcing | 2.00694307171 W m^-2 |

This diagnostic fits a linearized all-logit shortwave quadrature-weight
update on the current capped/post-weight base. It exact-evaluates every
ridge candidate and accepts only objective/TOA improvements that remain
inside the configured absolute surface-forcing cap.

## Ridge Summary

| Ridge lambda | Objective | TOA forcing | Surface forcing | Accepted |
|---:|---:|---:|---:|---:|
| 1e-06 | 13.3841256544 | 2.17744766106 | 4.01523769631 | false |
| 0.0001 | 13.3841269894 | 2.17744768941 | 4.01523809683 | false |
| 0.01 | 13.3842604266 | 2.17745052315 | 4.01527812798 | false |
| 1 | 13.3973246902 | 2.17772795946 | 4.01919740707 | false |
| 100 | 13.9756858058 | 2.19000435711 | 4.19270574174 | false |
| 10000 | 9.74757110026 | 2.2812210358 | 2.92427133008 | false |
