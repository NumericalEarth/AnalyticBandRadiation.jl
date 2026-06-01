# Reduced ecCKD Retained Current Component-Scale Optimizer

Status: **current_component_scale_optimizer_improved**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Basis | per_gpoint_static_h2o_rayleigh_component_scales |
| Basis count | 48 |
| Probe step | 0.0009765625 |
| Max log scale | 0.001953125 |
| Surface cap | 2.03 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 0.001 |
| Base objective | 6.97212754254 |
| Base TOA forcing | 2.09163826276 W m^-2 |
| Base surface forcing | 2.00694307171 W m^-2 |
| Base heating RMSE | 0.348575997036 K day^-1 |
| Best objective | 6.96933451127 |
| Best objective reduction | 0.00279303127267 |
| Best TOA forcing | 2.09080035338 W m^-2 |
| Best surface forcing | 2.01498794456 W m^-2 |
| Best heating RMSE | 0.347361493009 K day^-1 |
| Accepted | true |
| Accepted ridge lambda | 0.0001 |
| Accepted objective | 6.96933451127 |
| Accepted TOA forcing | 2.09080035338 W m^-2 |
| Accepted surface forcing | 2.01498794456 W m^-2 |
| Accepted heating RMSE | 0.347361493009 K day^-1 |

This diagnostic tests whether the component-scale parameter family that
helped earlier in the chain has any remaining cap-safe leverage on the
fully composed current base. The residual and acceptance guards match
the current heating-profile diagnostics.

## Ridge Summary

| Ridge lambda | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted |
|---:|---:|---:|---:|---:|---:|
| 1e-08 | 6.9761778822 | 2.09285336466 | 2.0190413847 | 0.346908303526 | false |
| 1e-06 | 6.9713510972 | 2.09140532916 | 2.01495187743 | 0.347400671763 | false |
| 0.0001 | 6.96933451127 | 2.09080035338 | 2.01498794456 | 0.347361493009 | true |
| 0.01 | 6.97133827408 | 2.09140148222 | 2.01624187984 | 0.34740708106 | false |
| 1 | 6.97132235646 | 2.09139670694 | 2.01712429666 | 0.347405342225 | false |
| 100 | 6.97605788408 | 2.09281736522 | 2.00963042938 | 0.347490749942 | false |
| 10000 | 6.97705522769 | 2.09311656831 | 2.0119604924 | 0.347268755764 | false |
| 1000000 | 7.01911343862 | 2.10573403159 | 2.01637435761 | 0.348157114283 | false |
