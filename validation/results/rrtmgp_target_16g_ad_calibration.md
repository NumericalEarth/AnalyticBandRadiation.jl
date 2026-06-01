# RRTMGP-Target 16-g AD Calibration

Status: `passed`

Reference model: RRTMGP clear-sky fluxes through AnalyticBandRadiationRRTMGPExt

Candidate model: 16-g EcCKDGasOpticsModel shortwave absorption/weight prototype

Loss: mean squared shortwave flux residual plus scaled shortwave heating-rate residual against package-native RRTMGP reference

Parameterization: 16 softmax shortwave weights, 16 H2O absorption log-scales, 16 CO2 absorption log-scales

Training method: `Enzyme reverse-mode gradient descent`

Initial loss: 34711.9494945

Final loss: 34703.9349766

Reactant used for training: `true`

Enzyme used for training: `true`

AD calibration evidence: Reactant compiled the RRTMGP-target 16-g loss and Enzyme reverse-mode gradients drove a decreasing gradient-descent calibration of the 16-g model parameters.

## Final Metrics

- Flux RMSE: 277.987326399 W m^-2
- Flux max abs: 453.859291359 W m^-2
- Heating-rate RMSE: 8.02946963642e-05 K s^-1
- Heating-rate max abs: 0.000172384816897 K s^-1
- TOA forcing error: 569.42820632 W m^-2
- Surface forcing error: 630.794993633 W m^-2

## Acceptance Rule

This artifact may pass only when Reactant and Enzyme are used to differentiate the RRTMGP-target loss and gradient descent trains the 16-g model parameters.
