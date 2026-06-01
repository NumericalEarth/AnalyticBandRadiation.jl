# ecCKD Candidate Transfer Smoke

Status: **candidate_transfer_smoke_ready**

## Blockers

None.

## Transfer Metrics

| Metric | Value |
|---|---:|
| Layers | 54 |
| LW g-points | 32 |
| SW g-points | 32 |
| LW optical-depth max | 0.5868147812958674 |
| SW optical-depth max | 0.14584953550654836 |
| SW Rayleigh optical-depth max | 0.20746887619244184 |
| LW original-objective broadband loss | 22144.533068625602 |
| SW original-objective broadband loss | 1947.8125470141588 |
| LW spectral-projection bands | 13 |
| SW spectral-projection bands | 13 |
| LW spectral-projection original-objective loss | 356.58059037574515 |
| SW spectral-projection original-objective loss | 330.40047827264766 |
| LW spectral-projection non-finite flux count | 0 |
| SW spectral-projection non-finite flux count | 0 |

This is a cloudless transfer smoke test through package gas optics and solvers. It now includes a g-point-to-CKDMIP-band spectral projection score, but full recovery still requires the differentiable ecCKD/equivalent spectral transfer used by the original objective.
