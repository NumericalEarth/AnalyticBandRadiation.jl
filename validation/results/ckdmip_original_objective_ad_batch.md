# CKDMIP Original Objective Optimizer Batch

Status: **optimizer_batch_ready**

CKDMIP data root: `/shared/home/greg/data/ckdmip`

## Blockers

None.

## Optimizer Probe

| Metric | Value |
|---|---:|
| Parameter count | 32 |
| Initial loss | 23.664279526707595 |
| Final loss after one accepted step | 14.767512921270171 |
| Loss reduction factor | 1.60245531206701 |
| Gradient norm | 5632.184439709027 |
| Finite-difference relative error | 0.0 |

This artifact uses real CKDMIP LW/SW rel-415 training samples and a compact original-objective flux-correction batch. It is an optimizer-readiness probe, not a published-model recovery; Enzyme coverage remains in `test/test_ecckd_original_objective_loss.jl`.
