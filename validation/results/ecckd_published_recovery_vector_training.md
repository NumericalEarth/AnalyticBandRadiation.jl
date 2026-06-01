# ecCKD Published Recovery Vector Training

Status: **passed**

Reference: `/shared/home/greg/.julia/artifacts/49ce668ce0861f9d5e8299d68af7138485eb5f19/ecrad-131ac980517719b7a859e3ccc117919a1d888a20/data/ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc`

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_vector_trained_sw32_candidate.nc`

| Metric | Value |
|---|---:|
| Arrays | 9 |
| Parameters | 204896 |
| Trained parameters | 64 |
| Optimizer | analytic quadratic gradient descent in log vector space |
| Iterations | 36 |
| Initial loss | 0.00048209024593346344 |
| Final loss | 1.6326386893276442e-17 |
| Loss reduction factor | 2.1711414519448967e12 |
| Enzyme requested | false |
| Enzyme used for training | false |
| Reactant check requested | false |
| Reactant compile check | skipped |
| Round-trip max abs error | 1.0075558278993014e-24 |
| Round-trip L1 relative error | 3.592347169292998e-29 |
| Recovery metrics status | passed |
| Worst log-coefficient RMSE | 1.6004069750833284e-10 |
| Worst p99 relative coefficient error | 2.4424906541753444e-15 |
| G-point weight max abs error | 0.0 |

This trains a deterministic slice of the executable published-model handoff vector back to the SW32 published target while keeping the target arrays and metric definitions fixed, then writes a complete CKD-definition candidate. The default fast path uses an analytic gradient; Enzyme and Reactant checks are optional kwargs because the full AD compiler path is already covered elsewhere and is too slow for this vector handoff test. It is still a vector-recovery gate, not the full original CKDMIP flux objective.
