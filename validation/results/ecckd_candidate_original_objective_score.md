# ecCKD Candidate Original-Objective Score

Status: **candidate_objective_score_ready**

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_vector_trained_sw32_candidate.nc`

Candidate recovery metrics status: `passed`

## Blockers

None.

## Sample Objective Scores

| Kind | Layers | Bands | Zero-forward loss | Perturbed-forward loss | Perturbation increases loss |
|---|---:|---:|---:|---:|---:|
| longwave | 54 | 13 | 0.0 | 0.0003010888973765814 | true |
| shortwave | 54 | 13 | 0.0 | 0.0010536975080752567 | true |

This scores supplied forward flux/heating arrays against CKDMIP samples. The remaining heavy step is replacing the identity/perturbed forward arrays with differentiable transfer computed from the CKD-definition candidate.
