# CKDMIP Original Objective Dataset

Status: **dataset_samples_ready**

CKDMIP data root: `/shared/home/greg/data/ckdmip`

Training flux schemas ready: 52 / 52

## Blockers

None.

## Representative Samples

| Kind | Scenario | Column | mu0 index | Layers | Bands | Sum layer weight | Self loss | Spectral boundary |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| longwave | `rel-415` | 1 |  | 54 | 13 | 1.0 | 0.0 | false |
| shortwave | `rel-415` | 1 | 1 | 54 | 13 | 1.0 | 0.0 | true |

This artifact proves the original-objective recovery path can read the CKDMIP LBL flux products, reproduce ecCKD layer weights, compute K s^-1 heating targets from flux divergence, and feed those arrays into the Julia CKD loss assembly.
