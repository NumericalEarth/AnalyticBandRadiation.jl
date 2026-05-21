# ecCKD Published Training Manifest

Status: **passed**

ecCKD source root: `/shared/home/greg/.julia/artifacts/7b210aef53e908cfe3c709945f0763c37ca82aaa/ecckd-6115f9b8e29a55cb0f48916857bdc77fec41badd`

Keep these ecCKD source files, script modes, data selections, gases, weights, regularization, and stopping settings fixed; replace only the optimizer implementation/settings in the Julia Reactant/Enzyme recovery pipeline.

## Config

| Field | Value |
|---|---|
| CKDMIP data dir template | `/hugetmp/parr/ckdmip` |
| Training dataset | `evaluation1` |
| Evaluation dataset | `evaluation2` |
| Train with both evaluation sets | `no` |
| MMM dataset | `mmm` |
| Idealized dataset | `idealized` |

## Source Files

| Path | Present |
|---|---:|
| `src/ecckd/optimize_lut.cpp` | true |
| `src/ecckd/solve_adept.cpp` | true |
| `src/ecckd/calc_cost_function_lw.cpp` | true |
| `src/ecckd/calc_cost_function_sw.cpp` | true |
| `test/config.h` | true |
| `test/find_g_points_lw.sh` | true |
| `test/find_g_points_sw.sh` | true |
| `test/create_lut_lw.sh` | true |
| `test/create_lut_sw.sh` | true |
| `test/optimize_lut_lw.sh` | true |
| `test/optimize_lut_sw.sh` | true |
| `test/run_ckd_lw.sh` | true |
| `test/run_ckd_sw.sh` | true |

## Optimization Scripts

### `test/optimize_lut_lw.sh`

- Selected common options: `COMMON_OPTIONS="prior_error=8.0 broadband_weight=0.8 flux_profile_weight=0.2 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 spectral_boundary_weight=0.1"`
- Modes: nwp, zero-minor1, zero-minor2, climate | all-in-one, all1, relative-base, relative-ch4, relative-minor, relative-n2o, relative-cfc11, relative-cfc

### `test/optimize_lut_sw.sh`

- Selected common options: `COMMON_OPTIONS="prior_error=2.0 broadband_weight=0.4 flux_weight=0.4 flux_profile_weight=0.1 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 max_iterations=2000 bounded_optimization=0"`
- Modes: nwp, climate-base, relative-base, relative-ch4, relative-n2o, relative-minor

Remaining external data requirement: Full CKDMIP line-by-line spectral absorption and LBL flux database under RH_CKDMIP_DATA_PATH.
