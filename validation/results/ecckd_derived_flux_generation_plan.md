# ecCKD Derived Flux Generation Plan

Status: **derived_flux_generation_required**

CKDMIP data root: `/shared/home/greg/data/ckdmip`

ecCKD source root: `/shared/home/greg/.julia/artifacts/7b210aef53e908cfe3c709945f0763c37ca82aaa/ecckd-6115f9b8e29a55cb0f48916857bdc77fec41badd`

The 5gas-* and rel-* flux products are generated ecCKD training targets, not public CKDMIP archive files. Generate them in a writable ecCKD working copy, then rerun the CKDMIP preflight.

## Progress

- Expected derived flux products: 18
- Final products present: 3
- Products with raw chunks present: 1
- Raw chunks present: 3/90
- Completed-equivalent raw chunks: 18/90
- Observed raw chunk rate: `8.046269021356558` chunks/hour
- Estimated raw chunk hours remaining: `10.81246472981241`

## Missing Derived Products

Missing derived flux products: 15

| Path | Domain | Scenario | Script |
|---|---|---|---|
| `evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_5gas-1120.h5` | `lw` | `5gas-1120` | `test/run_lw_lbl_evaluation.sh` |
| `evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_5gas-180.h5` | `lw` | `5gas-180` | `test/run_lw_lbl_evaluation.sh` |
| `evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_5gas-2240.h5` | `lw` | `5gas-2240` | `test/run_lw_lbl_evaluation.sh` |
| `evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_5gas-280.h5` | `lw` | `5gas-280` | `test/run_lw_lbl_evaluation.sh` |
| `evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_5gas-415.h5` | `lw` | `5gas-415` | `test/run_lw_lbl_evaluation.sh` |
| `evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_5gas-560.h5` | `lw` | `5gas-560` | `test/run_lw_lbl_evaluation.sh` |
| `evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_rel-1120.h5` | `lw` | `rel-1120` | `test/run_lw_lbl_evaluation.sh` |
| `evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_rel-2240.h5` | `lw` | `rel-2240` | `test/run_lw_lbl_evaluation.sh` |
| `evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_rel-560.h5` | `lw` | `rel-560` | `test/run_lw_lbl_evaluation.sh` |
| `evaluation1/sw_fluxes/ckdmip_evaluation1_sw_fluxes_rel-1120.h5` | `sw` | `rel-1120` | `test/run_sw_lbl_evaluation.sh` |
| `evaluation1/sw_fluxes/ckdmip_evaluation1_sw_fluxes_rel-180.h5` | `sw` | `rel-180` | `test/run_sw_lbl_evaluation.sh` |
| `evaluation1/sw_fluxes/ckdmip_evaluation1_sw_fluxes_rel-2240.h5` | `sw` | `rel-2240` | `test/run_sw_lbl_evaluation.sh` |
| `evaluation1/sw_fluxes/ckdmip_evaluation1_sw_fluxes_rel-280.h5` | `sw` | `rel-280` | `test/run_sw_lbl_evaluation.sh` |
| `evaluation1/sw_fluxes/ckdmip_evaluation1_sw_fluxes_rel-415.h5` | `sw` | `rel-415` | `test/run_sw_lbl_evaluation.sh` |
| `evaluation1/sw_fluxes/ckdmip_evaluation1_sw_fluxes_rel-560.h5` | `sw` | `rel-560` | `test/run_sw_lbl_evaluation.sh` |

## Raw Chunk Progress

| Path | Final present | Raw chunks | Missing chunks | Raw bytes |
|---|---:|---:|---:|---:|
| `evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_rel-560.h5` | false | 3/5 | 2 | 258835 |

## Required ecCKD Scripts

| Path | Present |
|---|---:|
| `test/run_lw_lbl_evaluation.sh` | true |
| `test/run_sw_lbl_evaluation.sh` | true |
| `test/copy_to_ckdmip_lw.sh` | true |
| `test/copy_to_ckdmip_sw.sh` | true |
| `test/config.h` | true |

## Concatenation Tool

- `ncrcat` present: true
- `ncrcat` path: `/shared/home/greg/.local/bin/ncrcat`
- Julia concat shim: true
- Note: ncrcat resolves to the Julia concat shim used for CKDMIP raw chunk assembly.

## Scenario Batches

- LW scenarios: `5gas-180 5gas-280 5gas-415 5gas-560 5gas-1120 5gas-2240 rel-180 rel-280 rel-415 rel-560 rel-1120 rel-2240`
- SW scenarios: `rel-180 rel-280 rel-415 rel-560 rel-1120 rel-2240`

## Suggested Working-Copy Commands

```sh
RH_ECCKD_DERIVED_FLUX_DRY_RUN=true RH_CKDMIP_DATA_PATH=/shared/home/greg/data/ckdmip RH_ECCKD_SOURCE_PATH=/shared/home/greg/.julia/artifacts/7b210aef53e908cfe3c709945f0763c37ca82aaa/ecckd-6115f9b8e29a55cb0f48916857bdc77fec41badd RH_ECCKD_LBL_WORKDIR=/shared/home/greg/ecckd-derived-flux-work bash validation/generate_ecckd_derived_fluxes.sh
RH_ECCKD_DERIVED_FLUX_DRY_RUN=false RH_CKDMIP_DATA_PATH=/shared/home/greg/data/ckdmip RH_ECCKD_SOURCE_PATH=/shared/home/greg/.julia/artifacts/7b210aef53e908cfe3c709945f0763c37ca82aaa/ecckd-6115f9b8e29a55cb0f48916857bdc77fec41badd RH_ECCKD_LBL_WORKDIR=/shared/home/greg/ecckd-derived-flux-work RH_CKDMIP_TOOL_DIR=/path/to/ckdmip/bin bash validation/generate_ecckd_derived_fluxes.sh
# The launcher patches run_lw_lbl_evaluation.sh to SCENARIOS="5gas-180 5gas-280 5gas-415 5gas-560 5gas-1120 5gas-2240 rel-180 rel-280 rel-415 rel-560 rel-1120 rel-2240".
# The launcher patches run_sw_lbl_evaluation.sh to SCENARIOS="rel-180 rel-280 rel-415 rel-560 rel-1120 rel-2240".
RH_CKDMIP_DATA_PATH="$RH_CKDMIP_DATA_PATH" julia --project=test validation/ckdmip_training_data_preflight.jl
```
