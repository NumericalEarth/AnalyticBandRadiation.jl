#!/usr/bin/env bash
set -euo pipefail

LW_SCENARIOS="5gas-180 5gas-280 5gas-415 5gas-560 5gas-1120 5gas-2240 rel-180 rel-280 rel-415 rel-560 rel-1120 rel-2240"
SW_SCENARIOS="rel-180 rel-280 rel-415 rel-560 rel-1120 rel-2240"

RH_ECCKD_DERIVED_FLUX_DRY_RUN=${RH_ECCKD_DERIVED_FLUX_DRY_RUN:-true}
RH_ECCKD_LBL_WORKDIR=${RH_ECCKD_LBL_WORKDIR:-"${HOME:-/tmp}/ecckd-derived-flux-work"}
PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

if [[ -z "${RH_CKDMIP_DATA_PATH:-}" ]]; then
    echo "Set RH_CKDMIP_DATA_PATH to the CKDMIP data root." >&2
    exit 2
fi

if [[ -z "${RH_ECCKD_SOURCE_PATH:-}" ]]; then
    RH_ECCKD_SOURCE_PATH=$(julia --project="${PROJECT_ROOT}" -e 'using Lightflux; print(Lightflux.ecckd_source_path(require=true))')
fi

if [[ "${RH_ECCKD_DERIVED_FLUX_DRY_RUN}" != "true" && -z "${RH_CKDMIP_TOOL_DIR:-}" ]]; then
    echo "Set RH_CKDMIP_TOOL_DIR to the CKDMIP install root or bin directory containing ckdmip_lw and ckdmip_sw, or leave RH_ECCKD_DERIVED_FLUX_DRY_RUN=true." >&2
    exit 2
fi

if [[ "${RH_ECCKD_DERIVED_FLUX_DRY_RUN}" != "true" ]]; then
    if [[ -x "${RH_CKDMIP_TOOL_DIR}/bin/ckdmip_lw" && -x "${RH_CKDMIP_TOOL_DIR}/bin/ckdmip_sw" ]]; then
        RH_CKDMIP_ROOT="${RH_CKDMIP_TOOL_DIR}"
    elif [[ -x "${RH_CKDMIP_TOOL_DIR}/ckdmip_lw" && -x "${RH_CKDMIP_TOOL_DIR}/ckdmip_sw" ]]; then
        RH_CKDMIP_ROOT=$(cd "${RH_CKDMIP_TOOL_DIR}/.." && pwd)
    else
        echo "RH_CKDMIP_TOOL_DIR must contain bin/ckdmip_lw and bin/ckdmip_sw, or be that bin directory." >&2
        exit 2
    fi
else
    RH_CKDMIP_ROOT="${RH_CKDMIP_TOOL_DIR:-/missing/ckdmip-tools}"
fi

mkdir -p "${RH_ECCKD_LBL_WORKDIR}"
rsync -a "${RH_ECCKD_SOURCE_PATH}/" "${RH_ECCKD_LBL_WORKDIR}/ecckd/"

cd "${RH_ECCKD_LBL_WORKDIR}/ecckd/test"

sed 's/@PACKAGE_VERSION@/1.2/g' version.h.in > version.h

sed -i \
    -e "s|^CKDMIP_DIR=.*|CKDMIP_DIR=${RH_CKDMIP_ROOT}|" \
    -e "s|^CKDMIP_DATA_DIR=.*|CKDMIP_DATA_DIR=${RH_CKDMIP_DATA_PATH}|" \
    -e "s|^WORK_DIR=.*|WORK_DIR=${RH_ECCKD_LBL_WORKDIR}/work|" \
    config.h

sed -i -E "0,/^SCENARIOS=\"[^\"]*\"/s//SCENARIOS=\"${LW_SCENARIOS}\"/" run_lw_lbl_evaluation.sh
sed -i -E "0,/^SCENARIOS=\"[^\"]*\"/s//SCENARIOS=\"${SW_SCENARIOS}\"/" run_sw_lbl_evaluation.sh
sed -i '/^module load nco$/s/^/# /' run_lw_lbl_evaluation.sh run_sw_lbl_evaluation.sh
sed -i 's/iverbose = 3/iverbose = 1/g' run_lw_lbl_evaluation.sh run_sw_lbl_evaluation.sh
sed -i -E \
    's|ncrcat -O \$OUTFILES \$OUTFILE|julia --project='"${PROJECT_ROOT}"'/test '"${PROJECT_ROOT}"'/validation/concat_ckdmip_flux_chunks.jl $OUTFILES $OUTFILE|' \
    run_lw_lbl_evaluation.sh run_sw_lbl_evaluation.sh
sed -i '/OUTFILES="$OUTFILES $OUTFILE"/a\
\
\tif [ -s "$OUTFILE" ]\
\tthen\
\t    echo "*** REUSING $OUTFILE ***"\
\t    continue\
\tfi' run_lw_lbl_evaluation.sh run_sw_lbl_evaluation.sh

echo "Prepared ecCKD working copy at ${RH_ECCKD_LBL_WORKDIR}/ecckd"
echo "LW scenarios: ${LW_SCENARIOS}"
echo "SW scenarios: ${SW_SCENARIOS}"

if [[ "${RH_ECCKD_DERIVED_FLUX_DRY_RUN}" == "true" ]]; then
    cat <<EOF
Dry run only. To generate the derived flux products, install/build CKDMIP tools,
then rerun with:

  RH_ECCKD_DERIVED_FLUX_DRY_RUN=false \\
  RH_CKDMIP_DATA_PATH=${RH_CKDMIP_DATA_PATH} \\
  RH_CKDMIP_TOOL_DIR=/path/to/ckdmip/bin \\
  bash validation/generate_ecckd_derived_fluxes.sh

The real run will execute:
  bash run_lw_lbl_evaluation.sh
  bash run_sw_lbl_evaluation.sh
  install generated ckdmip_evaluation1_*_fluxes_{5gas,rel}-*.h5 files into ${RH_CKDMIP_DATA_PATH}/evaluation1/{lw,sw}_fluxes
EOF
    exit 0
fi

bash run_lw_lbl_evaluation.sh
bash run_sw_lbl_evaluation.sh

bash "${PROJECT_ROOT}/validation/install_ecckd_derived_fluxes.sh"
