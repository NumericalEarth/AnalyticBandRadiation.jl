#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${RH_CKDMIP_DATA_PATH:-}" ]]; then
    echo "Set RH_CKDMIP_DATA_PATH to the CKDMIP data root." >&2
    exit 2
fi

RH_ECCKD_LBL_WORKDIR=${RH_ECCKD_LBL_WORKDIR:-"${HOME:-/tmp}/ecckd-derived-flux-work"}

LW_SOURCE_DIR="${RH_ECCKD_LBL_WORKDIR}/work/lw_lbl_fluxes"
SW_SOURCE_DIR="${RH_ECCKD_LBL_WORKDIR}/work/sw_lbl_fluxes"
LW_TARGET_DIR="${RH_CKDMIP_DATA_PATH}/evaluation1/lw_fluxes"
SW_TARGET_DIR="${RH_CKDMIP_DATA_PATH}/evaluation1/sw_fluxes"

mkdir -p "${LW_TARGET_DIR}" "${SW_TARGET_DIR}"

installed=0
for path in \
    "${LW_SOURCE_DIR}"/ckdmip_evaluation1_lw_fluxes_5gas-*.h5 \
    "${LW_SOURCE_DIR}"/ckdmip_evaluation1_lw_fluxes_rel-*.h5
do
    if [[ -e "${path}" ]]; then
        install -m 0644 "${path}" "${LW_TARGET_DIR}/"
        installed=$((installed + 1))
    fi
done

for path in "${SW_SOURCE_DIR}"/ckdmip_evaluation1_sw_fluxes_rel-*.h5
do
    if [[ -e "${path}" ]]; then
        install -m 0644 "${path}" "${SW_TARGET_DIR}/"
        installed=$((installed + 1))
    fi
done

echo "Installed ${installed} derived ecCKD flux product(s) into ${RH_CKDMIP_DATA_PATH}/evaluation1."
