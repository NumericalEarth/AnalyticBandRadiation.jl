#!/usr/bin/env bash
set -euo pipefail

# Opt-in helper for materializing the CKDMIP data tree required by the
# original ecCKD training scripts. This intentionally does not run from tests
# or package installation: the full line-by-line database is hundreds of GB to
# about 1 TB, so downloads must be user- or batch-job-owned.

BASE_URL="${CKDMIP_BASE_URL:-https://aux.ecmwf.int/ecpds/home/ckdmip}"
TARGET="${RH_CKDMIP_DATA_PATH:-}"
DRY_RUN="${CKDMIP_DRY_RUN:-true}"
WGET="${WGET:-wget}"
WGET_OPTS="${WGET_OPTS:---continue --tries=20 --timeout=60}"
JULIA="${JULIA:-julia}"
RUN_PREFLIGHT="${CKDMIP_RUN_PREFLIGHT:-false}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [[ -z "${TARGET}" ]]; then
    echo "Set RH_CKDMIP_DATA_PATH to a large target directory before running." >&2
    exit 2
fi

run_cmd() {
    if [[ "${DRY_RUN}" == "true" ]]; then
        printf 'DRY-RUN:'
        printf ' %q' "$@"
        printf '\n'
    else
        "$@"
    fi
}

download_file() {
    local source_path="$1"
    local destination_path="$2"
    local destination_dir
    destination_dir="$(dirname "${TARGET}/${destination_path}")"
    run_cmd mkdir -p "${destination_dir}"
    # shellcheck disable=SC2086
    run_cmd "${WGET}" ${WGET_OPTS} -O "${TARGET}/${destination_path}" "${BASE_URL}/${source_path}"
}

download_directory_flat() {
    local source_path="$1"
    local destination_path="$2"
    run_cmd mkdir -p "${TARGET}/${destination_path}"
    # shellcheck disable=SC2086
    run_cmd "${WGET}" ${WGET_OPTS} -r -np -nd -R 'index.html*' \
        -P "${TARGET}/${destination_path}" "${BASE_URL}/${source_path}/"
}

scenario_concentrations=(
    present preindustrial future glacialmax
    co2-180 co2-280 co2-560 co2-1120 co2-2240
    ch4-350 ch4-700 ch4-1200 ch4-2600 ch4-3500
    n2o-190 n2o-270 n2o-405 n2o-540
    cfc11-0 cfc11-2000 cfc12-0 cfc12-550
    co2-180-ch4-350 co2-180-ch4-3500
    co2-2240-ch4-350 co2-2240-ch4-3500
    co2-180-n2o-190 co2-180-n2o-540
    co2-2240-n2o-190 co2-2240-n2o-540
    ch4-350-n2o-190 ch4-350-n2o-540
    ch4-3500-n2o-190 ch4-3500-n2o-540
)

download_file "concentrations/ckdmip_mmm_concentrations.nc" \
    "mmm/conc/ckdmip_mmm_concentrations.nc"
download_file "concentrations/ckdmip_idealized_concentrations.nc" \
    "idealized/conc/ckdmip_idealized_concentrations.nc"

for dataset in evaluation1 evaluation2; do
    for scenario in "${scenario_concentrations[@]}"; do
        filename="ckdmip_${dataset}_concentrations_${scenario}.nc"
        download_file "concentrations/${filename}" "${dataset}/conc/${filename}"
    done
done

for dataset in mmm idealized evaluation1 evaluation2; do
    download_directory_flat "lw_spectra/${dataset}" "${dataset}/lw_spectra"
    download_directory_flat "sw_spectra/${dataset}" "${dataset}/sw_spectra"
done

for dataset in evaluation1 evaluation2; do
    download_directory_flat "lw_fluxes/${dataset}" "${dataset}/lw_fluxes"
    download_directory_flat "sw_fluxes/${dataset}" "${dataset}/sw_fluxes"
done

run_cmd mkdir -p "${TARGET}/mmm/sw_spectra_extras"
run_cmd ln -sf "../../evaluation1/sw_spectra/ckdmip_ssi.h5" \
    "${TARGET}/mmm/sw_spectra_extras/ckdmip_ssi.h5"

echo "CKDMIP download helper finished. Run:"
echo "  RH_CKDMIP_DATA_PATH=${TARGET} julia --project=test validation/ckdmip_training_data_preflight.jl"

if [[ "${RUN_PREFLIGHT}" == "true" ]]; then
    run_cmd env "RH_CKDMIP_DATA_PATH=${TARGET}" "${JULIA}" \
        --project="${REPO_ROOT}/test" \
        "${REPO_ROOT}/validation/ckdmip_training_data_preflight.jl"
fi
