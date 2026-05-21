#!/bin/bash
#SBATCH --job-name=rh-kmodel-coverage
#SBATCH --partition=gpu-prod
#SBATCH --gres=gpu:h100:1
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=/shared/home/greg/data/h100-kmodel-coverage-%j.log

# H100 RCEMIP-style benchmark across published ecCKD k-models, driven by the
# independent benchmarking project in this repo (which path-deps on the
# developing Breeze checkout). Five k-models at 512x512x128 / 262144 columns.

set -euo pipefail
export PATH="/usr/local/cuda-13.0/bin:${PATH}"

# Absolute path — Slurm copies this script to /var/spool/slurmd before
# executing, so $BASH_SOURCE[0] does not resolve to the repo.
repo_root=/shared/home/greg/Projects/AnalyticBandRadiation.jl
ecrad_data=/shared/home/greg/.julia/artifacts/49ce668ce0861f9d5e8299d68af7138485eb5f19/ecrad-131ac980517719b7a859e3ccc117919a1d888a20/data

cd "${repo_root}"

echo "=== H100 environment ==="
nvidia-smi || true
echo

run_k_model() {
    local label="$1"
    local lw_file="$2"
    local sw_file="$3"
    local result_dir="${repo_root}/benchmarking/results/rcemip_h100_${label}"

    echo
    echo "================================================================"
    echo "k-model run: ${label}   $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "  LW: ${lw_file}"
    echo "  SW: ${sw_file}"
    echo "  output: ${result_dir}"
    echo "================================================================"

    mkdir -p "${result_dir}"

    # `|| true` so the script keeps going if one configuration OOMs.
    RADIATIVE_HEATING_BACKEND=H100 \
    RADIATIVE_HEATING_GAS_MODEL_SOURCE=validated_ecCKD \
    RADIATIVE_HEATING_NX=512 \
    RADIATIVE_HEATING_NY=512 \
    RADIATIVE_HEATING_NZ=128 \
    RADIATIVE_HEATING_SAMPLES=5 \
    RADIATIVE_HEATING_RUN_RRTMGP=true \
    RADIATIVE_HEATING_OUTPUT_DIR="${result_dir}" \
    RADIATIVE_HEATING_ECCKD_LW_PATH="${ecrad_data}/${lw_file}" \
    RADIATIVE_HEATING_ECCKD_SW_PATH="${ecrad_data}/${sw_file}" \
        julia --project=benchmarking benchmarking/rcemip_h100_benchmark.jl \
        || echo "*** k-model ${label} FAILED — see log above; continuing ***"

    echo
}

# Five published-ecCKD k-model combinations. All use the same RCEMIP-style
# 512x512x128 / 262144 column workload, validated ecCKD provenance, and the
# Breeze update_radiation! entry point.
run_k_model "32x32_validated_ecckd" \
            "ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc" \
            "ecckd-1.0_sw_climate_rgb-32b_ckd-definition.nc"

run_k_model "32x64_validated_ecckd" \
            "ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc" \
            "ecckd-1.2_sw_climate_window-64b_ckd-definition.nc"

run_k_model "32x96_validated_ecckd" \
            "ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc" \
            "ecckd-1.4_sw_climate_vfine-96b_ckd-definition.nc"

run_k_model "64x32_validated_ecckd" \
            "ecckd-1.2_lw_climate_narrow-64b_ckd-definition.nc" \
            "ecckd-1.0_sw_climate_rgb-32b_ckd-definition.nc"

run_k_model "64x64_validated_ecckd" \
            "ecckd-1.2_lw_climate_narrow-64b_ckd-definition.nc" \
            "ecckd-1.2_sw_climate_window-64b_ckd-definition.nc"

echo
echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) — all k-model coverage runs done"
