# CKDMIP Training Data Download Plan

Status: **manual_or_batch_download_required**

This plan intentionally does not run downloads. The full CKDMIP LBL database is hundreds of GB to about 1 TB and should be materialized by an explicit user or batch job.

Set `RH_CKDMIP_DATA_PATH` to a large scratch/data directory before running any command.

Base URL: `https://aux.ecmwf.int/ecpds/home/ckdmip`

Tasks: 83

## Expected Layout

- `mmm/conc`
- `mmm/lw_spectra`
- `mmm/sw_spectra`
- `mmm/sw_spectra_extras`
- `idealized/conc`
- `idealized/lw_spectra`
- `idealized/sw_spectra`
- `evaluation1/conc`
- `evaluation1/lw_spectra`
- `evaluation1/sw_spectra`
- `evaluation1/lw_fluxes`
- `evaluation1/sw_fluxes`
- `evaluation2/conc`
- `evaluation2/lw_spectra`
- `evaluation2/sw_spectra`
- `evaluation2/lw_fluxes`
- `evaluation2/sw_fluxes`

## Commands

```bash
test -n "$RH_CKDMIP_DATA_PATH"
mkdir -p "$RH_CKDMIP_DATA_PATH/mmm/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/mmm/conc/ckdmip_mmm_concentrations.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_mmm_concentrations.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/idealized/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/idealized/conc/ckdmip_idealized_concentrations.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_idealized_concentrations.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_present.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_present.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_preindustrial.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_preindustrial.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_future.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_future.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_glacialmax.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_glacialmax.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_co2-180.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_co2-180.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_co2-280.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_co2-280.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_co2-560.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_co2-560.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_co2-1120.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_co2-1120.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_co2-2240.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_co2-2240.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_ch4-350.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_ch4-350.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_ch4-700.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_ch4-700.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_ch4-1200.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_ch4-1200.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_ch4-2600.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_ch4-2600.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_ch4-3500.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_ch4-3500.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_n2o-190.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_n2o-190.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_n2o-270.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_n2o-270.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_n2o-405.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_n2o-405.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_n2o-540.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_n2o-540.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_cfc11-0.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_cfc11-0.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_cfc11-2000.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_cfc11-2000.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_cfc12-0.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_cfc12-0.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_cfc12-550.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_cfc12-550.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_co2-180-ch4-350.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_co2-180-ch4-350.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_co2-180-ch4-3500.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_co2-180-ch4-3500.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_co2-2240-ch4-350.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_co2-2240-ch4-350.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_co2-2240-ch4-3500.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_co2-2240-ch4-3500.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_co2-180-n2o-190.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_co2-180-n2o-190.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_co2-180-n2o-540.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_co2-180-n2o-540.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_co2-2240-n2o-190.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_co2-2240-n2o-190.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_co2-2240-n2o-540.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_co2-2240-n2o-540.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_ch4-350-n2o-190.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_ch4-350-n2o-190.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_ch4-350-n2o-540.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_ch4-350-n2o-540.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_ch4-3500-n2o-190.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_ch4-3500-n2o-190.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation1/conc/ckdmip_evaluation1_concentrations_ch4-3500-n2o-540.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation1_concentrations_ch4-3500-n2o-540.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_present.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_present.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_preindustrial.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_preindustrial.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_future.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_future.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_glacialmax.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_glacialmax.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_co2-180.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_co2-180.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_co2-280.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_co2-280.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_co2-560.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_co2-560.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_co2-1120.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_co2-1120.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_co2-2240.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_co2-2240.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_ch4-350.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_ch4-350.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_ch4-700.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_ch4-700.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_ch4-1200.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_ch4-1200.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_ch4-2600.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_ch4-2600.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_ch4-3500.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_ch4-3500.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_n2o-190.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_n2o-190.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_n2o-270.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_n2o-270.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_n2o-405.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_n2o-405.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_n2o-540.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_n2o-540.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_cfc11-0.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_cfc11-0.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_cfc11-2000.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_cfc11-2000.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_cfc12-0.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_cfc12-0.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_cfc12-550.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_cfc12-550.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_co2-180-ch4-350.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_co2-180-ch4-350.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_co2-180-ch4-3500.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_co2-180-ch4-3500.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_co2-2240-ch4-350.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_co2-2240-ch4-350.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_co2-2240-ch4-3500.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_co2-2240-ch4-3500.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_co2-180-n2o-190.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_co2-180-n2o-190.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_co2-180-n2o-540.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_co2-180-n2o-540.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_co2-2240-n2o-190.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_co2-2240-n2o-190.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_co2-2240-n2o-540.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_co2-2240-n2o-540.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_ch4-350-n2o-190.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_ch4-350-n2o-190.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_ch4-350-n2o-540.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_ch4-350-n2o-540.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_ch4-3500-n2o-190.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_ch4-3500-n2o-190.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/conc" && wget -c -O "$RH_CKDMIP_DATA_PATH/evaluation2/conc/ckdmip_evaluation2_concentrations_ch4-3500-n2o-540.nc" "https://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/ckdmip_evaluation2_concentrations_ch4-3500-n2o-540.nc"
mkdir -p "$RH_CKDMIP_DATA_PATH/mmm/lw_spectra" && wget -r -np -nd -R 'index.html*' -P "$RH_CKDMIP_DATA_PATH/mmm/lw_spectra" "https://aux.ecmwf.int/ecpds/home/ckdmip/lw_spectra/mmm/"
mkdir -p "$RH_CKDMIP_DATA_PATH/mmm/sw_spectra" && wget -r -np -nd -R 'index.html*' -P "$RH_CKDMIP_DATA_PATH/mmm/sw_spectra" "https://aux.ecmwf.int/ecpds/home/ckdmip/sw_spectra/mmm/"
mkdir -p "$RH_CKDMIP_DATA_PATH/idealized/lw_spectra" && wget -r -np -nd -R 'index.html*' -P "$RH_CKDMIP_DATA_PATH/idealized/lw_spectra" "https://aux.ecmwf.int/ecpds/home/ckdmip/lw_spectra/idealized/"
mkdir -p "$RH_CKDMIP_DATA_PATH/idealized/sw_spectra" && wget -r -np -nd -R 'index.html*' -P "$RH_CKDMIP_DATA_PATH/idealized/sw_spectra" "https://aux.ecmwf.int/ecpds/home/ckdmip/sw_spectra/idealized/"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/lw_spectra" && wget -r -np -nd -R 'index.html*' -P "$RH_CKDMIP_DATA_PATH/evaluation1/lw_spectra" "https://aux.ecmwf.int/ecpds/home/ckdmip/lw_spectra/evaluation1/"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/sw_spectra" && wget -r -np -nd -R 'index.html*' -P "$RH_CKDMIP_DATA_PATH/evaluation1/sw_spectra" "https://aux.ecmwf.int/ecpds/home/ckdmip/sw_spectra/evaluation1/"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/lw_spectra" && wget -r -np -nd -R 'index.html*' -P "$RH_CKDMIP_DATA_PATH/evaluation2/lw_spectra" "https://aux.ecmwf.int/ecpds/home/ckdmip/lw_spectra/evaluation2/"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/sw_spectra" && wget -r -np -nd -R 'index.html*' -P "$RH_CKDMIP_DATA_PATH/evaluation2/sw_spectra" "https://aux.ecmwf.int/ecpds/home/ckdmip/sw_spectra/evaluation2/"
mkdir -p "$RH_CKDMIP_DATA_PATH/mmm/sw_spectra_extras" && ln -sf "../../evaluation1/sw_spectra/ckdmip_ssi.h5" "$RH_CKDMIP_DATA_PATH/mmm/sw_spectra_extras/ckdmip_ssi.h5"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/lw_fluxes" && wget -r -np -nd -R 'index.html*' -P "$RH_CKDMIP_DATA_PATH/evaluation1/lw_fluxes" "https://aux.ecmwf.int/ecpds/home/ckdmip/lw_fluxes/evaluation1/"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation1/sw_fluxes" && wget -r -np -nd -R 'index.html*' -P "$RH_CKDMIP_DATA_PATH/evaluation1/sw_fluxes" "https://aux.ecmwf.int/ecpds/home/ckdmip/sw_fluxes/evaluation1/"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/lw_fluxes" && wget -r -np -nd -R 'index.html*' -P "$RH_CKDMIP_DATA_PATH/evaluation2/lw_fluxes" "https://aux.ecmwf.int/ecpds/home/ckdmip/lw_fluxes/evaluation2/"
mkdir -p "$RH_CKDMIP_DATA_PATH/evaluation2/sw_fluxes" && wget -r -np -nd -R 'index.html*' -P "$RH_CKDMIP_DATA_PATH/evaluation2/sw_fluxes" "https://aux.ecmwf.int/ecpds/home/ckdmip/sw_fluxes/evaluation2/"
```
