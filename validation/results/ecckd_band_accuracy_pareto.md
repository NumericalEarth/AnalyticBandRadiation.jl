# ecCKD Band-Count Accuracy Pareto

Status: **passed**

- Accuracy points: 22
- Passing accuracy points: 1
- Published ecCKD inventory entries: 6
- Plot: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_band_accuracy_pareto.svg`
- CSV: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_band_accuracy_pareto.csv`

This artifact plots all currently available ecCKD reduced-accuracy rows. Published 64/96 definitions are inventoried and recovered by the teacher-student scan, but radiative accuracy for newly trained intermediate models remains future work.

## Pareto Front

| Total g-points | LW | SW | Passed | Worst boundary forcing error | Method |
|---:|---:|---:|---:|---:|---|
| 32 | 16 | 16 | false | 156.746 | adjacent official ecCKD g-point bins with spectral-weighted coefficient averages |
| 48 | 32 | 16 | false | 2.07524 | weighted greedy 16 shortwave g-point subset with boundary-aware table, component, structural, objective-probe, surface-probe, capped table, continuation, post-capped weight, post-weight surface-table, bounded weight, four current component-scale refits, selected current gas-pressure component scan refit, gas-pressure continuation refit, weighted gas-pressure continuation refit, and high-weight gas-pressure continuation refit |
| 64 | 32 | 32 | true | 0.0140335 | official ecCKD 32x32 baseline without shortwave reduction |

## Published Inventory

| File | Kind | Bands | G-points |
|---|---|---:|---:|
| `ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc` | longwave | 1 | 32 |
| `ecckd-1.0_sw_climate_rgb-32b_ckd-definition.nc` | shortwave | 5 | 32 |
| `ecckd-1.2_lw_climate_narrow-64b_ckd-definition.nc` | longwave | 13 | 64 |
| `ecckd-1.2_sw_climate_window-64b_ckd-definition.nc` | shortwave | 19 | 64 |
| `ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc` | shortwave | 5 | 32 |
| `ecckd-1.4_sw_climate_vfine-96b_ckd-definition.nc` | shortwave | 44 | 96 |
