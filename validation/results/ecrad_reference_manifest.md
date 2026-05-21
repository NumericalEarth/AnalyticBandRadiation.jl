# ecRad Reference Manifest

Status: **references_present_schema_valid**

This report defines the first required external ecRad/ecCKD reference artifacts. It does not run ecRad and does not validate accuracy by itself.

| Case | Status | Missing variables | Path |
|---|---|---:|---|
| clear_sky_tropical_column | present_schema_valid | 0 | `validation/reference/ecrad/clear_sky_tropical_column.nc` |
| all_sky_tropical_column | present_schema_valid | 0 | `validation/reference/ecrad/all_sky_tropical_column.nc` |
| ecckd_clear_sky_tropical_column | present_schema_valid | 0 | `validation/reference/ecrad/ecckd_clear_sky_tropical_column.nc` |
| ecckd_all_sky_tropical_column | present_schema_valid | 0 | `validation/reference/ecrad/ecckd_all_sky_tropical_column.nc` |
| rcemip_style_column_subset | present_schema_valid | 0 | `validation/reference/ecrad/rcemip_style_column_subset.nc` |
| ecckd_rcemip_style_column_subset | present_schema_valid | 0 | `validation/reference/ecrad/ecckd_rcemip_style_column_subset.nc` |

## Initial Hard Thresholds

| Metric | Threshold |
|---|---:|
| Flux RMSE | 1 W m^-2 |
| Flux max abs | 5 W m^-2 |
| Heating-rate RMSE | 0.05 K day^-1 |
| Heating-rate max abs | 0.5 K day^-1 |
| TOA forcing abs error | 0.3 W m^-2 |
| Surface forcing abs error | 0.3 W m^-2 |

Reference file instructions: `validation/reference/ecrad/README.md`

Missing reference files or invalid schemas block final ecRad/ecCKD validation and reduced-model acceptance.
