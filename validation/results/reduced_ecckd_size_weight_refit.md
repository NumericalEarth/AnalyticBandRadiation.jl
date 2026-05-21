# Reduced ecCKD Size Weight Refit

Status: **refit_still_failed**

This diagnostic refits only shortwave weights for larger reduced
official-g-point candidates. It tests whether the simple size-scan
failures are caused by missing g-points or by naive weights.

| Method | ng_lw | ng_sw | Passed | Initial objective | Refit objective | Worst TOA forcing error | Worst surface forcing error |
|---|---:|---:|---:|---:|---:|---:|---:|
| even_select | 32 | 16 | false | 221.891704602 | 221.84624119 | 15.1406161051 W m^-2 | 66.5538723569 W m^-2 |
| even_select | 32 | 20 | false | 494.111026801 | 494.110490438 | 38.6133639731 W m^-2 | 148.233147131 W m^-2 |
| even_select | 32 | 24 | false | 138.230050117 | 138.229269146 | 12.7082198868 W m^-2 | 41.4687807437 W m^-2 |
| even_select | 32 | 28 | false | 85.6938850481 | 85.6930761267 | 18.9034884384 W m^-2 | 25.707922838 W m^-2 |
| even_select | 32 | 30 | false | 44.7715051806 | 44.7702880713 | 6.66257178559 W m^-2 | 13.4310864214 W m^-2 |
| weighted_bins | 32 | 16 | false | 549.944116676 | 549.904826764 | 66.0460093072 W m^-2 | 164.971448029 W m^-2 |
| weighted_bins | 32 | 20 | false | 469.020034901 | 468.961155253 | 48.0854926084 W m^-2 | 140.688346576 W m^-2 |
| weighted_bins | 32 | 24 | false | 161.611035974 | 161.610424139 | 16.8469112018 W m^-2 | 48.4831272418 W m^-2 |
| weighted_bins | 32 | 28 | false | 83.9692902546 | 83.968474503 | 6.55295837346 W m^-2 | 25.1905423509 W m^-2 |
| weighted_bins | 32 | 30 | false | 45.9700100591 | 45.9616905991 | 4.48262977834 W m^-2 | 13.7885071797 W m^-2 |
