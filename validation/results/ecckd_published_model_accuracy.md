# Published ecCKD Model Accuracy

Status: **failed_threshold**

Reference scope: clean ecCKD cloudless/no-aerosol tropical and RCEMIP-style cases.

| Model | LW | SW | Passed | Worst TOA forcing | Worst surface forcing | Hard objective | Limiting metric |
|---|---:|---:|---:|---:|---:|---:|---|
| official ecCKD 1.0 32-LW x 32-SW climate model | 32 | 32 | true | 0.00806408713606 W m^-2 | 0.0140335470378 W m^-2 | 0.182186683249 | heating_rate_max_abs |
| official ecCKD 1.2 64-LW x 64-SW climate model | 64 | 64 | false | 2.96502926216 W m^-2 | 3.73485800638 W m^-2 | 28.6371129464 | heating_rate_max_abs |
| official ecCKD 1.2/1.4 64-LW x 96-SW climate/vfine model | 64 | 96 | false | 2.77131137562 W m^-2 | 3.59621389676 W m^-2 | 28.6607206965 | heating_rate_max_abs |

## Mixed-component isolation diagnostics

| Diagnostic | LW | SW | Passed | Worst TOA forcing | Worst surface forcing | Hard objective | Limiting metric |
|---|---:|---:|---:|---:|---:|---:|---|
| 32-LW reference with 64-SW published component | 32 | 64 | false | 1.03319239342 W m^-2 | 3.37183442967 W m^-2 | 11.2394480989 | surface_forcing |
| 32-LW reference with 96-SW published component | 32 | 96 | false | 0.989270364656 W m^-2 | 3.27806404982 W m^-2 | 10.9268801661 | surface_forcing |
| 64-LW published component with 32-SW reference | 64 | 32 | false | 3.56125541041 W m^-2 | 0.752565639935 W m^-2 | 28.102842578 | heating_rate_max_abs |

This artifact evaluates published full-accuracy ecCKD CKD-definition combinations directly against the same clean package-native ecRad reference cases used by the reduced-model gate.
The mixed-component rows are diagnostic only: they are not published ecCKD products, but they localize the current 64/96 parity failure by swapping one published component at a time against the passing 32x32 baseline.
