# Published ecCKD Model Accuracy

Status: **passed**

Reference scope: clean ecCKD cloudless/no-aerosol tropical and RCEMIP-style cases, using matched ecRad reference products for every promoted non-32x32 published combination.

| Model | LW | SW | Passed | Worst TOA forcing | Worst surface forcing | Hard objective | Limiting metric |
|---|---:|---:|---:|---:|---:|---:|---|
| official ecCKD 1.0 32-LW x 32-SW climate model | 32 | 32 | true | 0.00806408713606 W m^-2 | 0.0140335470378 W m^-2 | 0.182186683249 | heating_rate_max_abs |
| official ecCKD 1.0/1.2 32-LW x 64-SW climate/window model | 32 | 64 | true | 0.00805285825516 W m^-2 | 0.014141867794 W m^-2 | 0.1883533103 | heating_rate_max_abs |
| official ecCKD 1.0/1.4 32-LW x 96-SW climate/vfine model | 32 | 96 | true | 0.00801132011759 W m^-2 | 0.0141703247564 W m^-2 | 0.163908639602 | heating_rate_max_abs |
| official ecCKD 1.2/1.4 64-LW x 32-SW narrow/rgb model | 64 | 32 | true | 0.00808122718877 W m^-2 | 0.0137164799025 W m^-2 | 0.184823929466 | heating_rate_max_abs |
| official ecCKD 1.2 64-LW x 64-SW climate model | 64 | 64 | true | 0.00806999830786 W m^-2 | 0.0138248006588 W m^-2 | 0.155860684411 | heating_rate_max_abs |
| official ecCKD 1.2/1.4 64-LW x 96-SW climate/vfine model | 64 | 96 | true | 0.00802846017029 W m^-2 | 0.0138532576211 W m^-2 | 0.173093167827 | heating_rate_max_abs |

## Boundary compatibility

| Model | LW surface spectral matches | SW surface albedo matches | SW direct albedo matches | SW incoming spectral matches |
|---|---:|---:|---:|---:|
| official ecCKD 1.0 32-LW x 32-SW climate model | true | true | false | false |
| official ecCKD 1.0/1.2 32-LW x 64-SW climate/window model | true | true | false | false |
| official ecCKD 1.0/1.4 32-LW x 96-SW climate/vfine model | true | true | false | false |
| official ecCKD 1.2/1.4 64-LW x 32-SW narrow/rgb model | true | true | false | false |
| official ecCKD 1.2 64-LW x 64-SW climate model | true | true | false | false |
| official ecCKD 1.2/1.4 64-LW x 96-SW climate/vfine model | true | true | false | false |

## Mixed-component isolation diagnostics

| Diagnostic | LW | SW | Passed | Worst TOA forcing | Worst surface forcing | Hard objective | Limiting metric |
|---|---:|---:|---:|---:|---:|---:|---|
| 32-LW reference with 64-SW published component | 32 | 64 | false | 1.03319239342 W m^-2 | 3.37183442967 W m^-2 | 11.2394480989 | surface_forcing |
| 32-LW reference with 96-SW published component | 32 | 96 | false | 0.989270364656 W m^-2 | 3.27806404982 W m^-2 | 10.9268801661 | surface_forcing |
| 64-LW published component with 32-SW reference | 64 | 32 | false | 3.56125541041 W m^-2 | 0.752565639935 W m^-2 | 28.102842578 | heating_rate_max_abs |

## Boundary-projection experiment

| Diagnostic | LW | SW | Passed | Worst TOA forcing | Worst surface forcing | Hard objective | Limiting metric |
|---|---:|---:|---:|---:|---:|---:|---|
| 64-LW x 64-SW with 32-g boundary arrays projected by gpoint_fraction overlap | 64 | 64 | false | 114.609965762 W m^-2 | 819.85956528 W m^-2 | 2732.8652176 | surface_forcing |
| 64-LW x 96-SW with 32-g boundary arrays projected by gpoint_fraction overlap | 64 | 96 | false | 114.553404437 W m^-2 | 820.018417133 W m^-2 | 2733.39472378 | surface_forcing |
| 32-LW x 64-SW with 32-g SW boundary arrays projected by gpoint_fraction overlap | 32 | 64 | false | 0.869109925447 W m^-2 | 0.390494040902 W m^-2 | 2.89703308482 | toa_forcing |
| 32-LW x 96-SW with 32-g SW boundary arrays projected by gpoint_fraction overlap | 32 | 96 | false | 0.899594132991 W m^-2 | 0.463188269764 W m^-2 | 2.99864710997 | toa_forcing |
| 64-LW x 32-SW with 32-g LW boundary arrays projected by gpoint_fraction overlap | 64 | 32 | false | 114.567073586 W m^-2 | 820.023239425 W m^-2 | 2733.41079808 | surface_forcing |

This artifact evaluates published full-accuracy ecCKD CKD-definition combinations directly against the same clean package-native ecRad reference cases used by the reduced-model gate.
The mixed-component rows are old-reference diagnostics only: they evaluate promoted definition combinations against the original 32x32 package-native references to show why matched spectral-boundary products are required.
Boundary compatibility records whether the package-native reference files contain spectral boundary arrays with the same g-point count as the tested model; mismatches fall back to broadband boundary treatment in the current evaluator.
The boundary-projection experiment is also diagnostic only: it projects available 32-g boundary arrays onto 64/96 g-grids by official `gpoint_fraction` overlap to test whether missing boundary-grid compatibility is the dominant error source.
