# ecRad Flux Bias Diagnostics

This diagnostic reports signed candidate-minus-reference biases for `radiative_heating_*` variables. It complements the hard absolute-threshold gate by showing error direction.

## clear_sky_tropical_column

Path: `validation/reference/ecrad/clear_sky_tropical_column.nc`

### Boundary Net-Flux Bias

| Boundary | Mean bias | Min bias | Max bias | Mean abs | Units |
|---|---:|---:|---:|---:|---|
| toa | -0.0383487909363 | -1.8531127012 | 0.989114527618 | 0.811571246791 | W m^-2 |
| surface | -2.83274233097 | -4.01937447101 | -1.5897072919 | 2.83274233097 | W m^-2 |

### Variable Bias

| Variable | Mean bias | Min bias | Max bias | Mean abs | Units |
|---|---:|---:|---:|---:|---|
| `lw_up` | 1.46053301274 | -0.376889069477 | 6.04610532155 | 1.47893920471 | W m^-2 |
| `lw_down` | -0.0486199918057 | -5.74028576438 | 1.77006850157 | 0.37568354281 | W m^-2 |
| `sw_up` | -0.821674182577 | -1.87463611319 | 2.47627486961 | 1.04092677452 | W m^-2 |
| `sw_down` | 1.21616326564 | -3.94200124601 | 4.14992189704 | 1.91926584206 | W m^-2 |
| `heating_rate` | -0.359906391578 | -13.1696360171 | 9.47104824036 | 0.500386958453 | K day^-1 |

## all_sky_tropical_column

Path: `validation/reference/ecrad/all_sky_tropical_column.nc`

### Boundary Net-Flux Bias

| Boundary | Mean bias | Min bias | Max bias | Mean abs | Units |
|---|---:|---:|---:|---:|---|
| toa | 92.9218818733 | 20.3459945778 | 217.429653213 | 92.9218818733 | W m^-2 |
| surface | 31.3659925333 | -164.624931579 | 208.007750236 | 84.1149336946 | W m^-2 |

### Variable Bias

| Variable | Mean bias | Min bias | Max bias | Mean abs | Units |
|---|---:|---:|---:|---:|---|
| `lw_up` | -5.03896916792 | -48.8410838611 | 14.9481269496 | 7.78483721682 | W m^-2 |
| `lw_down` | 3.01524145333 | -8.73128530385 | 59.4812161281 | 3.70157506875 | W m^-2 |
| `sw_up` | -66.87996295 | -260.454726641 | 17.2111461205 | 68.2169149419 | W m^-2 |
| `sw_down` | -6.42930563959 | -425.797467574 | 212.667935854 | 37.4952536066 | W m^-2 |
| `heating_rate` | -0.120352838804 | -14.3229257322 | 45.2851087667 | 1.28412128527 | K day^-1 |

## ecckd_clear_sky_tropical_column

Path: `validation/reference/ecrad/ecckd_clear_sky_tropical_column.nc`

### Boundary Net-Flux Bias

| Boundary | Mean bias | Min bias | Max bias | Mean abs | Units |
|---|---:|---:|---:|---:|---|
| toa | -0.028149930378 | -0.0313832809406 | -0.0261437526692 | 0.028149930378 | W m^-2 |
| surface | -0.0099510373663 | -0.011339478423 | -0.00506376069495 | 0.0099510373663 | W m^-2 |

### Variable Bias

| Variable | Mean bias | Min bias | Max bias | Mean abs | Units |
|---|---:|---:|---:|---:|---|
| `lw_up` | 0.000200227925928 | -0.00178466796785 | 0.00248171115015 | 0.000742844987271 | W m^-2 |
| `lw_down` | -0.000471071725031 | -0.00610034046844 | 0.00225655876281 | 0.00127281774608 | W m^-2 |
| `sw_up` | 0.0194628957682 | -0.000768044720189 | 0.0317722276526 | 0.0194903053396 | W m^-2 |
| `sw_down` | -0.001637288035 | -0.00576806610809 | 9.8667775319e-05 | 0.00164083408944 | W m^-2 |
| `heating_rate` | -6.03899578596e-05 | -0.0540071082913 | 0.0911075171087 | 0.0013148115788 | K day^-1 |

## rcemip_style_column_subset

Path: `validation/reference/ecrad/rcemip_style_column_subset.nc`

### Boundary Net-Flux Bias

| Boundary | Mean bias | Min bias | Max bias | Mean abs | Units |
|---|---:|---:|---:|---:|---|
| toa | -0.122253463912 | -1.8531127012 | 0.989114527618 | 0.584519136816 | W m^-2 |
| surface | -0.496085853216 | -4.01937447101 | 2.48348938873 | 1.58758003138 | W m^-2 |

### Variable Bias

| Variable | Mean bias | Min bias | Max bias | Mean abs | Units |
|---|---:|---:|---:|---:|---|
| `lw_up` | 1.20003573522 | -1.15423000678 | 6.04610532155 | 1.21026365821 | W m^-2 |
| `lw_down` | 0.221301116135 | -5.74028576438 | 2.48348938873 | 0.426150882248 | W m^-2 |
| `sw_up` | -0.722287873748 | -2.93105723734 | 2.47627486961 | 0.818815992451 | W m^-2 |
| `sw_down` | 0.736240392647 | -3.94200124601 | 4.14992189704 | 1.30076219668 | W m^-2 |
| `heating_rate` | -0.264377945457 | -13.1696360171 | 9.47104824036 | 0.356282426409 | K day^-1 |

## ecckd_rcemip_style_column_subset

Path: `validation/reference/ecrad/ecckd_rcemip_style_column_subset.nc`

### Boundary Net-Flux Bias

| Boundary | Mean bias | Min bias | Max bias | Mean abs | Units |
|---|---:|---:|---:|---:|---|
| toa | -0.0266108921988 | -0.0719878208337 | 0.000549983924287 | 0.026698408905 | W m^-2 |
| surface | -0.0045266119556 | -0.011339478423 | 0.00207154079555 | 0.00540458759134 | W m^-2 |

### Variable Bias

| Variable | Mean bias | Min bias | Max bias | Mean abs | Units |
|---|---:|---:|---:|---:|---|
| `lw_up` | -0.000322139199521 | -0.00219969629183 | 0.00248171115015 | 0.000654038735735 | W m^-2 |
| `lw_down` | 0.000316707660218 | -0.00610034046844 | 0.00228364640671 | 0.000906182557018 | W m^-2 |
| `sw_up` | 0.0187162217699 | -0.000832297792755 | 0.0742817861557 | 0.0187430223354 | W m^-2 |
| `sw_down` | -0.000819463325102 | -0.00678073921642 | 0.00766861910938 | 0.00138156046158 | W m^-2 |
| `heating_rate` | -0.000130007244055 | -0.0577565077396 | 0.0911075171087 | 0.0010228924372 | K day^-1 |
