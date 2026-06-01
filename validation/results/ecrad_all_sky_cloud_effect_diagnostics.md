# ecRad All-Sky Cloud-Effect Diagnostics

This diagnostic compares cloud radiative effects, defined as total-sky net flux minus clear-sky net flux, between ecRad references and `radiative_heating_*` candidate variables.

## all_sky_tropical_column

Status: **diagnosed**

Path: `validation/reference/ecrad/all_sky_tropical_column.nc`

### Boundary Cloud Effect Error

| Boundary | Component | RMSE | Max abs | Mean bias | Mean abs | Units |
|---|---|---:|---:|---:|---:|---|
| toa | lw | 3.12722239965 | 9.14527207194 | -0.924335606512 | 1.64890504629 | W m^-2 |
| toa | sw | 10.632772438 | 22.5296647866 | 5.28430078653 | 7.05066177672 | W m^-2 |
| toa | total | 10.2203559211 | 22.1567535461 | 4.35996518002 | 7.19961814986 | W m^-2 |
| surface | lw | 1.65771896388 | 3.40752780308 | -0.420178312636 | 1.10178744054 | W m^-2 |
| surface | sw | 16.5438887109 | 39.6712747378 | 8.45064609069 | 10.7492623945 | W m^-2 |
| surface | total | 16.0912131979 | 39.6483432211 | 8.03046777806 | 10.6197508868 | W m^-2 |

### Profile Cloud Effect Error

| Component | RMSE | Max abs | Mean bias | Mean abs | Units |
|---|---:|---:|---:|---:|---|
| lw | 3.72638511579 | 20.0192400559 | -1.00554769876 | 1.95930439605 | W m^-2 |
| sw | 15.2303137845 | 55.6327097666 | 7.58694752113 | 9.39997217121 | W m^-2 |
| total | 13.6737079217 | 41.0500358212 | 6.58139982237 | 9.02529507873 | W m^-2 |
| heating_rate | 0.621702738881 | 7.77571187081 | -0.0694219693257 | 0.191859988107 | K day^-1 |

## ecckd_all_sky_tropical_column

Status: **diagnosed**

Path: `validation/reference/ecrad/ecckd_all_sky_tropical_column.nc`

### Boundary Cloud Effect Error

| Boundary | Component | RMSE | Max abs | Mean bias | Mean abs | Units |
|---|---|---:|---:|---:|---:|---|
| toa | lw | 0.0433570381216 | 0.11927976108 | 0.0268216077742 | 0.0268216077742 | W m^-2 |
| toa | sw | 0.00411567929992 | 0.00852997698621 | 0.0028310599709 | 0.00283106148203 | W m^-2 |
| toa | total | 0.0459546717159 | 0.12181664447 | 0.0296526677451 | 0.0296526677451 | W m^-2 |
| surface | lw | 0.260672050223 | 0.738309685835 | 0.156024909177 | 0.156024909177 | W m^-2 |
| surface | sw | 0.0031723861552 | 0.00624643419815 | 0.00223758239063 | 0.00223825173012 | W m^-2 |
| surface | total | 0.262591402409 | 0.740687884607 | 0.158262491567 | 0.158262491567 | W m^-2 |

### Profile Cloud Effect Error

| Component | RMSE | Max abs | Mean bias | Mean abs | Units |
|---|---:|---:|---:|---:|---|
| lw | 0.0961445798646 | 0.738309685835 | 0.0450488267791 | 0.0450493536772 | W m^-2 |
| sw | 0.00376641287361 | 0.00875834214355 | 0.00260725284653 | 0.00261232819798 | W m^-2 |
| total | 0.0978679724984 | 0.740687884607 | 0.0476560796257 | 0.0476569665398 | W m^-2 |
| heating_rate | 0.00771176159438 | 0.112223865109 | -0.00104610655551 | 0.00243130570169 | K day^-1 |
