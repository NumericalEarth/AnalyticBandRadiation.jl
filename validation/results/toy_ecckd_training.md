# Toy ecCKD Training

Status: **pass**

This is a deterministic reduced gas-optics training fixture. It exercises
loss construction, finite-difference gradient consistency, and loss
reduction for fixed-topology ecCKD-style parameters. It is not a CKDMIP
production training run. When Enzyme and Reactant are available in the
active environment, it also checks an Enzyme gradient through a pure
two-stream radiative loss and Reactant compilation of a toy loss. The full
mutating radiative-transfer training path is still open.

| Metric | Value |
|---|---:|
| Initial loss | 0.0113011861681 |
| Final loss | 6.62185338745e-05 |
| Loss ratio | 0.00585943217724 |
| Initial flux RMSE | 0.106306902776 W m^-2 |
| Final flux RMSE | 0.00813742167428 W m^-2 |
| Flux RMSE ratio | 0.0765465032068 |
| Initial heating RMSE | 0.00169086417626 K day^-1 |
| Final heating RMSE | 0.000300394647974 K day^-1 |
| Heating RMSE ratio | 0.177657467816 |
| Gradient relative difference | 6.24261421237e-10 |
| Gradient threshold | 0.0001 |
| Enzyme status | passed |
| Enzyme check | pure_two_stream_radiative_loss |
| Enzyme relative error | 4.27732803781e-10 |
| Reactant status | passed |

Initial parameters: `[0.07, 0.003, 0.01, 0.0004]`

Final parameters: `[0.13969503976166664, 0.009484740972396512, 0.017783592097912335, 0.0015998708590926256]`

Target parameters: `[0.14, 0.006, 0.02, 0.0009]`

Configuration: `(deterministic_seed = "none_deterministic_fixture", iterations = 12, gradient_step = 1.0e-5, gradient_check_step = 5.0e-6, initial_line_search_step = 0.5, minimum_line_search_step = 1.0e-8, loss_flux_weight = 1.0, loss_heating_rate_weight = 0.01, atmosphere = (pressure_layers_pa = [13000.0, 42500.0, 80000.0], pressure_interfaces_pa = [1000.0, 25000.0, 60000.0, 100000.0], temperature_layers_k = [230.0, 262.0, 290.0], h2o = [0.0015, 0.006, 0.018], co2 = 0.00042, surface_temperature_k = 300.0, surface_albedo = 0.08, cos_zenith = 0.5))`
