# Reduced ecCKD Post-Table Weight Refit

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Iterations | 2000 |
| Base objective | 8.40968478599 |
| Refit objective | 8.41012609185 |
| Objective reduction | -0.000441305855171 |
| Accepted | false |
| Approximate hard-gate objective | 9.60413137131 |
| Base worst TOA forcing | 2.5229054358 W m^-2 |
| Base worst surface forcing | 2.45622471063 W m^-2 |
| Refit worst TOA forcing | 2.52303782755 W m^-2 |
| Refit worst surface forcing | 2.4540439232 W m^-2 |

This diagnostic refits only the 16 shortwave weights after applying the
current coefficient scales, pressure-band table moves, and best available
targeted active table-entry moves. The refit is accepted only if the exact
hard-gate objective improves.
