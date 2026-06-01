# Reduced ecCKD Support Expansion Refit Scan

Status: **failed_threshold**

| Field | Value |
|---|---:|
| Search objective | support_expansion_multistart_hardgate_weight_refit |
| Iterations | 4000 |
| p-norm | 16 |
| Candidates | 5 |
| Starts per candidate | 4 |
| Best label | `subset_hardgate_plus_3_6` |
| Best objective | 50.3454530164 |
| Best start | `official_normalized` |
| Best g-points | 18 |
| Best passed | false |
| Best added g-points | `3, 6` |
| Best indices | `2, 3, 4, 6, 7, 10, 11, 12, 14, 16, 18, 21, 22, 27, 28, 30, 31, 32` |

This diagnostic reruns the best 17/18-g expansion candidates with
longer hard-gate weight optimization and multiple warm starts, including
zero-padded and epsilon-padded seed weights. It tests whether the
bounded support-expansion rejection was caused by the initial weight
refit rather than by the support itself.

| Candidate | Added | ng | Seed objective | Best start | Best objective | TOA forcing | Surface forcing | Heating RMSE | Passed |
|---|---|---:|---:|---|---:|---:|---:|---:|---:|
| `subset_hardgate_plus_3_6` | `3, 6` | 18 | 34.7591363002 | `official_normalized` | 50.3454530164 | 7.48945833657 | 9.90790435354 | 2.51727265082 | false |
| `subset_hardgate_plus_5_6` | `5, 6` | 18 | 34.7591363002 | `official_normalized` | 50.5675822715 | 7.37858185002 | 9.18356632065 | 2.52837911358 | false |
| `canonical_plus_7_11` | `7, 11` | 18 | 142.876457579 | `official_normalized` | 60.6963087702 | 5.51386117679 | 8.66775276885 | 3.03481543851 | false |
| `subset_hardgate_plus_1` | `1` | 17 | 34.7591363002 | `official_normalized` | 60.7530603356 | 4.90953518784 | 8.41748364323 | 3.03765301678 | false |
| `canonical_plus_7` | `7` | 17 | 142.876457579 | `official_normalized` | 72.8234511377 | 6.96582209781 | 15.2870501281 | 3.64117255688 | false |

## subset_hardgate_plus_3_6

| Start | Approximate objective | Exact objective | TOA forcing | Surface forcing | Heating RMSE | Passed |
|---|---:|---:|---:|---:|---:|---:|
| `official_normalized` | 33.8828400384 | 50.3454530164 | 7.48945833657 | 9.90790435354 | 2.51727265082 | false |
| `uniform` | 1841.23202994 | 1743.9014877 | 280.477482336 | 488.596778142 | 87.1950743851 | false |
| `seed_zero_padded` | 34.1309678313 | 67.5356720642 | 4.42624738196 | 6.83464314061 | 3.37678360321 | false |
| `seed_epsilon_padded` | 33.9471743151 | 67.1730404157 | 4.40099191901 | 6.72482943691 | 3.35865202079 | false |

## subset_hardgate_plus_5_6

| Start | Approximate objective | Exact objective | TOA forcing | Surface forcing | Heating RMSE | Passed |
|---|---:|---:|---:|---:|---:|---:|
| `official_normalized` | 31.5493658905 | 50.5675822715 | 7.37858185002 | 9.18356632065 | 2.52837911358 | false |
| `uniform` | 1844.19921697 | 1743.89644356 | 280.766729643 | 489.506659693 | 87.194822178 | false |
| `seed_zero_padded` | 34.130967369 | 67.5356694789 | 4.42626707099 | 6.83468752059 | 3.37678347395 | false |
| `seed_epsilon_padded` | 33.9471778043 | 67.1730152451 | 4.40135966259 | 6.72565801259 | 3.35865076226 | false |

## canonical_plus_7_11

| Start | Approximate objective | Exact objective | TOA forcing | Surface forcing | Heating RMSE | Passed |
|---|---:|---:|---:|---:|---:|---:|
| `official_normalized` | 31.4044803158 | 60.6963087702 | 5.51386117679 | 8.66775276885 | 3.03481543851 | false |
| `uniform` | 22387.3666258 | 19299.8214371 | 323.92118916 | 518.409596085 | 964.991071856 | false |
| `seed_zero_padded` | 142.875041632 | 137.261056182 | 16.099664444 | 41.1783168545 | 4.79533095661 | false |
| `seed_epsilon_padded` | 142.392369803 | 136.802158865 | 16.0624670824 | 41.0406476596 | 4.74076136213 | false |

## subset_hardgate_plus_1

| Start | Approximate objective | Exact objective | TOA forcing | Surface forcing | Heating RMSE | Passed |
|---|---:|---:|---:|---:|---:|---:|
| `official_normalized` | 31.6753543289 | 60.7530603356 | 4.90953518784 | 8.41748364323 | 3.03765301678 | false |
| `uniform` | 1890.27176572 | 1855.50344733 | 296.408264016 | 497.419153467 | 92.7751723667 | false |
| `seed_zero_padded` | 34.1345424911 | 67.5430178292 | 4.42665560546 | 6.83699368413 | 3.37715089146 | false |
| `seed_epsilon_padded` | 34.0372947481 | 67.3515228257 | 4.40892887995 | 6.7703983813 | 3.36757614129 | false |

## canonical_plus_7

| Start | Approximate objective | Exact objective | TOA forcing | Surface forcing | Heating RMSE | Passed |
|---|---:|---:|---:|---:|---:|---:|
| `official_normalized` | 54.7127286855 | 72.8234511377 | 6.96582209781 | 15.2870501281 | 3.64117255688 | false |
| `uniform` | 23706.4408201 | 20438.4004058 | 344.13980766 | 550.658968455 | 1021.92002029 | false |
| `seed_zero_padded` | 142.875166131 | 137.26117443 | 16.099672934 | 41.178352329 | 4.79534740564 | false |
| `seed_epsilon_padded` | 142.599750025 | 136.999129643 | 16.0766092725 | 41.099738893 | 4.76556602169 | false |

