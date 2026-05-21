# Reduced ecCKD Importance-Guided Group Scan

Status: **all_importance_groups_failed**

This diagnostic keeps high-impact leave-one-out shortwave g-points as
singleton bins, merges the remaining g-points into 16 total groups,
then refits nonnegative shortwave group weights against the normalized
hard-gate objective.

Importance source: `validation/results/reduced_ecckd_leave_one_out_scan.json`

Grouping rule: critical_singleton_kN keeps the N largest leave-one-out objective g-points as singleton bins, then splits all remaining g-points into 16 - N spectral-order groups; criticalN_redundantM also keeps the M smallest leave-one-out objective g-points as singleton bins before splitting the rest.

| Label | Passed | Initial objective | Refit objective | Worst TOA forcing | Worst surface forcing | Groups |
|---|---:|---:|---:|---:|---:|---|
| critical_singleton_k4 | false | 391.107002554 | 391.01556279 | 50.7899217381 W m^-2 | 117.304668837 W m^-2 | `[1,2] [3,4,5] [6,7] [8,9] [10] [11,12,13] [14,15] [16,17] [18,19,20] [21,22] [23,24] [25,29,30] [26] [27] [28] [31,32]` |
| critical_singleton_k6 | false | 413.146663935 | 413.13980092 | 51.3372046246 W m^-2 | 123.941940276 W m^-2 | `[1,2,3] [4,5] [6,8,9] [7] [10] [11,12] [13,14,15] [16,17,18] [19,20] [21,22,23] [24,25] [26] [27] [28] [29,30,31] [32]` |
| critical_singleton_k8 | false | 346.438622867 | 346.425547582 | 49.8364556994 W m^-2 | 103.927664275 W m^-2 | `[1,2,3] [4] [5,6,8] [7] [9,11,12] [10] [13,14,15] [16] [17,18,19] [20,21,22] [23,24,25] [26] [27] [28] [29,30,31] [32]` |
| critical_singleton_k10 | false | 332.59943837 | 332.347531088 | 34.095542173 W m^-2 | 99.7042593264 W m^-2 | `[1,2,3,5] [4] [6,8,9] [7] [10] [11,12,14,15] [13] [16] [17,18,19,20] [21,22,23] [24,25,30,31] [26] [27] [28] [29] [32]` |
| critical_singleton_k12 | false | 278.461302178 | 278.21969455 | 30.655727141 W m^-2 | 83.4659083651 W m^-2 | `[1,2,3,5,6] [4] [7] [8,9,11,15,17] [10] [12] [13] [14] [16] [18,19,20,21,22] [23,24,25,30,31] [26] [27] [28] [29] [32]` |
| critical8_redundant2_singletons | false | 443.355425041 | 443.203681829 | 72.6849606946 W m^-2 | 132.961104549 W m^-2 | `[1,2,3,5] [4] [6,8,9] [7] [10] [11,12,13,14] [15,18,19,20] [16] [17] [21,22,24] [23] [25,29,30,31] [26] [27] [28] [32]` |
| critical8_redundant4_singletons | false | 431.531214424 | 431.476419407 | 60.5595360858 W m^-2 | 129.442925822 W m^-2 | `[1,2,3,5,6] [4] [7] [8,9,11,12,13] [10] [14,15,18,19,20] [16] [17] [21,25,29,30,31] [22] [23] [24] [26] [27] [28] [32]` |
| critical10_redundant2_singletons | false | 318.815769048 | 318.661828063 | 32.2751067957 W m^-2 | 95.598548419 W m^-2 | `[1,2,3,5,6] [4] [7] [8,9,11,12,14] [10] [13] [15,18,19,20,21] [16] [17] [22,24,25,30,31] [23] [26] [27] [28] [29] [32]` |
| critical10_redundant4_singletons | false | 345.337989704 | 345.281587661 | 36.7900457857 W m^-2 | 103.584476298 W m^-2 | `[1,2,3,5,6,8,9,11,12] [4] [7] [10] [13] [14,15,18,19,20,21,25,30,31] [16] [17] [22] [23] [24] [26] [27] [28] [29] [32]` |

Best grouped candidate: `critical_singleton_k12` with objective `278.21969455`.
