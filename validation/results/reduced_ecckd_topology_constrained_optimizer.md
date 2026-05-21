# Reduced ecCKD Topology Constrained Optimizer

Status: **failed_threshold**

This diagnostic applies the constrained active-entry table optimizer to
alternate 16-g shortwave topologies from the subset-search artifact.
It is separate from the canonical reduced model artifact and does not
overwrite the accepted canonical constrained moves.

| Topology | Accepted | Objective | TOA forcing error | Surface forcing error |
|---|---:|---:|---:|---:|
| subset_search_official_weight_greedy | true | 25.007407771 | 7.34730685554 W m^-2 | 6.97942314031 W m^-2 |
| subset_search_boundary_weight_30 | false | 23.3015949902 | 6.99047849707 W m^-2 | 6.84847043029 W m^-2 |
| subset_search_boundary_weight_10 | true | 24.7122008243 | 7.41366024729 W m^-2 | 7.38001611787 W m^-2 |
