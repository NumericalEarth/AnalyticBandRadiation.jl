# RadiativeHeating RCEMIP benchmark

- backend: H100
- grid: 512 × 512 × 128 (262144 columns)
- samples: 5 (post-warmup median)
- gas model: official_ecCKD_64_lw_64_sw (validated_ecCKD)
- ng_lw / ng_sw: 64 / 64

- RadiativeHeating update median: 1603.626 ms
- RRTMGP update median: 5870.138 ms
- speedup vs RRTMGP: 3.66x
