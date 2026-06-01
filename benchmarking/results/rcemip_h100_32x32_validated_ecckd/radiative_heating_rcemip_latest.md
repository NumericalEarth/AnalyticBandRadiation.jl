# RadiativeHeating RCEMIP benchmark

- backend: H100
- grid: 512 × 512 × 128 (262144 columns)
- samples: 5 (post-warmup median)
- gas model: official_ecCKD_32_lw_32_sw (validated_ecCKD)
- ng_lw / ng_sw: 32 / 32

- RadiativeHeating update median: 794.743 ms
- RRTMGP update median: 5863.345 ms
- speedup vs RRTMGP: 7.38x
