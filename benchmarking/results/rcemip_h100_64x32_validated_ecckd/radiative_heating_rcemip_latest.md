# RadiativeHeating RCEMIP benchmark

- backend: H100
- grid: 512 × 512 × 128 (262144 columns)
- samples: 5 (post-warmup median)
- gas model: official_ecCKD_64_lw_32_sw (validated_ecCKD)
- ng_lw / ng_sw: 64 / 32

- RadiativeHeating update median: 1352.696 ms
- RRTMGP update median: 5869.765 ms
- speedup vs RRTMGP: 4.34x
