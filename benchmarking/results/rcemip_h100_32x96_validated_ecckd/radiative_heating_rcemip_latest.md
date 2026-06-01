# RadiativeHeating RCEMIP benchmark

- backend: H100
- grid: 512 × 512 × 128 (262144 columns)
- samples: 5 (post-warmup median)
- gas model: official_ecCKD_32_lw_96_sw (validated_ecCKD)
- ng_lw / ng_sw: 32 / 96

- RadiativeHeating update median: 1279.349 ms
- RRTMGP update median: 5867.678 ms
- speedup vs RRTMGP: 4.59x
