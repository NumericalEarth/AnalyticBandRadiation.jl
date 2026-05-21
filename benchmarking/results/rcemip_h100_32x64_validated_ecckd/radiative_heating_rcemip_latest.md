# RadiativeHeating RCEMIP benchmark

- backend: H100
- grid: 512 × 512 × 128 (262144 columns)
- samples: 5 (post-warmup median)
- gas model: official_ecCKD_32_lw_64_sw (validated_ecCKD)
- ng_lw / ng_sw: 32 / 64

- RadiativeHeating update median: 1047.307 ms
- RRTMGP update median: 5864.235 ms
- speedup vs RRTMGP: 5.60x
