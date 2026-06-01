# Reduced ecCKD Leave-One-Out Weight Coordinate Scan

Status: **weight_coordinate_scan_improved**

| Field | Value |
|---|---:|
| Omitted SW g-point | 23 |
| Model | 32x31 |
| Coordinate count | 16 |
| Exact rows | 128 |
| Boundary cap | 0.3 W m^-2 |
| Base objective | 1.62442314379 |
| Base boundary forcing | 0.0751392004466 W m^-2 |
| Base heating RMSE | 0.0812211571893 K day^-1 |
| Best objective | 1.47156896537 |
| Best objective reduction | 0.152854178418 |
| Best boundary forcing | 0.441470689611 W m^-2 |
| Best heating RMSE | 0.0733566893747 K day^-1 |
| Accepted | true |
| Accepted objective | 1.48556457006 |
| Accepted boundary forcing | 0.254747858736 W m^-2 |
| Accepted heating RMSE | 0.074278228503 K day^-1 |

This diagnostic exact-evaluates normalized one-coordinate weight moves
around the objective-best 32x31 leave-one-out weight-refit row. It
tests whether the approximate hard-gate weight refit left any simple
exact descent direction before changing support.

## Best Exact Rows

| Rank | G-point | Base weight | Log scale | Objective | Reduction | Boundary forcing | Heating RMSE | Accepted |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 29 | 0.0647002790773 | -0.0316227766017 | 1.47156896537 | 0.152854178418 | 0.441470689611 | 0.0733566893747 | false |
| 2 | 28 | 0.139581641199 | -0.01 | 1.48556457006 | 0.138858573726 | 0.254747858736 | 0.074278228503 | true |
| 3 | 3 | 0.043544717568 | -0.0316227766017 | 1.48708726537 | 0.137335878416 | 0.290531753936 | 0.0743543632685 | true |
| 4 | 5 | 0.0424649479649 | -0.0316227766017 | 1.49014993449 | 0.134273209299 | 0.260077221292 | 0.0745074967244 | true |
| 5 | 10 | 0.13176827178 | -0.01 | 1.49216192069 | 0.132261223099 | 0.304045589656 | 0.0746080960344 | false |
| 6 | 27 | 0.135698558686 | -0.01 | 1.5043461086 | 0.120077035186 | 0.278278499068 | 0.07521730543 | true |
| 7 | 11 | 0.0351984879441 | -0.0316227766017 | 1.51534487091 | 0.109078272872 | 0.266670468128 | 0.0757672435457 | true |
| 8 | 7 | 0.058483331877 | -0.0316227766017 | 1.51830758923 | 0.106115554561 | 0.455492276768 | 0.0721775136923 | false |
| 9 | 26 | 0.114114524831 | -0.01 | 1.52265615647 | 0.10176698732 | 0.230794777085 | 0.0761328078233 | true |
| 10 | 1 | 0.0264374455328 | -0.0316227766017 | 1.54167878332 | 0.0827443604698 | 0.0790966612042 | 0.0770839391658 | true |
| 11 | 8 | 0.0239363472401 | -0.0316227766017 | 1.54874901987 | 0.0756741239171 | 0.272192850216 | 0.0774374509935 | true |
| 12 | 6 | 0.0220986304147 | -0.0316227766017 | 1.55486576938 | 0.0695573744049 | 0.392964797056 | 0.0777432884691 | false |
