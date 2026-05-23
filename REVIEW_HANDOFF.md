# Radiation Goal Handoff

## Current status

Focus has shifted to completing the ecCKD recovery objective side in this repo without touching Breeze code paths.

Latest audit and requirement tracking live in:
- `validation/results/recovery_goal_audit.md`
- `validation/results/recovery_goal_audit.json`
- `validation/recovery_goal_audit.jl`

`validation/results/recovery_goal_audit.md` currently reports:
- Teacher-student recovery status: passed
- Objective reconstruction status: `ready_to_reconstruct_original_objective`
- Official training data preflight status: `ready_for_original_ecckd_objective`
- Partial requirements: 3
- Unmet requirements: 3
- Training recovery target status: `partial`

## What is implemented in the codebase

1. New humidity-aware reduced-SW table optimization workflows were added:
   - probe
   - aggregate descent
   - minimax descent (standalone artifact)
   - weighted descent
   - constrained descent

2. Supporting files added
   - scripts:
     - `validation/ecckd_candidate_table_humidity_tribin_probe.jl`
     - `validation/ecckd_candidate_table_humidity_tribin_descent.jl`
     - `validation/ecckd_candidate_table_humidity_tribin_weighted_descent.jl`
     - `validation/ecckd_candidate_table_humidity_tribin_constrained_descent.jl`
   - tests:
     - `test/test_ecckd_candidate_table_humidity_tribin_probe.jl`
     - `test/test_ecckd_candidate_table_humidity_tribin_descent.jl`
     - `test/test_ecckd_candidate_table_humidity_tribin_weighted_descent.jl`
     - `test/test_ecckd_candidate_table_humidity_tribin_constrained_descent.jl`
   - runner wiring:
     - `test/runtests.jl`
     - `test/test_recovery_goal_audit.jl`
     - `validation/recovery_goal_audit.jl`

3. Artifacts produced in `validation/results/`:
   - `ecckd_candidate_table_humidity_tribin_probe.json/.md`
   - `ecckd_candidate_table_humidity_tribin_descent.json/.md`
   - `ecckd_candidate_table_humidity_tribin_weighted_descent.json/.md`
   - `ecckd_candidate_table_humidity_tribin_constrained_descent.json/.md`
   - corresponding candidate table files: `ecckd_table_*_sw32_candidate.nc`

4. PR-facing documentation was updated during the run:
   - `PR_WORK_SUMMARY.md`
   - `RUNNING_REVIEW.md`

## Numeric outcome to date

- Probe step improved aggregate by `1.0003671142776498` and improved worst-ratio by `1.0000436882331465`.
- Aggregate descent improved aggregate by `1.0016410217592027`; worst-ratio worsened slightly.
- Weighted descent reproduced the aggregate descent outcome using scalar humidity-weight variants.
- Constrained descent did not accept any move under no-worsening strict policy, then selected tolerance-aware variants that improve aggregate but still violate the current acceptance policy when constrained.
- Final constrained status in the audit is currently marked as `humidity_tribin_constrained_descent_no_descent`.

## Blocking items (current)

1. Hard recovery objective remains far from target for the reduced training objective despite optimizer probes/iterative refinement.
2. Constraint handling in humidity-sensitive minimax balance remains unresolved; aggregate improvements continue to trade off against worst-scenario ratio.
3. The final remaining work is to implement a stronger optimization strategy while keeping objective/optimizer settings as required.

## Immediate next action

Run/continue training-recovery pipeline with a joint optimization strategy (e.g., tighter humidity/pressure-aware parameterization or multi-objective weighted schedule) and feed its outputs into:
- `validation/ecckd_official_training.jl` path
- `validation/results/{ecckd_original_objective_terms, ckdmip_original_objective_*}.json`
- `validation/results/recovery_goal_audit.json`

Then update:
- `validation/recovery_goal_audit.jl` with new status lines
- `validation/results/recovery_goal_audit.md`
- `validation/results/recovery_goal_audit.json`
- `PR_WORK_SUMMARY.md`
- `RUNNING_REVIEW.md`

## No-Breeze constraint

Do not touch any Breeze repository files for this pass. Continue all changes in this repo only:
- `validation/*`
- `test/*`
- documentation files in this repo

## Recommended immediate implementation sequence

1. Add a constrained optimization mode that optimizes a combined scalar objective for aggregate loss and worst-ratio penalty (soft penalty rather than hard no-worsening gate).
2. Re-run constrained experiment with a conservative penalty schedule and persist:
   - `validation/results/ecckd_candidate_table_humidity_tribin_constrained_descent.json`
   - `validation/results/ecckd_candidate_table_humidity_tribin_constrained_descent.md`
3. Re-run audit wiring checks in `validation/recovery_goal_audit.jl` and update the partial/goal rows only after objective status changes.
