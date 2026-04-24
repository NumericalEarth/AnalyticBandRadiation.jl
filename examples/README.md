# Examples

Standalone scripts that exercise `AnalyticBandRadiation` in richer settings
than the docs examples. They live in their own Julia environment
(`examples/Project.toml`) so the main package stays dep-light.

## `rrtmgp_comparison.jl`

Compare Williams (2026) clear-sky longwave heating rates against RRTMGP on a
shared idealized tropical column.

Breeze's `RadiativeTransferModel(grid, ClearSkyOptics(), …)` provides the
RRTMGP reference; we extract the same `T(z)`, `qᵥ(z)`, and reference-state
pressure profile and feed them to
[`AnalyticBandRadiation.solve_longwave!`](../src/longwave/williams_longwave.jl).

### Run it

```bash
julia --project=examples examples/rrtmgp_comparison.jl
```

The first run will download RRTMGP lookup-table artifacts via `NCDatasets`.
After that the script prints OLR for both schemes and saves a figure next
to itself.

### What to expect

- OLR agreement within ~10–20 W/m² (Williams is a 41-band clear-sky scheme
  with analytic absorption coefficients; RRTMGP is full-spectrum correlated-k).
- Qualitatively similar heating-rate profile: cooling dominant in the
  mid-troposphere where water vapour is abundant, weak cooling above.
- Williams under-cools in the stratosphere by design: it has no ozone line,
  and the CO₂ 15 μm band is a single Lorentzian rather than RRTMGP's full
  correlated-k g-point stack.

The point of the comparison is *not* bit-equivalence — it's to confirm
Williams stays within the factor-2 band of a reference line-by-line-derived
scheme, which is the claim of the paper.
