# Future work

This note scopes two extensions: an **all-sky longwave** variant of the
Williams (2026) scheme and a **two-band shortwave** variant of the current
SPEEDY one-band shortwave.

## All-sky longwave

The Williams (2026) Simple Spectral Model is defined for clear skies only.
Its optical-depth increment [`williams_delta_tau`](src/longwave/williams_absorption.jl)
sums H₂O line, H₂O continuum, and CO₂ contributions; nothing for clouds.

Three levels of cloud coupling, roughly in order of effort:

### 1. Gray emissivity overlay (least effort, recommended starting point)

Treat each layer's cloud as a gray body with emissivity

    ε_c(k) = 1 − exp(−κ_LW · LWP_k)

where `LWP_k = ρ q_l Δz` is the layer liquid water path (plus an analogous
ice term with its own mass absorption coefficient). Blend clear/cloudy fluxes
by the layer cloud fraction:

    Δτ_layer = Δτ_gas(ν̃)     with probability (1 − c_k)
    Δτ_layer = Δτ_gas(ν̃) + Δτ_cloud(ν̃, LWP, IWP)  with probability c_k

and propagate *either one* pair of upward/downward sweeps over the
fractional-weighted layer transmittances. This is the approach taken by the
SPEEDY shortwave already (`OneBandShortwave` treats cloud as a reflector at
`cloud_top`), and by intermediate-complexity schemes generally.

**New inputs**:

- Per-layer cloud fraction, liquid water content, ice water content.
- Two new absorption constants in [`WilliamsLongwave`](src/longwave/williams_longwave.jl):
  `κ_cloud_liquid`, `κ_cloud_ice` — wavenumber-independent if we stay
  "simple", or wavenumber-dependent piecewise constants (split at the
  8–12 μm window; ice much more absorbing in the window, liquid much more
  absorbing outside).

**Where it plugs in**:

- Extend the new optical-depth signature:
  `williams_delta_tau(k, ν̃, temperature, humidity, cloud_liquid, cloud_ice,
   cloud_fraction, surface_pressure, geometry, scheme, gravity)`
- Extend `ColumnProfile` with `cloud_liquid`, `cloud_ice`, `cloud_fraction`
  vectors (optional; clear-sky callers pass `fill(0, nlayers)`).
- Add `κ_cloud_liquid`, `κ_cloud_ice` constants to the scheme; wavenumber-split
  them via a helper analogous to `h2o_cont_kappa_ref`.
- Tests: low-cloud-over-warm-surface case must have reduced OLR and near-zero
  downward surface flux reduction relative to clear sky (cloud radiative
  effect > 0 W/m²); sign of TOA cloud forcing depends on surface vs. cloud-top
  temperature.

**Effort**: ~2 days for a working implementation + tests. The column sweeps
need no structural change — `exp(-Δτ)` just includes the cloud contribution.

### 2. Maximum-random overlap

Cloud layers are not independent between levels. The standard fix is the
Geleyn-Hollingsworth (1979) maximum-random overlap: contiguous cloudy layers
overlap maximally; separated cloudy layers overlap randomly. Implementable
in the upward/downward sweeps by splitting each spectral iteration into
(clear, cloudy) streams and combining at the end.

**Effort**: ~1 week on top of (1). Adds an outer "stream" loop of length 2
around the wavenumber loop, effectively doubling the cost. Brings the scheme
to parity with schemes like Morcrette (RRTMG) for the overlap treatment.

### 3. Scattering two-stream

Upgrade Schwarzschild to a full two-stream that carries the single-scattering
albedo ω and asymmetry factor g (δ-Eddington, practical pattern from
Toon et al. 1989). For LW this is usually unnecessary (scattering by
liquid/ice is ≪ absorption at infrared wavelengths) but it is required for
bright stratiform ice clouds and for the shortwave (see below).

**Effort**: ~2 weeks. The column solver becomes a tridiagonal system per
spectral bin; no longer a simple running-scalar sweep.

**Recommendation**: ship (1) first. It covers the dominant greenhouse-like
cloud effect; (2) can follow once we need it, and (3) is overkill for
intermediate-complexity modelling.

---

## Two-band shortwave

The current [`OneBandShortwave`](src/shortwave/one_band_shortwave.jl) is
derived from Fortran SPEEDY with the visible/near-IR band weights (`fband =
0.95 / 0.05`) already folded into the absorption constants at source time.
Splitting back into two bands restores physical realism in two places:

- **Cloud reflection** is dominated by the visible band; cloud absorption is
  dominated by the near-IR band.
- **Water vapour** barely absorbs visible light but absorbs strongly in the
  near-IR.
- Surface albedo is much higher in the visible over snow/ice — a two-band
  split lets us correctly represent snow-albedo feedback.

### Implementation plan

**New scheme type**:

```julia
struct TwoBandShortwave{C, TV, TNIR, R} <: AbstractShortwaveScheme
    clouds::C
    visible_transmissivity::TV
    near_infrared_transmissivity::TNIR
    radiative_transfer::R
end
```

**Band partition**: `f_visible = 0.521`, `f_near_ir = 0.479` (Stephens 1978)
or the SPEEDY `fband = 0.95 / 0.05` split (matches the current scheme's
collapsed tuning). Make the fraction a scheme parameter.

**Per-band transmissivity**:

- Visible: `a_dry ≈ 0`, `a_wv ≈ 0`, `a_cloud ≈ 0.02–0.05` (cloud droplets
  barely absorb in the visible).
- Near-IR: `a_dry ≈ 0.03`, `a_wv ≈ 15` (per kg/kg per 10⁵ Pa),
  `a_cloud ≈ 0.3` (droplets are strong near-IR absorbers).

These are the Fortran SPEEDY `absdry`, `abswv1/abswv2`, `abscl1/abscl2`
pre-weighting values, which we can recover from the `0.95/0.05` combined
weights stored in the current `BackgroundShortwaveTransmissivity` defaults.

**Per-band cloud albedo**: visible ≈ 0.7, near-IR ≈ 0.2. Expose as scheme
parameters rather than a single `cloud_albedo`.

**Solver**: run the current downward/upward sweeps twice, sum the fluxes and
the tendency contributions. Cost is ~2× the current scheme; still trivially
cheap compared to LW (which is 41×).

**Surface albedo**: accept `visible_albedo` and `near_infrared_albedo` in
[`ColumnSurface`](src/column_views.jl) (or keep broadband and internally
split 0.7/0.3 over snow for starters).

**Effort**: ~3 days including tests for (a) matching the current `OneBandShortwave`
output when fed the band-merged parameters, (b) increased SW at surface over
snow when two-band is active, (c) cloud-albedo response larger in visible.

### What NOT to do

- Don't add ozone bands: the current single scalar `ozone_absorption = 0.01`
  is already a 1-band approximation for a concentrated stratospheric UV
  sink. A visible-band ozone (Chappuis) would add ~0.02 ppm-scale detail
  not worth the parameter.
- Don't add Rayleigh scattering as a separate band: for a simple scheme, its
  effect (~6% of TOA visible flux) is absorbed into the "visible cloud"
  albedo tuning.

---

## Summary table

| Extension | Effort | New inputs | New scheme fields | New tests |
|---|---|---|---|---|
| Gray-emissivity all-sky LW | ~2 days | cloud_fraction, LWP, IWP per layer | `κ_cloud_liquid`, `κ_cloud_ice` (optionally ν̃-split) | LW cloud radiative effect > 0 |
| Max-random overlap | +1 week | (same) | cloud-overlap toggle | Stream-combined fluxes = independent-layers fluxes when cloud fraction is 0 or 1 |
| Scattering two-stream | +2 weeks | ω, g per band | tridiagonal solver state | Conservative-scattering limit; clear-sky reduces to Schwarzschild |
| Two-band SW | ~3 days | `visible_albedo`, `near_infrared_albedo` (optional) | per-band transmissivity params + cloud albedo params | Snow albedo feedback; band-merged ≡ one-band |

Recommended order: **two-band SW first** (smallest, most widely useful),
then **gray-emissivity all-sky LW** (unlocks cloud-radiation interaction for
the dominant longwave effect).
