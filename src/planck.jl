"""
    planck_wavenumber(T, ν̃)

Spectral Planck radiance at temperature `T` [K] and wavenumber `ν̃` [cm⁻¹],
in units of W m⁻² sr⁻¹ (cm⁻¹)⁻¹.

The two-stream hemispherical flux source term is `π × planck_wavenumber(T, ν̃)`.
Integrating `π × planck_wavenumber(T, ν̃)` over all wavenumbers recovers the
Stefan–Boltzmann law σ T⁴.

Constants follow CODATA 2018.
"""
@inline function planck_wavenumber(T::NF, ν̃::NF) where NF
    h   = NF(6.62607015e-34)   # Planck constant  [J s]
    c   = NF(2.99792458e8)     # speed of light   [m s⁻¹]
    k_B = NF(1.380649e-23)     # Boltzmann        [J K⁻¹]
    ν_m = ν̃ * 100              # cm⁻¹ → m⁻¹
    # Radiance per-m⁻¹, multiply by 100 to convert to per-cm⁻¹.
    return NF(100) * 2 * h * ν_m^3 * c^2 / (exp(h * c * ν_m / (k_B * T)) - 1)
end

@inline planck_wavenumber(T, ν̃) = planck_wavenumber(promote(T, ν̃)...)
