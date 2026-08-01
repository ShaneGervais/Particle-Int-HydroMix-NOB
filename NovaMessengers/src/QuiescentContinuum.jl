module QuiescentContinuum

using DataFrames
using ..MesaIO

export PlanckSource, wd_quiescent_source, spectral_luminosity_ev, photon_energy_grid_ev

const H_ERG_S = 6.62607015e-27     # erg s (exact, 2019 SI)
const K_B_ERG_K = 1.380649e-16     # erg/K (exact, 2019 SI)
const C_CM_S = 2.99792458e10       # cm/s (exact)
const EV_TO_ERG = 1.602176634e-12  # erg/eV (exact)

"""
    PlanckSource(teff_K, radius_cm)

A single blackbody emitter: effective temperature (K) and radius (cm).
Both the WD ([`wd_quiescent_source`](@ref)) and the companion (an
external parameter, see that function's docstring) are represented this
way -- the simplest possible photosphere model, appropriate for a
baseline "what does the system look like before anything interesting
happens" continuum, not a substitute for real stellar atmosphere
modeling.
"""
struct PlanckSource
    teff_K::Float64
    radius_cm::Float64
end

"""
    wd_quiescent_source(run::MesaRun; history_row=1) -> PlanckSource

The WD's own quiescent photosphere, read directly from MESA's history
data (`log_Teff`, `radius_cm`) -- unlike the companion, this is not an
external assumption, since MESA already models the WD surface. Defaults
to `history_row=1` (the first saved model, i.e. before the TNR); pass a
later row to see how the WD's own photosphere evolves (e.g. its
post-burst re-expansion, before ejecta-related effects dominate).
"""
function wd_quiescent_source(run::MesaRun; history_row::Integer=1)
    h = history(run).data
    teff_K = 10.0^h.log_Teff[history_row]
    radius_cm = h.radius_cm[history_row]
    return PlanckSource(teff_K, radius_cm)
end

"""
    spectral_luminosity_ev(source::PlanckSource, E_ev::Real) -> Float64

Blackbody spectral luminosity density (erg/s/eV) at photon energy `E_ev`
(eV), from the Planck function. In frequency space,

    B_nu(T) = (2 h nu^3 / c^2) / (exp(h nu / kT) - 1)         [erg/s/cm^2/sr/Hz]
    L_nu    = 4 pi^2 R^2 * B_nu(T)                             [erg/s/Hz]

(flux = pi*B_nu integrated over the outward hemisphere, times surface
area 4 pi R^2 -- integrating `L_nu` over all `nu` recovers the
Stefan-Boltzmann law `L = 4 pi R^2 sigma T^4`, a useful self-check).
Substituting `E = h*nu` and converting `d(erg)/d(eV)`:

    L_E(E) = 8 pi^2 R^2 E^3 / (h^3 c^2 [exp(E/kT) - 1]) * EV_TO_ERG^2

(one factor of `EV_TO_ERG` from `E^3` being evaluated in erg inside the
Planck form, one from converting the per-erg result to per-eV).
"""
function spectral_luminosity_ev(source::PlanckSource, E_ev::Real)
    E_ev > 0 || error("photon energy must be positive, got $E_ev")
    E_erg = E_ev * EV_TO_ERG
    x = E_erg / (K_B_ERG_K * source.teff_K)
    # expm1, not exp(x)-1: x can be tiny (far-IR/radio against an eV-scale
    # kT) or huge (X-ray tail), and the naive difference loses precision
    # at both ends the way klein_nishina_factor's log1p already guards
    # against in Transport.jl.
    denom = expm1(x)
    prefactor = 8 * pi^2 * source.radius_cm^2 * E_erg^3 / (H_ERG_S^3 * C_CM_S^2)
    return prefactor / denom * EV_TO_ERG
end

"""
    spectral_luminosity_ev(sources::AbstractVector{PlanckSource}, E_ev::Real) -> Float64

Total spectral luminosity density (erg/s/eV) summed over several
independent blackbody sources (e.g. WD + companion) at the same energy
-- valid because blackbody emission from non-overlapping, non-interacting
photospheres simply adds.
"""
function spectral_luminosity_ev(sources::AbstractVector{PlanckSource}, E_ev::Real)
    return sum(spectral_luminosity_ev(s, E_ev) for s in sources)
end

"""
    photon_energy_grid_ev(E_min_ev, E_max_ev; n=200) -> Vector{Float64}

Log-spaced photon-energy grid (eV) spanning `[E_min_ev, E_max_ev]` --
the natural sampling for a blackbody continuum viewed across many
decades of energy (radio through X-ray), matching the log-spaced
quadrature already used in `ShockAcceleration._log_trapz` for the same
reason (a smooth function spanning many decades, not a locally-varying
one a linear grid would resolve better).
"""
function photon_energy_grid_ev(E_min_ev::Real, E_max_ev::Real; n::Integer=200)
    return exp.(range(log(E_min_ev), log(E_max_ev); length=n))
end

end # module QuiescentContinuum
