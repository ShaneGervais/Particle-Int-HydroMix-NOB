module ExtendedMessengers

using CSV
using DataFrames
using ..MesaIO
using ..NuclearDecay
using ..Transport
using ..TrajectoryPostProcessing: _linterp

export read_mass_fraction_history, extended_decay_rate_per_gram,
       extended_positron_rate_per_gram, extended_gamma_rate_per_gram,
       freezeout_gamma_lightcurve, freezeout_gamma_energy_lightcurve

const AMU_G = 1.66053906892e-24  # grams (CODATA)
const MEV_TO_ERG = 1.602176634e-6
const YEAR_S = 365.25 * 24 * 3600

"""
    read_mass_fraction_history(path) -> DataFrame

Parse ReacNetJl's own `mass_fractions.csv` output (written by
`ReacNetJl.run_ppn` / `TrajectoryPostProcessing.postprocess_trajectory`):
one row per solver-accepted time, columns `time_s, T9, rho, <species>...`.
Column names already match this project's own isotope-Symbol convention
(e.g. `na22`, `be7`), so no name translation is needed. `time_s` here is
relative to the trajectory's own t=0 (the pristine epoch, see
[`TrajectoryPostProcessing.postprocess_trajectory`](@ref)), not
`star_age` -- callers that need to line this up with MESA's other light
curves must rebase it (see [`freezeout_gamma_lightcurve`](@ref)'s
`t0_star_age_s`).
"""
function read_mass_fraction_history(path::AbstractString)
    return CSV.read(path, DataFrame)
end

"""
    extended_decay_rate_per_gram(mass_fraction_history, isotope::Symbol) -> (time_s, rate)

Decay rate PER GRAM of the tracked zone's material (decays/second/gram),
from ReacNetJl's own time-resolved mass fraction X(t):

    rate(t) = lambda * X(t) / (A * AMU_G)

Deliberately "per gram," not a whole-star total: unlike
[`MessengerProduction.decay_rate`](@ref) (which sums MESA's own
`total_mass_<iso>` over the *whole star*), this comes from a single-zone
ReacNetJl trajectory that only knows the *fraction* of that one zone's
material, not how many grams of the star that zone represents.
Converting to an actual photon/neutrino count needs an assumed total
processed mass -- see [`freezeout_gamma_lightcurve`](@ref)'s
`total_mass_g` argument, which makes that assumption an explicit,
visible input rather than a hidden default.
"""
function extended_decay_rate_per_gram(mass_fraction_history::AbstractDataFrame, isotope::Symbol)
    col = isotope
    col in propertynames(mass_fraction_history) || error("mass fraction history has no column $col")
    iso = EXTENDED_ISOTOPES[isotope]
    lambda = decay_constant(iso)
    A = mass_number(isotope)
    N_per_gram = mass_fraction_history[!, col] ./ (A * AMU_G)
    return (time_s=mass_fraction_history.time_s, rate=lambda .* N_per_gram)
end

"""
    extended_positron_rate_per_gram(mass_fraction_history, isotope) -> (time_s, rate)

Positron production rate per gram:

    rate_e+(t) = decay_rate_per_gram(t) * positron_branching

NOT `= decay_rate_per_gram(t)` the way [`MessengerProduction.positron_rate`](@ref)
assumes for the 7 (idealized pure beta+) isotopes. This is exactly the
correction be7 needs: `positron_branching = 0`, so it produces zero
positrons -- pure electron capture, not "1 per decay."
"""
function extended_positron_rate_per_gram(mass_fraction_history::AbstractDataFrame, isotope::Symbol)
    r = extended_decay_rate_per_gram(mass_fraction_history, isotope)
    br = EXTENDED_ISOTOPES[isotope].positron_branching
    return (time_s=r.time_s, rate=r.rate .* br)
end

"""
    extended_gamma_rate_per_gram(mass_fraction_history, isotope) -> (time_s, rate, line_energy_mev)

Prompt daughter-de-excitation gamma-line production rate per gram:

    rate_gamma(t) = decay_rate_per_gram(t) * gamma_branching

at the isotope's characteristic `gamma_line_mev` -- a physically distinct
channel from e+e- annihilation (511 keV), not routed through
[`MessengerProduction.annihilation_photon_rate`](@ref).
"""
function extended_gamma_rate_per_gram(mass_fraction_history::AbstractDataFrame, isotope::Symbol)
    r = extended_decay_rate_per_gram(mass_fraction_history, isotope)
    iso = EXTENDED_ISOTOPES[isotope]
    return (time_s=r.time_s, rate=r.rate .* iso.gamma_branching, line_energy_mev=iso.gamma_line_mev)
end

"""
    freezeout_gamma_lightcurve(run::MesaRun, mass_fraction_history::AbstractDataFrame,
                                isotope::Symbol, mass_coordinate::Real, total_mass_g::Real;
                                t0_star_age_s::Real) -> (age, rate)

Observable freeze-out gamma-line escape rate (photons/second) for
`isotope` (na22's 1.275 MeV or al26's 1.809 MeV line):

    rate(t_p) = gamma_rate_per_gram(t_p - t0) * total_mass_g * P_escape(mass_coordinate, p)

evaluated at each of MESA's saved profile snapshots `p` (matching
[`SignalSynthesis.gamma_lightcurve`](@ref)'s cadence), where:

- `gamma_rate_per_gram` ([`extended_gamma_rate_per_gram`](@ref)) is
  linearly interpolated from ReacNetJl's own (generally irregular,
  solver-driven) time grid onto each profile's absolute `star_age`
  (converted via `t0_star_age_s`, the pristine epoch's own `star_age`
  in seconds -- see `TrajectoryPostProcessing.postprocess_trajectory`);
- production is scaled by the caller-supplied `total_mass_g` (grams of
  material assumed to share this trajectory's history) -- an explicit
  approximation, since single-zone post-processing has no built-in
  notion of "how much of the star is like this," unlike the whole-star
  `total_mass_<iso>` bookkeeping the original 7 isotopes use. A
  reasonable starting estimate is the envelope mass above the tracked
  zone (from MESA's own `mass` profile column) or a shock-model shell
  mass (`ShockAcceleration.shell_mass`), depending on which phase the
  observation is for.
- `P_escape` ([`Transport.escape_probability_gamma`](@ref)) is
  evaluated at `mass_coordinate` using MESA's *real* envelope
  density/composition structure at profile `p` -- this part is not an
  approximation the way the mass scaling is, since the envelope
  structure itself doesn't depend on which isotope is decaying in it.

REAL FINDING from running this against the CO WD baseline: `P_escape`
at the tracked (peak-temperature) zone's mass coordinate is
indistinguishable from 0 at *every* saved profile, so this function
currently returns ~0 photons/second throughout the whole run for na22/
al26. This is not a bug -- confirmed by inspecting the profile data
directly: the tracked zone sits ~1e-6 Msun below where escape
probability rises from 0 to 1, and that entire transition happens
within the outermost ~1e-13 Msun of the star (profile 97's last 6
zones, checked directly). That's a *much* stricter test than
`SignalSynthesis.gamma_lightcurve`'s existing "~0 near burst peak"
finding: that function sums production x escape over *every* zone,
so it still picks up whatever the already-transparent outermost sliver
contributes, even while the deep zones (including this one) contribute
~0. A single tracked zone has no such outer-zone contribution to fall
back on. Getting a non-zero freeze-out signal out of this pipeline
needs either a MESA run that actually resolves homologous ejection
(the whole envelope becoming transparent, not just its outermost
sliver -- the same gap `ShockAcceleration.jl` has), or explicitly
modeling na22/al26 transport to the transparent outer layers rather
than assuming it stays at the deep zone where it was made.

Returns `age` in seconds (`star_age`, same convention as
`SignalSynthesis`'s other light curves) so this can be plotted or
composited directly alongside `gamma_lightcurve`/`neutrino_lightcurve`.
"""
function freezeout_gamma_lightcurve(run::MesaRun, mass_fraction_history::AbstractDataFrame,
    isotope::Symbol, mass_coordinate::Real, total_mass_g::Real; t0_star_age_s::Real)
    prod = extended_gamma_rate_per_gram(mass_fraction_history, isotope)

    pidx = profiles_index(run)
    hist = history(run).data
    rsun_cm = history(run).header.rsun

    n = nrow(pidx)
    age = Vector{Float64}(undef, n)
    rate = Vector{Float64}(undef, n)
    for (i, row) in enumerate(eachrow(pidx))
        j = findfirst(==(row.model_number), hist.model_number)
        star_age_s = j === nothing ? NaN : hist.star_age[j] * YEAR_S
        age[i] = star_age_s

        t_relative = star_age_s - t0_star_age_s
        rate_per_gram = _linterp(prod.time_s, prod.rate, t_relative)

        prof = profile(run, row.profile_number).data
        esc = escape_probability_gamma(prof, prod.line_energy_mev; rsun_cm=rsun_cm)
        order = sortperm(prof.mass)
        P_escape = _linterp(prof.mass[order], esc[order], mass_coordinate)

        rate[i] = rate_per_gram * total_mass_g * P_escape
    end
    order = sortperm(age)
    return (age=age[order], rate=rate[order])
end

"""
    freezeout_gamma_energy_lightcurve(run, mass_fraction_history, isotope, mass_coordinate,
                                       total_mass_g; t0_star_age_s) -> (age, rate)

Observable freeze-out line luminosity (erg/second):

    L_gamma(t) = freezeout_gamma_lightcurve(t) * gamma_line_mev * MEV_TO_ERG
"""
function freezeout_gamma_energy_lightcurve(run::MesaRun, mass_fraction_history::AbstractDataFrame,
    isotope::Symbol, mass_coordinate::Real, total_mass_g::Real; t0_star_age_s::Real)
    g = freezeout_gamma_lightcurve(run, mass_fraction_history, isotope, mass_coordinate, total_mass_g;
        t0_star_age_s=t0_star_age_s)
    line_energy_mev = EXTENDED_ISOTOPES[isotope].gamma_line_mev
    return (age=g.age, rate=g.rate .* line_energy_mev .* MEV_TO_ERG)
end

end # module ExtendedMessengers
