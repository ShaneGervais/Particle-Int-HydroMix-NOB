module MessengerProduction

using DataFrames
using ..MesaIO
using ..NuclearDecay

export isotope_number, production_rate, positron_rate, neutrino_energy_loss_rate,
       annihilation_photon_rate, ne18_production_rate

const AMU_G = 1.66053906892e-24  # grams (CODATA)
const YEAR_S = 365.25 * 24 * 3600

"""
    isotope_number(run::MesaRun, isotope::Symbol) -> Vector{Float64}

Number of nuclei of `isotope` in the whole star at each saved history model,
from the `total_mass_<isotope>` column (Msun, converted using this run's own
`msun` header value) and the isotope's mass number (amu approximation).
"""
function isotope_number(run::MesaRun, isotope::Symbol)
    h = history(run)
    col = Symbol("total_mass_", isotope)
    col in propertynames(h.data) || error("history data has no column $col")
    A = mass_number(isotope)
    return (h.data[!, col] .* h.header.msun) ./ (A * AMU_G)
end

"""
    production_rate(run::MesaRun, isotope::Symbol) -> (age, rate)

Reconstruct the production rate (nuclei/second) of `isotope` from
finite-differencing its whole-star inventory in history.data, correcting
for the isotope's own decay: dN/dt = production - lambda*N, so
production = dN/dt + lambda*N.

This is only physically meaningful when the isotope's half-life is resolved
by the saved history cadence (true for all 7 messenger isotopes given
history_interval=1 near burst peak, except ne18 -- its 1.7s half-life can be
shorter than a single hydro timestep, so use [`ne18_production_rate`](@ref)
instead). Returns rates at the midpoints between consecutive saved ages.
"""
function production_rate(run::MesaRun, isotope::Symbol)
    N = isotope_number(run, isotope)
    age = history(run).data.star_age .* YEAR_S
    lambda = decay_constant(MESSENGER_ISOTOPES[isotope])
    return _finite_diff_production(age, N, lambda)
end

"""
    _finite_diff_production(age, N, lambda) -> (age, rate)

Numeric core of [`production_rate`](@ref), decoupled from MESA file I/O so
it can be unit-tested against analytic decay/production solutions directly.
`age` and `N` are whole-star age (seconds) and nuclei-count time series;
`lambda` is the isotope's decay constant (1/second).
"""
function _finite_diff_production(age::AbstractVector, N::AbstractVector, lambda::Real)
    n = length(age)
    age_mid = Vector{Float64}(undef, n - 1)
    rate = Vector{Float64}(undef, n - 1)
    for i in 1:(n - 1)
        dt = age[i + 1] - age[i]
        dNdt = dt > 0 ? (N[i + 1] - N[i]) / dt : 0.0
        Nmid = 0.5 * (N[i] + N[i + 1])
        age_mid[i] = 0.5 * (age[i] + age[i + 1])
        rate[i] = dNdt + lambda * Nmid
    end
    return (age=age_mid, rate=rate)
end

"""
    positron_rate(run::MesaRun, isotope::Symbol) -> (age, rate)

Positron production rate (positrons/second) -- one per beta+ decay, so
identical to [`production_rate`](@ref) for these isotopes (all beta+, see
`NuclearDecay.MESSENGER_ISOTOPES`).
"""
positron_rate(run::MesaRun, isotope::Symbol) = production_rate(run, isotope)

"""
    neutrino_energy_loss_rate(run::MesaRun, isotope::Symbol) -> (age, rate)

Neutrino energy loss rate (MeV/second) from `isotope`'s beta+ decay:
production rate times the average neutrino energy per decay (`q_neu_mev`).
"""
function neutrino_energy_loss_rate(run::MesaRun, isotope::Symbol)
    r = production_rate(run, isotope)
    q = MESSENGER_ISOTOPES[isotope].q_neu_mev
    return (age=r.age, rate=r.rate .* q)
end

"""
    annihilation_photon_rate(run::MesaRun, isotope::Symbol) -> (age, rate)

511 keV annihilation photon production rate (photons/second): two photons
per positron, assuming in-situ annihilation (site-of-annihilation / true
in-flight annihilation is a Phase 3 transport question, deferred here).
"""
function annihilation_photon_rate(run::MesaRun, isotope::Symbol)
    r = positron_rate(run, isotope)
    return (age=r.age, rate=2 .* r.rate)
end

"""
    ne18_production_rate(run::MesaRun) -> (age, rate)

ne18's 1.7s half-life is shorter than a hydro timestep near burst peak, so
finite-differencing its total_mass is not physically meaningful (it's in
secular equilibrium within a single saved step). Instead, sum the screened
f17(p,g)ne18 reaction rate (`screened_rate_r_f17_pg_ne18`, already
zone-integrated to reactions/second via MESA's own `dm(k)` weighting) over
all zones, at each saved profile snapshot. Coarser in time than
[`production_rate`](@ref) -- limited to profile save points, not every
history row.
"""
function ne18_production_rate(run::MesaRun)
    pidx = profiles_index(run)
    hist = history(run).data
    age = Vector{Float64}(undef, nrow(pidx))
    rate = Vector{Float64}(undef, nrow(pidx))
    for (i, row) in enumerate(eachrow(pidx))
        prof = profile(run, row.profile_number).data
        rate[i] = sum(prof.screened_rate_r_f17_pg_ne18)
        j = findfirst(==(row.model_number), hist.model_number)
        age[i] = j === nothing ? NaN : hist.star_age[j] * YEAR_S
    end
    order = sortperm(age)
    return (age=age[order], rate=rate[order])
end

end # module MessengerProduction
