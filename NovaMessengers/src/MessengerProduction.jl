module MessengerProduction

using DataFrames
using ..MesaIO
using ..NuclearDecay

export isotope_number, decay_rate, positron_rate, neutrino_energy_loss_rate,
       annihilation_photon_rate, DECAY_REACTIONS, zone_decay_rate,
       zone_annihilation_photon_rate, reaction_decay_rate,
       FORMATION_REACTIONS, zone_formation_rate, formation_rate

const AMU_G = 1.66053906892e-24  # grams (CODATA)
const YEAR_S = 365.25 * 24 * 3600

"""
    DECAY_REACTIONS::Dict{Symbol,Symbol}

The single weak-decay reaction (as named in the CO WD baseline's
`cno_extras_o18_to_mg26_plus_fe56.net`) responsible for each messenger
isotope's beta+ decay -- one positron + one neutrino per reaction. This
*is* the messenger emission channel; matches `NuclearDecay.MESSENGER_ISOTOPES`'
`daughter` field exactly.
"""
const DECAY_REACTIONS = Dict{Symbol,Symbol}(
    :n13 => :r_n13_wk_c13,
    :o14 => :r_o14_wk_n14,
    :o15 => :r_o15_wk_n15,
    :f17 => :r_f17_wk_o17,
    :f18 => :r_f18_wk_o18,
    :ne18 => :r_ne18_wk_f18,
    :ne19 => :r_ne19_wk_f19,
)

"""
    FORMATION_REACTIONS::Dict{Symbol,Vector{Symbol}}

The proton/alpha-capture reaction(s) that *form* each messenger isotope
(some have more than one channel, e.g. f18 has three). This is a genuinely
different quantity from [`DECAY_REACTIONS`](@ref): during hot-CNO breakout
burning, proton captures onto an isotope can run faster than its own
decay, consuming it before it gets a chance to emit a messenger particle.
Useful for nucleosynthesis-flow bookkeeping (how fast is this isotope
being synthesized), but NOT the messenger emission rate -- use
[`zone_decay_rate`](@ref) / [`decay_rate`](@ref) for that.
"""
const FORMATION_REACTIONS = Dict{Symbol,Vector{Symbol}}(
    :n13 => [:r_c12_pg_n13],
    :o14 => [:r_n13_pg_o14],
    :o15 => [:r_n14_pg_o15],
    :f17 => [:r_o16_pg_f17, :r_o14_ap_f17],
    :f18 => [:r_o17_pg_f18, :r_n14_ag_f18, :r_o15_ap_f18],
    :ne18 => [:r_f17_pg_ne18],
    :ne19 => [:r_o15_ag_ne19, :r_f18_pg_ne19],
)

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
    decay_rate(run::MesaRun, isotope::Symbol) -> (age, rate)

Whole-star messenger emission rate (decays/second, = positrons/second =
neutrinos/second) of `isotope`, at history.data's full time cadence:
`lambda * N(t)`, directly -- no finite-differencing needed, since the
decay rate is exactly the decay constant times the current inventory by
definition. (An earlier version of this function tried to reconstruct a
"production rate" by finite-differencing `total_mass` and correcting for
decay alone; that conflated isotope *formation* with isotope *decay* and
silently ignored further consumption via proton capture, which during
hot-CNO burning can dominate over decay -- see [`FORMATION_REACTIONS`](@ref)
for why those are a different quantity.)
"""
function decay_rate(run::MesaRun, isotope::Symbol)
    N = isotope_number(run, isotope)
    age = history(run).data.star_age .* YEAR_S
    lambda = decay_constant(MESSENGER_ISOTOPES[isotope])
    return (age=age, rate=lambda .* N)
end

"""
    positron_rate(run::MesaRun, isotope::Symbol) -> (age, rate)

Positron production rate (positrons/second) -- one per beta+ decay, so
identical to [`decay_rate`](@ref) for these isotopes (all beta+, see
`NuclearDecay.MESSENGER_ISOTOPES`).
"""
positron_rate(run::MesaRun, isotope::Symbol) = decay_rate(run, isotope)

"""
    neutrino_energy_loss_rate(run::MesaRun, isotope::Symbol) -> (age, rate)

Neutrino energy loss rate (MeV/second) from `isotope`'s beta+ decay:
decay rate times the average neutrino energy per decay (`q_neu_mev`).
"""
function neutrino_energy_loss_rate(run::MesaRun, isotope::Symbol)
    r = decay_rate(run, isotope)
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
    _sum_zone_reaction_rates(profile_data, reaction_names) -> Vector{Float64}

Sum the `screened_rate_<reaction>` profile columns for `reaction_names`,
zone by zone. Each column is already integrated over that zone's mass
(MESA's own `dm(k)` weighting), so no further normalization is needed.
Shared numeric core of [`zone_decay_rate`](@ref) and
[`zone_formation_rate`](@ref).
"""
function _sum_zone_reaction_rates(profile_data::AbstractDataFrame, reaction_names)
    total = zeros(Float64, nrow(profile_data))
    for r in reaction_names
        col = Symbol("screened_rate_", r)
        col in propertynames(profile_data) || error("profile data has no column $col")
        total .+= profile_data[!, col]
    end
    return total
end

"""
    zone_decay_rate(profile_data::AbstractDataFrame, isotope::Symbol) -> Vector{Float64}
    zone_decay_rate(run::MesaRun, profile_number::Integer, isotope::Symbol) -> Vector{Float64}

Per-zone messenger emission rate (decays/second) of `isotope`, from its
single weak-decay channel ([`DECAY_REACTIONS`](@ref)). This is the
per-zone quantity meant to be multiplied against Phase 3's
`escape_probability_gamma` (also per-zone, same profile) before summing
over zones -- unlike [`decay_rate`](@ref), which only has a whole-star
total available and can't be weighted by depth. The `MesaRun` form reads
the requested profile; the `AbstractDataFrame` form takes an
already-loaded one, so callers combining several isotopes / a Transport
escape probability from the same snapshot (e.g. `SignalSynthesis`) don't
re-read the same multi-MB profile file over and over.
"""
function zone_decay_rate(profile_data::AbstractDataFrame, isotope::Symbol)
    reaction = get(DECAY_REACTIONS, isotope) do
        error("no known decay reaction for isotope $isotope")
    end
    return _sum_zone_reaction_rates(profile_data, (reaction,))
end

function zone_decay_rate(run::MesaRun, profile_number::Integer, isotope::Symbol)
    return zone_decay_rate(profile(run, profile_number).data, isotope)
end

"""
    zone_annihilation_photon_rate(profile_data::AbstractDataFrame, isotope::Symbol) -> Vector{Float64}
    zone_annihilation_photon_rate(run::MesaRun, profile_number::Integer, isotope::Symbol) -> Vector{Float64}

Per-zone 511 keV annihilation photon production rate (photons/second):
two photons per positron, assuming in-situ annihilation (site-of-annihilation
/ true in-flight annihilation is a further transport refinement, deferred).
"""
function zone_annihilation_photon_rate(profile_data::AbstractDataFrame, isotope::Symbol)
    return 2 .* zone_decay_rate(profile_data, isotope)
end

function zone_annihilation_photon_rate(run::MesaRun, profile_number::Integer, isotope::Symbol)
    return zone_annihilation_photon_rate(profile(run, profile_number).data, isotope)
end

"""
    reaction_decay_rate(run::MesaRun, isotope::Symbol) -> (age, rate)

Whole-star messenger emission rate (decays/second) of `isotope`, from
summing [`zone_decay_rate`](@ref) over all zones at each saved profile
snapshot. Coarser in time than [`decay_rate`](@ref) (limited to profile
save points, not every history row), but computed from an entirely
independent path -- MESA's local, screening-corrected reaction rate,
rather than whole-star mass bookkeeping -- so the two are a genuine
cross-check on each other.
"""
function reaction_decay_rate(run::MesaRun, isotope::Symbol)
    pidx = profiles_index(run)
    hist = history(run).data
    age = Vector{Float64}(undef, nrow(pidx))
    rate = Vector{Float64}(undef, nrow(pidx))
    for (i, row) in enumerate(eachrow(pidx))
        rate[i] = sum(zone_decay_rate(run, row.profile_number, isotope))
        j = findfirst(==(row.model_number), hist.model_number)
        age[i] = j === nothing ? NaN : hist.star_age[j] * YEAR_S
    end
    order = sortperm(age)
    return (age=age[order], rate=rate[order])
end

"""
    zone_formation_rate(profile_data::AbstractDataFrame, isotope::Symbol) -> Vector{Float64}
    zone_formation_rate(run::MesaRun, profile_number::Integer, isotope::Symbol) -> Vector{Float64}

Per-zone formation (synthesis) rate (nuclei/second) of `isotope`, summed
over all of its production channels ([`FORMATION_REACTIONS`](@ref)). NOT
the messenger emission rate -- see [`zone_decay_rate`](@ref) for that.
"""
function zone_formation_rate(profile_data::AbstractDataFrame, isotope::Symbol)
    reactions = get(FORMATION_REACTIONS, isotope) do
        error("no known formation reaction(s) for isotope $isotope")
    end
    return _sum_zone_reaction_rates(profile_data, reactions)
end

function zone_formation_rate(run::MesaRun, profile_number::Integer, isotope::Symbol)
    return zone_formation_rate(profile(run, profile_number).data, isotope)
end

"""
    formation_rate(run::MesaRun, isotope::Symbol) -> (age, rate)

Whole-star formation (synthesis) rate (nuclei/second) of `isotope`, from
summing [`zone_formation_rate`](@ref) over all zones at each saved profile
snapshot. A nucleosynthesis-flow diagnostic, not the messenger emission
rate: comparing this to [`reaction_decay_rate`](@ref) shows how much of a
freshly-formed isotope gets consumed by further captures before it has a
chance to decay (the hallmark of hot-CNO breakout burning).
"""
function formation_rate(run::MesaRun, isotope::Symbol)
    pidx = profiles_index(run)
    hist = history(run).data
    age = Vector{Float64}(undef, nrow(pidx))
    rate = Vector{Float64}(undef, nrow(pidx))
    for (i, row) in enumerate(eachrow(pidx))
        rate[i] = sum(zone_formation_rate(run, row.profile_number, isotope))
        j = findfirst(==(row.model_number), hist.model_number)
        age[i] = j === nothing ? NaN : hist.star_age[j] * YEAR_S
    end
    order = sortperm(age)
    return (age=age[order], rate=rate[order])
end

end # module MessengerProduction
