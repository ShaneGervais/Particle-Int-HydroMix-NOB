module ReactionEnergetics

using CSV
using DataFrames
using ..MesaIO

export REACTION_Q_VALUES_MEV, MEV_TO_ERG, TRACKED_REACTIONS,
       zone_reaction_energy_rate, reaction_energy_rate, zone_energy_breakdown

const MEV_TO_ERG = 1.602176634e-6  # erg/MeV (exact: 2019 SI redefinition fixes the elementary charge)

function _load_reaction_q_values()
    path = joinpath(@__DIR__, "..", "data", "reaction_q_values.csv")
    df = CSV.read(path, DataFrame)
    return Dict(Symbol(row.reaction) => row.q_value_mev for row in eachrow(df))
end

"""
    REACTION_Q_VALUES_MEV::Dict{Symbol,Float64}

Q-values (MeV) for the 18 reactions the CO WD baseline's live network has a
`screened_rate` profile column for (see `profile_columns.list`): the 7
weak decays in `MessengerProduction.DECAY_REACTIONS` and the 11
proton/alpha captures in `MessengerProduction.FORMATION_REACTIONS`.

Computed directly from MESA's own `data/chem_data/isotopes.data` mass
excesses (Q = sum(mass excess, reactants) - sum(mass excess, products),
minus 2*m_e*c^2 for the beta+ decays, using `Transport.ELECTRON_REST_ENERGY_MEV`),
NOT copied from an external compilation -- so these stay numerically
self-consistent with exactly the isotope masses this MESA version used
internally, the same reasoning as `NuclearDecay`'s `weak_info_subset.csv`.
Cross-checked: every decay Q here exceeds
`NuclearDecay.MESSENGER_ISOTOPES[iso].q_neu_mev` (the average neutrino
share alone), as it must.

One value is flagged rather than silently trusted: `r_o15_ap_f18`
(^15O(alpha,p)^18F) comes out endothermic (Q < 0) from this mass table.
That's not impossible (an endothermic channel can still run via a
resonance/thermal tail), but since this feeds a nuclear-rate sensitivity
study, treat it as "verify independently before relying on it," not a
settled number -- see the `note` column in `data/reaction_q_values.csv`.
"""
const REACTION_Q_VALUES_MEV = _load_reaction_q_values()

"""
    TRACKED_REACTIONS::Vector{Symbol}

All 18 reactions with a known Q-value here, sorted for reproducible
column ordering in [`zone_energy_breakdown`](@ref).
"""
const TRACKED_REACTIONS = sort(collect(keys(REACTION_Q_VALUES_MEV)))

"""
    zone_reaction_energy_rate(profile_data, reaction::Symbol) -> Vector{Float64}
    zone_reaction_energy_rate(run::MesaRun, profile_number::Integer, reaction::Symbol) -> Vector{Float64}

Per-zone energy generation rate (erg/second) from `reaction`: its
`screened_rate_<reaction>` profile column (reactions/second, already
integrated over the zone's mass -- MESA's own `dm(k)` weighting, same
convention `MessengerProduction` reads) times its Q-value.
"""
function zone_reaction_energy_rate(profile_data::AbstractDataFrame, reaction::Symbol)
    q = get(REACTION_Q_VALUES_MEV, reaction) do
        error("no known Q-value for reaction $reaction")
    end
    col = Symbol("screened_rate_", reaction)
    col in propertynames(profile_data) || error("profile data has no column $col")
    return profile_data[!, col] .* (q * MEV_TO_ERG)
end

function zone_reaction_energy_rate(run::MesaRun, profile_number::Integer, reaction::Symbol)
    return zone_reaction_energy_rate(profile(run, profile_number).data, reaction)
end

"""
    reaction_energy_rate(run::MesaRun, reaction::Symbol) -> (age, rate)

Whole-star energy generation rate (erg/second) from `reaction`, summed
over all zones at each saved profile snapshot -- profile cadence, not
every history row (same tradeoff as `MessengerProduction.reaction_decay_rate`,
which this mirrors).
"""
function reaction_energy_rate(run::MesaRun, reaction::Symbol)
    pidx = profiles_index(run)
    hist = history(run).data
    age = Vector{Float64}(undef, nrow(pidx))
    rate = Vector{Float64}(undef, nrow(pidx))
    for (i, row) in enumerate(eachrow(pidx))
        rate[i] = sum(zone_reaction_energy_rate(run, row.profile_number, reaction))
        j = findfirst(==(row.model_number), hist.model_number)
        age[i] = j === nothing ? NaN : hist.star_age[j] * (365.25 * 24 * 3600)
    end
    order = sortperm(age)
    return (age=age[order], rate=rate[order])
end

"""
    zone_energy_breakdown(profile_data::AbstractDataFrame) -> DataFrame

Per-zone energy generation rate (erg/second) from every tracked reaction,
one column per reaction (e.g. `:r_c12_pg_n13`, in `TRACKED_REACTIONS`
order), plus a `:total` column summing them. This is the reaction-resolved
decomposition of (a slice of) what MESA's own `eps_nuc` profile column
already gives you in bulk -- useful for seeing which specific reaction
dominates energy release at a given depth/time, and for comparing that
`:total` against `eps_nuc` as a sanity check (it will be a fraction of
`eps_nuc`, not all of it: `eps_nuc` includes every reaction in the live
network, e.g. the pp chain and triple-alpha, not just these 18).
"""
function zone_energy_breakdown(profile_data::AbstractDataFrame)
    df = DataFrame([r => zone_reaction_energy_rate(profile_data, r) for r in TRACKED_REACTIONS])
    df.total = sum(eachcol(df))
    return df
end

end # module ReactionEnergetics
