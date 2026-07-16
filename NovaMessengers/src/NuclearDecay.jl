module NuclearDecay

using CSV
using DataFrames

export DecayIsotope, MESSENGER_ISOTOPES, decay_constant, mass_number

"""
    DecayIsotope(name, daughter, halflife_s, q_neu_mev, mode)

Weak-decay properties for one messenger isotope. `q_neu_mev` is the average
neutrino energy per decay (MeV). No MESA-run dependency -- values are frozen
from `~/mesa-25.12.1/data/rates_data/weak_info.list` (the same table MESA
itself uses internally), so the physics here stays numerically consistent
with what produced the coupled hydro run.
"""
struct DecayIsotope
    name::Symbol
    daughter::Symbol
    halflife_s::Float64
    q_neu_mev::Float64
    mode::Symbol
end

"""
    decay_constant(iso::DecayIsotope) -> Float64

Decay constant lambda = ln(2) / halflife, in 1/second.
"""
decay_constant(iso::DecayIsotope) = log(2) / iso.halflife_s

"""
    mass_number(isotope::Symbol) -> Int

Mass number A parsed from a MESA-style isotope name (e.g. `:f18` -> 18).
"""
function mass_number(isotope::Symbol)
    s = String(isotope)
    digits_str = match(r"\d+$", s)
    digits_str === nothing && error("could not parse mass number from isotope name $isotope")
    return parse(Int, digits_str.match)
end

function _load_weak_info_subset()
    path = joinpath(@__DIR__, "..", "data", "weak_info_subset.csv")
    df = CSV.read(path, DataFrame)
    pairs = [
        Symbol(row.isotope) => DecayIsotope(
            Symbol(row.isotope), Symbol(row.daughter), row.halflife_s, row.q_neu_mev, Symbol(row.mode),
        ) for row in eachrow(df)
    ]
    return Dict(pairs)
end

"""
    MESSENGER_ISOTOPES::Dict{Symbol,DecayIsotope}

The 7 beta+ emitters tracked by the CO WD baseline's live network
(n13, o14, o15, f17, f18, ne18, ne19), keyed by isotope symbol.
"""
const MESSENGER_ISOTOPES = _load_weak_info_subset()

end # module NuclearDecay
