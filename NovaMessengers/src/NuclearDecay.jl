module NuclearDecay

using CSV
using DataFrames

export DecayIsotope, MESSENGER_ISOTOPES, EXTENDED_ISOTOPES, decay_constant, mass_number

"""
    DecayIsotope(name, daughter, halflife_s, q_neu_mev, mode,
                 positron_branching, gamma_line_mev, gamma_branching)

Weak-decay properties for one messenger isotope.

- `q_neu_mev`: average neutrino energy per decay (MeV), combined over
  whatever mix of decay branches the isotope actually has.
- `mode`: `:beta+` (pure positron emitter), `:electron_capture` (pure EC,
  no positron at all), or `:mixed` (both channels populated).
- `positron_branching`: fraction of decays that emit a positron --
  `1.0` for a pure beta+ emitter, `0.0` for pure EC, else the actual
  branching fraction. This is what [`MessengerProduction.positron_rate`](@ref)
  must scale by; treating every decay as "1 positron" (an idealization
  valid for the 7 MESA-tracked isotopes, whose EC side-branches are a
  few tenths of a percent or less) is silently wrong once EC is a
  double-digit-percent branch, as it is for na22/al26, or the whole
  story, as it is for be7.
- `gamma_line_mev` / `gamma_branching`: a characteristic *prompt daughter
  de-excitation* gamma line and its branching ratio per decay, when the
  daughter is left in an excited state (e.g. na22 -> ne22* -> ne22 + 1.275 MeV).
  This is a physically distinct emission mechanism from e+e- annihilation
  (511 keV) and gets its own fields rather than being folded into the
  511 keV machinery -- 0.0/0.0 for isotopes with no such line (all 7 in
  [`MESSENGER_ISOTOPES`](@ref) decay directly to their ground state in
  this simplified treatment).

No MESA-run dependency -- [`MESSENGER_ISOTOPES`](@ref) values are frozen
from `~/mesa-25.12.1/data/rates_data/weak_info.list` (the same table MESA
itself uses internally), so the physics here stays numerically consistent
with what produced the coupled hydro run. [`EXTENDED_ISOTOPES`](@ref)
values are NOT MESA-sourced (MESA's live net doesn't track these isotopes
at all) -- see that constant's own docstring for provenance and caveats.
"""
struct DecayIsotope
    name::Symbol
    daughter::Symbol
    halflife_s::Float64
    q_neu_mev::Float64
    mode::Symbol
    positron_branching::Float64
    gamma_line_mev::Float64
    gamma_branching::Float64
end

"""
    decay_constant(iso::DecayIsotope) -> Float64

Decay constant, in 1/second, from the isotope's half-life:

    lambda = ln(2) / halflife_s

the rate constant in the exponential decay law N(t) = N0 * exp(-lambda*t)
that [`MessengerProduction.decay_rate`](@ref) (`= lambda * N(t)`) is
built from.
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

function _load_decay_isotopes(csv_name::AbstractString)
    path = joinpath(@__DIR__, "..", "data", csv_name)
    df = CSV.read(path, DataFrame)
    pairs = [
        Symbol(row.isotope) => DecayIsotope(
            Symbol(row.isotope), Symbol(row.daughter), row.halflife_s, row.q_neu_mev, Symbol(row.mode),
            row.positron_branching, row.gamma_line_mev, row.gamma_branching,
        ) for row in eachrow(df)
    ]
    return Dict(pairs)
end

"""
    MESSENGER_ISOTOPES::Dict{Symbol,DecayIsotope}

The 7 beta+ emitters tracked by the CO WD baseline's live network
(n13, o14, o15, f17, f18, ne18, ne19), keyed by isotope symbol. Every
function defaulting to `isotopes=keys(MESSENGER_ISOTOPES)` (e.g.
`SignalSynthesis.neutrino_lightcurve`) relies on each of these having a
`total_mass_<iso>` column in MESA's own `history.data` -- deliberately
NOT merged with [`EXTENDED_ISOTOPES`](@ref), which don't.
"""
const MESSENGER_ISOTOPES = _load_decay_isotopes("weak_info_subset.csv")

"""
    EXTENDED_ISOTOPES::Dict{Symbol,DecayIsotope}

Isotopes beyond MESA's live network (na22, al26, be7), keyed by isotope
symbol. Their abundances come from `TrajectoryPostProcessing`'s ReacNetJl
bridge, not from MESA's `history.data`/`profile*.data` -- there is no
`total_mass_na22` column to read, so these are NOT in
[`MESSENGER_ISOTOPES`](@ref) and are not picked up by that dict's
default-argument callers. Use `ExtendedMessengers` for their production/
light-curve functions.

Provenance: standard nuclear data tables (ENSDF/NNDC-style
compilations), NOT frozen from MESA (its live net doesn't carry these
species at all, so there is no MESA-internal source to match). The
`q_neu_mev` values in particular are combined beta+/EC branch-weighted
averages assembled from textbook branching ratios and endpoint energies,
not independently re-derived here the way `MESSENGER_ISOTOPES`' Q-values
are cross-checked in `ReactionEnergetics.jl` -- treat these as good to
the quoted precision (2-3 significant figures), not exact.

`al26` here is the long-lived ground state (7.17e5 yr) only. Al-26 also
has a short-lived isomeric state (6.3 s, pure beta+ to Mg-26 ground, no
gamma line) that ReacNetJl's own output tracks separately (as `al*6`) --
deliberately not added here: it came out at ~1e-16 mass fraction in this
run's own end-to-end test (vs. ~1e-6 for the ground state), and its
`al*6`-style name doesn't fit `mass_number`'s parsing convention (see
that function) without a dedicated exception, so it's deferred rather
than force-fit for a channel that's negligible in this run anyway.
"""
const EXTENDED_ISOTOPES = _load_decay_isotopes("extended_isotopes.csv")

end # module NuclearDecay
