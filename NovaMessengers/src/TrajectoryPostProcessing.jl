module TrajectoryPostProcessing

using Printf
using DataFrames
using ReacNetJl
using ..MesaIO
using ..NuclearDecay: mass_number

export LIVE_NET_ISOTOPES, peak_temperature_zone, zone_trajectory,
       write_trajectory_file, zone_initial_abundances,
       write_initial_abundance_file, pristine_profile_number, postprocess_trajectory

"""
    LIVE_NET_ISOTOPES::Vector{Symbol}

Every isotope in `cno_extras_o18_to_mg26_plus_fe56.net` (the CO WD
baseline's live network: `basic.net` + `add_hot_cno` + `add_cno_extras` +
`add_o18_and_ne22` + explicit `mg26`/`fe56`) that has its own tracked
abundance -- i.e. the species `add_abundances` in `profile_columns.list`
exposes as bare-name profile columns (`h1`, `c12`, ... not
`x_mass_fraction_h1`). This is the complete initial-composition vector
available to hand to ReacNetJl; anything beyond it (na23, al27, mg25,
22Na, 26Al, 7Be, ...) has no MESA-tracked abundance in this net and is
left absent from the written abundance file, which ReacNetJl treats as
zero -- physically correct for isotopes this CO WD's dredge-up doesn't
supply in the first place, not a placeholder to fill in later.
"""
const LIVE_NET_ISOTOPES = Symbol[
    :h1, :he3, :he4, :c12, :c13, :n13, :n14, :n15, :o14, :o15, :o16, :o17,
    :o18, :f17, :f18, :f19, :ne18, :ne19, :ne20, :ne22, :mg22, :mg24, :mg26, :fe56,
]

"""
ReacNetJl's STARLIB/Iliadis2002 rate tables are tabulated over
T9 in [1e-3, 10] (1e6-1e10 K) -- no rate library tabulates below ~1e6 K
because nuclear reaction rates there are utterly negligible regardless
of the exact temperature, so there is nothing physically lost by
flooring a trajectory's sub-table epochs to (just above) the table's
own minimum (see [`write_trajectory_file`](@ref)). `REACNETJL_T9_MIN`
sits 1% above the table's literal floor of 1e-3, not exactly on it:
the solver interpolates T9 *between* saved trajectory nodes, and a node
sitting exactly on the boundary can produce an interpolated value a
few ULPs below it (confirmed: `T9=0.0009999999999999998` from two nodes
both >= 0.001) -- physically inert either way, so the margin costs
nothing.
"""
const REACNETJL_T9_MIN = 1.01e-3
const REACNETJL_T9_MAX = 10.0

const ELEMENT_Z = Dict(
    "h" => 1, "he" => 2, "c" => 6, "n" => 7, "o" => 8, "f" => 9,
    "ne" => 10, "na" => 11, "mg" => 12, "al" => 13, "fe" => 26,
)

const YEAR_S = 365.25 * 24 * 3600

"""
Linear interpolation of `ys` at `x`, given `xs` strictly increasing.
Clamps to the nearest endpoint outside `[xs[1], xs[end]]` rather than
erroring -- the tracked mass coordinate is fixed across profiles whose
mass grids shift slightly (accretion, remeshing), so a few-ULP excursion
at either edge is expected, not a bug.
"""
function _linterp(xs::AbstractVector, ys::AbstractVector, x::Real)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    k = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    t = (x - xs[k]) / (xs[k + 1] - xs[k])
    return ys[k] + t * (ys[k + 1] - ys[k])
end

"""
    peak_temperature_zone(run::MesaRun) -> (profile_number, mass_coordinate, temperature_K)

The mass coordinate (Msun) that reaches the highest temperature at any
saved profile snapshot in the run -- the standard single-zone choice for
nova nucleosynthesis post-processing (the fluid parcel best representing
peak hot-CNO burning conditions). Returned alongside the profile number
where that peak occurs and the peak temperature itself (from the
profile's own `temperature` column, K, not `logT`, avoiding a log/delog
round trip).
"""
function peak_temperature_zone(run::MesaRun)
    pidx = profiles_index(run)
    best_profile = -1
    best_mass = NaN
    best_T = -Inf
    for row in eachrow(pidx)
        prof = profile(run, row.profile_number).data
        k = argmax(prof.temperature)
        if prof.temperature[k] > best_T
            best_profile = row.profile_number
            best_mass = prof.mass[k]
            best_T = prof.temperature[k]
        end
    end
    best_profile == -1 && error("run has no profiles")
    return (profile_number=best_profile, mass_coordinate=best_mass, temperature_K=best_T)
end

"""
    zone_trajectory(run::MesaRun, mass_coordinate::Real) -> (time_s, T9, rho)

Single-zone T9(t)/rho(t) trajectory at a fixed mass coordinate (Msun),
built by linearly interpolating each saved profile's `mass`/`temperature`/
`logRho` columns at that coordinate. Interpolates on mass rather than
following a raw zone index: profile zones are numbered surface-to-center
and MESA's zone count/spacing changes between saved snapshots (adaptive
remeshing), so a fixed zone *index* does not track the same fluid parcel
over time the way a fixed mass coordinate does. `time_s` comes from
`history.data`'s `star_age` (years -> seconds), sorted strictly
increasing -- ReacNetJl's own trajectory reader requires that.
"""
function zone_trajectory(run::MesaRun, mass_coordinate::Real)
    pidx = profiles_index(run)
    hist = history(run).data
    n = nrow(pidx)
    time_s = Vector{Float64}(undef, n)
    T9 = Vector{Float64}(undef, n)
    rho = Vector{Float64}(undef, n)
    for (i, row) in enumerate(eachrow(pidx))
        prof = profile(run, row.profile_number).data
        order = sortperm(prof.mass)  # zones run surface->center: mass decreasing, reverse for interpolation
        m = prof.mass[order]
        T9[i] = _linterp(m, prof.temperature[order], mass_coordinate) / 1.0e9
        rho[i] = 10.0^_linterp(m, prof.logRho[order], mass_coordinate)
        j = findfirst(==(row.model_number), hist.model_number)
        time_s[i] = j === nothing ? NaN : hist.star_age[j] * YEAR_S
    end
    order = sortperm(time_s)
    return (time_s=time_s[order], T9=T9[order], rho=rho[order])
end

"""
    write_trajectory_file(path, time_s, T9, rho) -> path

Write a ReacNetJl-format `trajectory.input`: `AGEUNIT=SEC`/`TUNIT=T9K`/
`RHOUNIT=CGS` header (this project's own natural units, so no further
conversion happens at read time), then `time  T9  rho` rows, at full
`%.17e` precision (the number of significant decimal digits needed to
round-trip any Float64 exactly). This run's saved profiles span from
~1e2 s to ~1e11-1e12 s (seconds-scale burst dynamics riding on top of a
multi-millennium pre-outburst accretion history in the same trajectory),
so anything less than full precision silently collapses distinct,
closely-spaced burst-era timestamps to identical printed values once
rounded against that much larger baseline -- which then fails
ReacNetJl's strict-monotonicity check on read-back, not obviously
connected to a formatting choice from the error alone.

`T9` is clamped to [`REACNETJL_T9_MIN`, `REACNETJL_T9_MAX`] before
writing: this run's early, pre-ignition epochs sit as low as T9~8e-5
(~79,000 K) at the tracked zone, below every rate table's tabulated
floor (~1e6 K) -- physically inert (nuclear rates there are ~0
regardless of the exact sub-floor value), so flooring changes nothing
about the burn but avoids an otherwise-opaque out-of-table-range error
deep in the solver.
"""
function write_trajectory_file(path::AbstractString, time_s::AbstractVector{<:Real},
    T9::AbstractVector{<:Real}, rho::AbstractVector{<:Real})
    n = length(time_s)
    (n == length(T9) == length(rho)) || error("time_s, T9, rho must be the same length")
    all(diff(time_s) .> 0) || error("time_s must be strictly increasing")
    T9_clamped = clamp.(T9, REACNETJL_T9_MIN, REACNETJL_T9_MAX)
    open(path, "w") do io
        println(io, "AGEUNIT = SEC")
        println(io, "TUNIT   = T9K")
        println(io, "RHOUNIT = CGS")
        for i in 1:n
            @printf(io, "%.17e  %.17e  %.17e\n", time_s[i], T9_clamped[i], rho[i])
        end
    end
    return path
end

"""
    zone_initial_abundances(run::MesaRun, profile_number::Integer, mass_coordinate::Real) -> Dict{String,Float64}

Mass-fraction composition at `mass_coordinate`, interpolated the same way
[`zone_trajectory`](@ref) does, from the [`LIVE_NET_ISOTOPES`](@ref)
abundance columns of the given profile -- pass the *earliest* profile in
the run (pre-burst, pristine composition) to get the actual starting
point of a trajectory, not the profile at peak temperature. Keyed by
isotope name in the STARLIB-style convention ReacNetJl's
`normalize_species_name` expects, which is already this project's own
isotope-Symbol convention (`String(iso)`, no translation needed).

Requires `add_abundances` to be enabled in `profile_columns.list` (see
that file's own comment) and the run to have been regenerated after that
edit -- errors clearly if the columns aren't present rather than
silently returning an empty/wrong composition.
"""
function zone_initial_abundances(run::MesaRun, profile_number::Integer, mass_coordinate::Real)
    prof = profile(run, profile_number).data
    cols = propertynames(prof)
    isos = filter(iso -> iso in cols, LIVE_NET_ISOTOPES)
    isempty(isos) && error(
        "no live-net isotope abundance columns found in profile $profile_number -- " *
        "profile_columns.list needs `add_abundances` enabled and this run rerun " *
        "(see mesa_work/wd_nova_burst_co/profile_columns.list)",
    )
    order = sortperm(prof.mass)
    m = prof.mass[order]
    return Dict(String(iso) => _linterp(m, prof[!, iso][order], mass_coordinate) for iso in isos)
end

function _element_symbol(iso::AbstractString)
    m = match(r"^([a-zA-Z]+)(\d+)$", iso)
    m === nothing && error("could not parse element symbol from isotope name \"$iso\"")
    return lowercase(m.captures[1])
end

"""
    write_initial_abundance_file(path, abundances::Dict{String,Float64}) -> path

Write a ReacNetJl-format `initial_abundance.dat`: one `Z  sym A  X` row
per isotope (mass fraction), matching `nova_cases/*/initial_abundance.dat`'s
own layout. `ELEMENT_Z` only needs to cover elements that actually appear
in [`LIVE_NET_ISOTOPES`](@ref) -- extend it if the live net ever grows.
"""
function write_initial_abundance_file(path::AbstractString, abundances::Dict{String,Float64})
    open(path, "w") do io
        for (iso, X) in sort(collect(abundances))
            sym = _element_symbol(iso)
            A = mass_number(Symbol(iso))
            Z = get(ELEMENT_Z, sym) do
                error("no known atomic number for element \"$sym\" (isotope $iso) -- add it to ELEMENT_Z")
            end
            @printf(io, "%3d %-4s%3d  %.10E\n", Z, sym, A, X)
        end
    end
    return path
end

"""
    pristine_profile_number(run::MesaRun, mass_coordinate::Real) -> Integer

The profile at which `mass_coordinate`'s own hydrogen mass fraction (`h1`,
interpolated the same way [`zone_trajectory`](@ref) does) is at its
*maximum* across the whole run -- NOT simply the earliest saved profile.

This matters because of a real finding from this run's own data: at the
CO WD baseline's peak-temperature zone, the *first* saved profile has
h1 ~ 6e-4 (99.9% he4) while the burst-peak profile has h1 ~ 0.68 -- i.e.
hydrogen at that fixed mass coordinate *increases* between the first
snapshot and the burst peak, which pure nuclear burning cannot do.
That's convective mixing: the TNR's convective zone grows and dredges
H-rich material down to mass coordinates that were H-poor before mixing
reached them. Using the literal first profile as "the initial
composition" would hand ReacNetJl an artifact of pre-mixing local
bookkeeping, not the fuel that's actually about to burn. Taking the
run's own H-maximum instead is a standard proxy in single-zone nova
post-processing for "composition right as convective mixing completes
and burning begins" -- H can only decrease from that point by nuclear
consumption (this single-zone approximation does not itself model
further mixing pulses after that point; MESA's live net already
captures those self-consistently in its own `total_mass_<iso>` history,
which remains the cross-check).

Among profiles at (or within `rtol` of) that maximum, the *latest* one
chronologically is picked, not the first. Confirmed from this run's own
data: h1 at the tracked zone hits 0.70 at profile 13 (~25 hours into a
run that starts from an artificial relaxed initial model) and then sits
flat at 0.70 for the next ~22,000 years (profiles 13-41) before
compressional heating starts consuming it -- picking the *first*
occurrence (profile 13) is still a physically valid pristine
composition, but leaves the resulting trajectory needing to traverse
~10,000+ physically near-inert years for no benefit even with a tight
`rtol` (confirmed: both the literal first occurrence and a 0.1%
tolerance still exhaust ReacNetJl's `max_steps` in the millions before
converging).

The default `rtol=5e-3` (0.5%) was chosen empirically from this run: h1
does not decline monotonically once heating starts -- it visibly
oscillates (0.694, 0.695, 0.694, 0.697, ...) over the final ~1.5 years
before peak, consistent with intermittent convective overturn pulses
still bringing some fresh H down to this coordinate even as burning
consumes it. `findlast` at 0.5% lands on the *last* such upward wobble,
which happens to sit right at the edge of the explosive rise -- cutting
the remaining trajectory span from ~10,000 years to ~1.5 years (about
4 orders of magnitude), the difference between the solver not
converging and converging.
"""
function pristine_profile_number(run::MesaRun, mass_coordinate::Real; rtol::Real=5.0e-3)
    pidx = profiles_index(run)
    sort!(pidx, :model_number)
    profile_numbers = Int[]
    h1_values = Float64[]
    for row in eachrow(pidx)
        prof = profile(run, row.profile_number).data
        :h1 in propertynames(prof) || error(
            "profile $(row.profile_number) has no h1 column -- needs `add_abundances` " *
            "enabled in profile_columns.list and a rerun",
        )
        order = sortperm(prof.mass)
        push!(profile_numbers, row.profile_number)
        push!(h1_values, _linterp(prof.mass[order], prof.h1[order], mass_coordinate))
    end
    isempty(profile_numbers) && error("run has no profiles")
    max_h1 = maximum(h1_values)
    threshold = (1 - rtol) * max_h1
    k = findlast(>=(threshold), h1_values)
    best_profile = profile_numbers[k]
    return best_profile
end

"""
    postprocess_trajectory(run::MesaRun; output_dir, rates=:starlib, kwargs...) -> NamedTuple

End-to-end bridge from a MESA run to a ReacNetJl post-processing result:
1. Locate the peak-temperature zone ([`peak_temperature_zone`](@ref)).
2. Find its composition at H-maximum ([`zone_initial_abundances`](@ref)
   at [`pristine_profile_number`](@ref) -- see that function's docstring
   for why this is not simply the earliest profile).
3. Extract the T9(t)/rho(t) trajectory ([`zone_trajectory`](@ref))
   **from that pristine epoch onward only**, rebased so it starts at
   t=0 there. This run's saved history spans ~1e2 s to ~1e11 s (a
   multi-millennium pre-outburst accretion tail before the tracked zone
   is even done mixing), and that whole leading span is both physically
   outside what "evolve this starting composition forward" means and
   expensive for the solver to step through for no reason (confirmed:
   including it hits ReacNetJl's default `max_steps=1_000_000` without
   converging) -- trimming to the pristine epoch onward fixes both.
4. Write both in ReacNetJl's own file formats under `output_dir`.
5. Call `ReacNetJl.run_ppn` directly, in-process -- no CSV round trip on
   the way back; `result.final_mass_fractions` is a `Dict{String,Float64}`
   already in the same isotope-name convention this project uses.

`kwargs` are forwarded to `run_ppn` (e.g. `screening`, `neutron_captures`,
`max_steps`). Returns `(mass_coordinate, trajectory, initial_abundances,
result)` so callers know which zone/time window the returned abundances
describe, not just the final numbers.
"""
function postprocess_trajectory(run::MesaRun; output_dir::AbstractString, rates::Symbol=:starlib,
    pristine_rtol::Real=5.0e-3, kwargs...)
    peak = peak_temperature_zone(run)
    pn = pristine_profile_number(run, peak.mass_coordinate; rtol=pristine_rtol)
    abund = zone_initial_abundances(run, pn, peak.mass_coordinate)

    hist = history(run).data
    pidx = profiles_index(run)
    pristine_model_number = pidx[findfirst(==(pn), pidx.profile_number), :model_number]
    j = findfirst(==(pristine_model_number), hist.model_number)
    t0 = hist.star_age[j] * YEAR_S

    full_traj = zone_trajectory(run, peak.mass_coordinate)
    keep = full_traj.time_s .>= t0
    traj = (time_s=full_traj.time_s[keep] .- t0, T9=full_traj.T9[keep], rho=full_traj.rho[keep])

    mkpath(output_dir)
    traj_path = joinpath(output_dir, "trajectory.input")
    abund_path = joinpath(output_dir, "initial_abundance.dat")
    write_trajectory_file(traj_path, traj.time_s, traj.T9, traj.rho)
    write_initial_abundance_file(abund_path, abund)

    # run_ppn's own dt_max default is duration>100s ? 20.0 : 0.05 (seconds) --
    # a heuristic clearly tuned for hours-to-days-scale nova trajectories.
    # Ours spans up to ~years (this run's tracked zone has a slow pre-peak
    # buildup on top of the fast burst itself), so that 20s ceiling forces
    # well over 2e6 steps just to traverse the quiet stretches at maximum
    # stride -- already past max_steps before any physics-driven shrinking
    # during the actual burst. Scale dt_max to this trajectory's own
    # duration instead, letting max_fractional_change (not this ceiling)
    # be what controls resolution during the fast phase. An explicit
    # `dt_max` in kwargs still wins (merge order below).
    duration = traj.time_s[end] - traj.time_s[1]
    default_dt_max = max(20.0, duration / 500)
    run_ppn_kwargs = merge((dt_max=default_dt_max,), NamedTuple(kwargs))

    result = ReacNetJl.run_ppn(traj_path, abund_path; rates=rates, output_dir=output_dir, run_ppn_kwargs...)
    return (mass_coordinate=peak.mass_coordinate, trajectory=traj, initial_abundances=abund, result=result)
end

end # module TrajectoryPostProcessing
