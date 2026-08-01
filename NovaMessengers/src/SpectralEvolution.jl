module SpectralEvolution

using ..TrajectoryPostProcessing: _linterp

export LineChannel, ContinuumChannel, SpectrumSnapshot, composite_spectrum, spectrum_timeline

const MEV_TO_EV = 1.0e6
const GEV_TO_EV = 1.0e9

"""
    LineChannel(label, energy_ev, age_s, rate)

A discrete-energy messenger channel -- a neutrino mean-energy line
(from [`NuclearDecay.DecayIsotope`](@ref)'s `q_neu_mev`) or a gamma-ray
decay line (511 keV / 1.275 MeV / 1.809 MeV / 478 keV) -- carrying its
own precomputed light curve. `age_s`/`rate` come from whichever
producing function is appropriate (e.g. [`MessengerProduction.decay_rate`](@ref),
[`SignalSynthesis.gamma_lightcurve`](@ref),
[`ExtendedMessengers.freezeout_gamma_lightcurve`](@ref)) -- built once
(potentially expensive, since some of those trace back to a ReacNetJl
solve), then interpolated cheaply at however many times a movie needs
via [`composite_spectrum`](@ref).
"""
struct LineChannel
    label::String
    energy_ev::Float64
    age_s::Vector{Float64}
    rate::Vector{Float64}
end

"""
    ContinuumChannel(label, luminosity_ev)

A continuum messenger channel: `luminosity_ev(age_s, E_ev) -> erg/s/eV`,
a closure wrapping whichever closed-form physics is appropriate (e.g.
`(t, E) -> QuiescentContinuum.spectral_luminosity_ev(source, E)` for the
WD/companion baseline, `(t, E) -> ShockAcceleration.shock_bremsstrahlung_luminosity_ev(p, t, E)`
for shock X-rays). Unlike [`LineChannel`](@ref), these are cheap
closed-form evaluations -- no precomputed table needed, just called
directly at whatever `(age_s, E_ev)` a snapshot asks for.
"""
struct ContinuumChannel
    label::String
    luminosity_ev::Function
end

"""
    SpectrumSnapshot(age_s, line_labels, line_energies_ev, line_rates,
                      energy_grid_ev, continuum_labels, continuum_L_ev)

The full multi-messenger spectrum at one instant: discrete lines
(`line_energies_ev`/`line_rates`, particle or photon rate at each line's
energy) plus continua (`continuum_L_ev[i, j]` = channel `i`'s spectral
luminosity at `energy_grid_ev[j]`). Lines and continua are kept as
separate arrays rather than merged onto one grid because they're
different physical quantities (a rate at a single energy vs. a spectral
density) and, for the neutrino lines specifically, a different *particle
type* than every other channel here -- summing them together would be
physically meaningless, even though plotting them on the same log-log
energy axis (see `examples/spectrum_movie.jl`) is standard practice in
real multi-messenger spectra.
"""
struct SpectrumSnapshot
    age_s::Float64
    line_labels::Vector{String}
    line_energies_ev::Vector{Float64}
    line_rates::Vector{Float64}
    energy_grid_ev::Vector{Float64}
    continuum_labels::Vector{String}
    continuum_L_ev::Matrix{Float64}
end

function _line_rate_at(lc::LineChannel, age_s::Real)
    isempty(lc.age_s) && return 0.0
    return _linterp(lc.age_s, lc.rate, age_s)
end

"""
    composite_spectrum(age_s, lines, continua, energy_grid_ev) -> SpectrumSnapshot

Assemble one [`SpectrumSnapshot`](@ref) at `age_s` (seconds, `star_age`
convention) from a set of [`LineChannel`](@ref)s (interpolated at
`age_s`:

    rate_i(age_s) = _linterp(line_i.age_s, line_i.rate, age_s)

) and [`ContinuumChannel`](@ref)s (evaluated directly:

    continuum_L_ev[i, j] = continuum_i.luminosity_ev(age_s, energy_grid_ev[j])

). Cheap enough to call once per movie frame even though the underlying
line light curves may themselves have been expensive to produce.
"""
function composite_spectrum(age_s::Real, lines::AbstractVector{LineChannel},
    continua::AbstractVector{ContinuumChannel}, energy_grid_ev::AbstractVector{<:Real})
    line_labels = [lc.label for lc in lines]
    line_energies = [lc.energy_ev for lc in lines]
    line_rates = [_line_rate_at(lc, age_s) for lc in lines]

    continuum_labels = [cc.label for cc in continua]
    continuum_L = Matrix{Float64}(undef, length(continua), length(energy_grid_ev))
    for (i, cc) in enumerate(continua), (j, E) in enumerate(energy_grid_ev)
        continuum_L[i, j] = cc.luminosity_ev(age_s, E)
    end

    return SpectrumSnapshot(
        Float64(age_s), line_labels, line_energies, line_rates,
        collect(Float64, energy_grid_ev), continuum_labels, continuum_L,
    )
end

"""
    spectrum_timeline(ages_s, lines, continua, energy_grid_ev) -> Vector{SpectrumSnapshot}

[`composite_spectrum`](@ref) at each time in `ages_s` -- the actual
"movie" data: one frame per requested time, spanning as much of
quiescence -> TNR -> post-TNR -> freeze-out -> shock as the supplied
`lines`/`continua` cover.
"""
function spectrum_timeline(ages_s::AbstractVector{<:Real}, lines::AbstractVector{LineChannel},
    continua::AbstractVector{ContinuumChannel}, energy_grid_ev::AbstractVector{<:Real})
    return [composite_spectrum(t, lines, continua, energy_grid_ev) for t in ages_s]
end

end # module SpectralEvolution
