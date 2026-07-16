module SignalSynthesis

using DataFrames
using ..MesaIO
using ..NuclearDecay
using ..MessengerProduction
using ..Transport

export neutrino_lightcurve, neutrino_energy_lightcurve, gamma_lightcurve,
       gamma_energy_lightcurve

const MEV_TO_ERG = 1.602176634e-6
const YEAR_S = 365.25 * 24 * 3600

"""
    neutrino_lightcurve(run::MesaRun; isotopes=keys(MESSENGER_ISOTOPES)) -> (age, rate)

Total whole-star neutrino emission rate (neutrinos/second) vs. `star_age`,
summed over `isotopes` (default: all 7 tracked messenger isotopes).
Neutrinos free-stream (`Transport.escape_probability_neutrino` is always
1), so Phase 2's [`decay_rate`](@ref) *is* the observable signal here,
unmodified by transport -- a direct, real-time probe of the burning zone,
at `history.data`'s full time cadence.
"""
function neutrino_lightcurve(run::MesaRun; isotopes=keys(MESSENGER_ISOTOPES))
    age = decay_rate(run, first(isotopes)).age
    total = zeros(Float64, length(age))
    for iso in isotopes
        total .+= decay_rate(run, iso).rate
    end
    return (age=age, rate=total)
end

"""
    neutrino_energy_lightcurve(run::MesaRun; isotopes=keys(MESSENGER_ISOTOPES)) -> (age, rate)

Total neutrino energy-loss rate (erg/second) vs. `star_age`, summed over
`isotopes`. Same free-streaming logic as [`neutrino_lightcurve`](@ref).
"""
function neutrino_energy_lightcurve(run::MesaRun; isotopes=keys(MESSENGER_ISOTOPES))
    age = decay_rate(run, first(isotopes)).age
    total = zeros(Float64, length(age))
    for iso in isotopes
        total .+= neutrino_energy_loss_rate(run, iso).rate .* MEV_TO_ERG
    end
    return (age=age, rate=total)
end

"""
    gamma_lightcurve(run::MesaRun, line_energy_mev::Real; isotopes=keys(MESSENGER_ISOTOPES)) -> (age, rate)

Observable annihilation-line photon escape rate (photons/second) at
`line_energy_mev`, vs. `star_age`: Phase 2's per-zone photon production
([`zone_annihilation_photon_rate`](@ref), summed over `isotopes`) weighted
zone-by-zone by Phase 3's [`escape_probability_gamma`](@ref), then summed
over zones, at each saved profile snapshot. This is the plot that directly
answers "does the signal preserve or reshape the underlying nuclear
information": compare its shape/timing against [`neutrino_lightcurve`](@ref)
from the same burst.

Coarser in time than `neutrino_lightcurve` (limited to profile save
points, not every history row), because it needs the zone-resolved
density/composition structure Phase 3's escape probability depends on.
Currently only 511 keV is physically populated by the tracked isotope
set (all 7 are beta+ emitters); other line energies are supported for
when the NuPPN trajectory extension adds isotopes with their own decay
lines (22Na at 1.275 MeV, 26Al at 1.809 MeV).
"""
function gamma_lightcurve(run::MesaRun, line_energy_mev::Real; isotopes=keys(MESSENGER_ISOTOPES))
    pidx = profiles_index(run)
    hist = history(run).data
    rsun_cm = history(run).header.rsun

    age = Vector{Float64}(undef, nrow(pidx))
    rate = Vector{Float64}(undef, nrow(pidx))
    for (i, row) in enumerate(eachrow(pidx))
        prof = profile(run, row.profile_number).data
        esc = escape_probability_gamma(prof, line_energy_mev; rsun_cm=rsun_cm)

        zone_total = zeros(Float64, nrow(prof))
        for iso in isotopes
            zone_total .+= zone_annihilation_photon_rate(prof, iso)
        end
        rate[i] = sum(zone_total .* esc)

        j = findfirst(==(row.model_number), hist.model_number)
        age[i] = j === nothing ? NaN : hist.star_age[j] * YEAR_S
    end
    order = sortperm(age)
    return (age=age[order], rate=rate[order])
end

"""
    gamma_energy_lightcurve(run::MesaRun, line_energy_mev::Real; isotopes=keys(MESSENGER_ISOTOPES)) -> (age, rate)

Observable line luminosity (erg/second) at `line_energy_mev`: `gamma_lightcurve`
converted from photons/second to erg/second.
"""
function gamma_energy_lightcurve(run::MesaRun, line_energy_mev::Real; isotopes=keys(MESSENGER_ISOTOPES))
    g = gamma_lightcurve(run, line_energy_mev; isotopes=isotopes)
    return (age=g.age, rate=g.rate .* line_energy_mev .* MEV_TO_ERG)
end

end # module SignalSynthesis
