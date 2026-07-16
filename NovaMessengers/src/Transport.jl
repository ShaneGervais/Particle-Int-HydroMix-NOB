module Transport

using DataFrames
using ..MesaIO

export ELECTRON_REST_ENERGY_MEV, escape_probability_neutrino, klein_nishina_factor,
       compton_opacity, optical_depth_gamma, escape_probability_gamma

const ELECTRON_REST_ENERGY_MEV = 0.510998950  # MeV (CODATA electron rest mass energy)

"""
    escape_probability_neutrino(args...; kwargs...) -> 1.0

Neutrinos free-stream: their mean free path is enormously larger than the
star, so the escape probability is always 1 regardless of where they're
produced -- the [`MessengerProduction`](@ref) rate *is* the observable
signal, unmodified. Kept as an explicit function (accepting and ignoring
any arguments) so this physics decision stays visible in the code, rather
than photons and neutrinos being silently treated the same way or
neutrinos being skipped entirely.
"""
escape_probability_neutrino(args...; kwargs...) = 1.0

"""
    klein_nishina_factor(E_mev) -> Float64

Ratio of the total Klein-Nishina Compton cross section to the Thomson
cross section, at photon energy `E_mev` (MeV). Equals 1 in the low-energy
(Thomson) limit and falls off above the electron rest energy (511 keV) --
this is why different messenger lines (511 keV annihilation vs 1.275 MeV /
1.809 MeV decay lines) don't escape with the same probability from the
same depth, even though they're all "gamma rays."
"""
function klein_nishina_factor(E_mev::Real)
    E_mev > 0 || error("photon energy must be positive, got $E_mev")
    x = E_mev / ELECTRON_REST_ENERGY_MEV
    # log1p(2x), not log(1+2x): the numerator is a difference of two O(x)
    # terms that must resolve an O(x^3) result, so forming "1+2x" before
    # taking the log (rounding x into the mantissa of a number near 1)
    # loses precision catastrophically for x below ~1e-4.
    l = log1p(2x)
    term1 = (1 + x) / x^3 * (2x * (1 + x) / (1 + 2x) - l)
    term2 = l / (2x)
    term3 = (1 + 3x) / (1 + 2x)^2
    return 0.75 * (term1 + term2 - term3)
end

"""
    compton_opacity(x_h1, E_mev) -> Float64

Compton scattering opacity (cm^2/g) for a photon of energy `E_mev` in
material with hydrogen mass fraction `x_h1`. Built from the free-electron
(Thomson) opacity `0.2*(1+X)` -- the standard stellar-structure formula
for electron-scattering opacity -- times the Klein-Nishina reduction
factor at this energy. This is deliberately NOT MESA's own `log_opacity`
profile column: that's a Rosseland mean appropriate for the thermal
radiation field (mixing in bound-free and other processes), not a single
discrete MeV-scale line photon, whose dominant interaction here is Compton
scattering.
"""
function compton_opacity(x_h1::Real, E_mev::Real)
    kappa_thomson = 0.2 * (1 + x_h1)
    return kappa_thomson * klein_nishina_factor(E_mev)
end

"""
    optical_depth_gamma(profile_data, E_mev; rsun_cm) -> Vector{Float64}

Compton optical depth from each zone out to the stellar surface, for a
photon of energy `E_mev`, using the profile's own `radius` (Rsun),
`logRho`, and `x_mass_fraction_H` columns. `profile_data` must be ordered
from zone 1 at the surface inward (MESA's own convention) -- one value is
returned per zone. `rsun_cm` is required explicitly (rather than defaulted)
so callers pull it from the run's own history header (`history(run).header.rsun`)
and stay consistent with the MESA version that produced the data.
"""
function optical_depth_gamma(profile_data::AbstractDataFrame, E_mev::Real; rsun_cm::Real)
    nz = nrow(profile_data)
    r_cm = profile_data.radius .* rsun_cm
    rho = 10.0 .^ profile_data.logRho
    kappa = compton_opacity.(profile_data.x_mass_fraction_H, E_mev)

    tau = Vector{Float64}(undef, nz)
    tau[1] = 0.0
    for k in 2:nz
        dr = abs(r_cm[k - 1] - r_cm[k])
        contrib = 0.5 * (rho[k - 1] * kappa[k - 1] + rho[k] * kappa[k]) * dr
        tau[k] = tau[k - 1] + contrib
    end
    return tau
end

"""
    escape_probability_gamma(profile_data, E_mev; rsun_cm) -> Vector{Float64}
    escape_probability_gamma(run::MesaRun, profile_number, E_mev) -> Vector{Float64}

`exp(-tau_gamma)` at each zone: the first-order (pure attenuation, no
scattering redistribution) escape probability for a photon of energy
`E_mev` produced at that zone. The `MesaRun` form reads the requested
profile and pulls `rsun_cm` from the run's own history header.
"""
function escape_probability_gamma(profile_data::AbstractDataFrame, E_mev::Real; rsun_cm::Real)
    return exp.(-optical_depth_gamma(profile_data, E_mev; rsun_cm=rsun_cm))
end

function escape_probability_gamma(run::MesaRun, profile_number::Integer, E_mev::Real)
    rsun_cm = history(run).header.rsun
    profile_data = profile(run, profile_number).data
    return escape_probability_gamma(profile_data, E_mev; rsun_cm=rsun_cm)
end

end # module Transport
