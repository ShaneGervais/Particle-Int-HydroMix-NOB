module ShockAcceleration

using ..Transport: ELECTRON_REST_ENERGY_MEV
using ..MesaIO

export ShockModelParams, calibrate_shock_params,
       wind_mass_loss_rate, wind_velocity, shell_mass, shell_velocity, shock_velocity,
       shock_radius, shock_power, shock_temperature, postshock_density, shell_temperature,
       column_density, shell_density,
       pion_creation_timescale, cosmic_ray_luminosity, cosmic_ray_energy,
       gamma_ray_luminosity, gamma_ray_luminosity_calorimetric,
       radiative_cooling_length_ratio, downstream_thickness, postshock_field,
       max_proton_energy_gev, max_proton_energy_gev_closedform,
       proton_spectrum, gamma_ray_spectrum,
       bethe_heitler_cross_section, bethe_heitler_optical_depth,
       optical_radiation_energy_density, gamma_gamma_optical_depth,
       bremsstrahlung_emissivity_nu, shock_bremsstrahlung_luminosity_nu,
       shock_bremsstrahlung_luminosity_ev

# Physical constants (cgs). MSUN_G matches this project's own MESA run's
# header `msun` value (frozen-from-MESA principle, same as elsewhere in
# this package) rather than a generic external solar mass.
const MSUN_G = 1.9884098706980504e33
const PROTON_MASS_G = 1.67262192369e-24        # CODATA
const BOLTZMANN_ERG_K = 1.380649e-16           # exact (2019 SI)
const C_CM_S = 2.99792458e10                   # exact (SI definition)
const SIGMA_SB = 5.670374419e-5                # erg cm^-2 s^-1 K^-4, exact-derived
const ELEMENTARY_CHARGE_ESU = 4.80320425e-10   # esu (Gaussian-cgs)
const GEV_TO_ERG = 1.602176634e-3              # exact (2019 SI elementary charge)
const DAY_S = 86400.0
const KM_S_TO_CM_S = 1.0e5
const G_CGS = 6.67430e-8                       # cm^3 g^-1 s^-2 (CODATA)

# Model-specific constants, given explicitly in Diesing & Metzger (2026),
# "A Unified Model for Shock Interaction and gamma-Ray Emission in
# Classical Novae" (references/model_shock_classical_novae.pdf):
const SIGMA_PP_CM2 = 5.0e-26   # proton-proton cross section, text near Eq (19)
const LAMBDA0 = 2.0e-27        # free-free cooling coefficient, Lambda(T) = LAMBDA0*sqrt(T), Eq (26)
const SIGMA_T_CM2 = 6.65e-25   # Thomson cross section, text near Eq (37)
const ALPHA_FS = 1.0 / 137.0   # fine structure constant, "alpha ~= 1/137" per text near Eq (37)

# Proton rest energy (GeV), derived from this module's own PROTON_MASS_G --
# self-consistent with the proton mass used everywhere else in this module,
# rather than a separately-copied literature value (~0.938 GeV).
const PROTON_REST_ENERGY_GEV = PROTON_MASS_G * C_CM_S^2 / GEV_TO_ERG

# rL[cm] = LARMOR_COEFF_CM * E[GeV] / B[G] for an ultra-relativistic
# proton (rL = pc/(eB) ~ E/(eB)), derived from GEV_TO_ERG/ELEMENTARY_CHARGE_ESU
# rather than copied from the paper's own rounded "3.3e6" (Eq (30)) --
# matches to the stated precision (3.3356e6 vs their 3.3e6), independently
# re-derived here for one more significant figure.
const LARMOR_COEFF_CM = GEV_TO_ERG / ELEMENTARY_CHARGE_ESU

"""
    ShockModelParams(; Menv_Msun=1e-4, tau_days=20.0, vf_km_s=6000.0,
                        fX=5e-5, delta_s_over_Rs=1e-2, fOmega=0.3,
                        xiCR=0.03, xiB=0.01, kappa=0.1, Mwd_Msun=1.0)

Parameters for Diesing & Metzger's parameterized ("toy") model of shock
interaction and gamma-ray emission in classical novae. Fiducial defaults
and ranges are their Table 1. Stored internally in cgs; construct with
the natural-unit keywords above.

This model describes the *ejection* phase of a nova outflow (a fast wind
colliding with a slower shell released earlier) -- it is deliberately
NOT derived from `mesa_work/wd_nova_burst_co`'s own `LOGS/` output,
because that run's ~480-step window captures the thermonuclear runaway
and initial decline but terminates (`extras_check_model`) before the
wind-driven mass-ejection phase these parameters (`Menv`, `tau`, `vf`)
actually describe -- `star_mass` in that run's history data changes by
only ~1e-5 Msun over the whole run. Use this module standalone with
literature-typical or hand-chosen parameters until/unless a MESA run
extending into the ejection phase is available to calibrate it from.

Fields (all cgs): `Menv` (g), `tau` (s), `vf` (cm/s), `fX`
(X-ray efficiency L_X/L_sh, dimensionless), `delta_s_over_Rs`
(shell thickness / shock radius, dimensionless), `fOmega` (shell
covering fraction, dimensionless), `xiCR` (shock -> cosmic ray
acceleration efficiency, dimensionless), `xiB` (shock ram pressure ->
B-field energy fraction, dimensionless), `kappa` (fraction of cosmic ray
energy radiated as gamma-rays per p-p interaction, dimensionless), `Mwd` (g).
"""
struct ShockModelParams
    Menv::Float64
    tau::Float64
    vf::Float64
    fX::Float64
    delta_s_over_Rs::Float64
    fOmega::Float64
    xiCR::Float64
    xiB::Float64
    kappa::Float64
    Mwd::Float64
end

function ShockModelParams(;
    Menv_Msun::Real=1e-4, tau_days::Real=20.0, vf_km_s::Real=6000.0,
    fX::Real=5e-5, delta_s_over_Rs::Real=1e-2, fOmega::Real=0.3,
    xiCR::Real=0.03, xiB::Real=0.01, kappa::Real=0.1, Mwd_Msun::Real=1.0,
)
    return ShockModelParams(
        Menv_Msun * MSUN_G, tau_days * DAY_S, vf_km_s * KM_S_TO_CM_S,
        fX, delta_s_over_Rs, fOmega, xiCR, xiB, kappa, Mwd_Msun * MSUN_G,
    )
end

"""
    calibrate_shock_params(run::MesaRun; ignition_mass_coordinate, tau_days=20.0,
                            fX=5e-5, delta_s_over_Rs=1e-2, fOmega=0.3, xiCR=0.03,
                            xiB=0.01, kappa=0.1, history_row=1) -> ShockModelParams

Builds [`ShockModelParams`](@ref) with `Mwd_Msun`, `Menv_Msun`, and
`vf_km_s` set from this run's own MESA data instead of Diesing &
Metzger's fiducial defaults:

    Mwd_Msun  = star_mass(history_row)
    Menv_Msun = max(Mwd_Msun - ignition_mass_coordinate, 0.0)
    vf_km_s   = sqrt(2 G Mwd_g / R_cm) / KM_S_TO_CM_S      (surface escape velocity)

- `Mwd_Msun`: this run's own WD mass, at `history_row` (default 1,
  pre-TNR).
- `Menv_Msun`: the mass of material *above* wherever the TNR actually
  ignites is the natural proxy for "how much envelope this event could
  plausibly eject" -- everything below that depth is unburned WD
  substrate, not part of the flash. Pass
  `TrajectoryPostProcessing.peak_temperature_zone(run).mass_coordinate`
  for `ignition_mass_coordinate`.
- `vf_km_s`: real nova ejecta velocities are observed to scale with the
  underlying WD's escape velocity (more compact/massive WDs -> faster
  ejecta), so the WD's own surface escape velocity is a physically
  motivated anchor, not an arbitrary substitute for the paper's
  literature-typical 6000 km/s.

`tau_days`/`fX`/`delta_s_over_Rs`/`fOmega`/`xiCR`/`xiB`/`kappa` remain
fiducial (paper defaults, overridable via keyword) -- nothing in a
single-zone TNR run informs the ejection *timescale* or the shock
microphysics efficiencies, so this is a partial calibration, not a full
replacement for an actual MESA-resolved ejection run (see this module's
own struct docstring for why `wd_nova_burst_co` can't supply that).
"""
function calibrate_shock_params(run::MesaRun; ignition_mass_coordinate::Real,
    tau_days::Real=20.0, fX::Real=5e-5, delta_s_over_Rs::Real=1e-2, fOmega::Real=0.3,
    xiCR::Real=0.03, xiB::Real=0.01, kappa::Real=0.1, history_row::Integer=1)
    h = history(run).data
    # Use the run's own maximum recorded mass, not history_row's, for the
    # Menv subtraction: star_mass grows slightly over the run (ongoing
    # accretion), so a mass_coordinate measured from a *later* profile
    # (as peak_temperature_zone's typically is) can exceed an earlier
    # history_row's star_mass, making Menv spuriously clamp to zero --
    # confirmed: this produced Menv=0 -> wind_mass_loss_rate=0 ->
    # radiative_cooling_length_ratio=Inf (division by zero) ->
    # downstream_thickness=Inf -> a 0*Inf = NaN bremsstrahlung volume.
    Mwd_Msun = maximum(h.star_mass)
    Menv_Msun = max(Mwd_Msun - ignition_mass_coordinate, 0.0)
    R_cm = h.radius_cm[history_row]
    Mwd_g = Mwd_Msun * MSUN_G
    v_esc_cm_s = sqrt(2 * G_CGS * Mwd_g / R_cm)
    vf_km_s = v_esc_cm_s / KM_S_TO_CM_S
    return ShockModelParams(;
        Menv_Msun=Menv_Msun, tau_days=tau_days, vf_km_s=vf_km_s,
        fX=fX, delta_s_over_Rs=delta_s_over_Rs, fOmega=fOmega,
        xiCR=xiCR, xiB=xiB, kappa=kappa, Mwd_Msun=Mwd_Msun,
    )
end

# --- Hydrodynamics (Eq 1-14) ------------------------------------------------

"""wind_mass_loss_rate(p, t) -> g/s. Eq (1): Mdot_w(t) = (Menv/tau) exp(-t/tau)."""
wind_mass_loss_rate(p::ShockModelParams, t::Real) = (p.Menv / p.tau) * exp(-t / p.tau)

"""wind_velocity(p, t) -> cm/s. Eq (2): v_w(t) = [1 - exp(-t/tau)] vf."""
wind_velocity(p::ShockModelParams, t::Real) = (1 - exp(-t / p.tau)) * p.vf

"""shell_mass(p, t) -> g. Eq (3): Ms(t) = [1 - exp(-t/tau)] Menv."""
shell_mass(p::ShockModelParams, t::Real) = (1 - exp(-t / p.tau)) * p.Menv

"""shell_velocity(p, t) -> cm/s. Eq (5): vs(t) = v_w(t)/2 (momentum-conserving reverse shock)."""
shell_velocity(p::ShockModelParams, t::Real) = wind_velocity(p, t) / 2

"""shock_velocity(p, t) -> cm/s. Eq (7): v_sh = v_w - v_s = v_w/2 = v_s."""
shock_velocity(p::ShockModelParams, t::Real) = shell_velocity(p, t)

"""
shock_radius(p, t) -> cm. Eq (6): Rs(t) = (vf/2)[t - tau + tau*exp(-t/tau)].
Domain: t > 0 (Rs -> 0 as t -> 0, as t^2; several downstream quantities
divide by Rs^2 and diverge in that limit -- physically sensible, since an
infinitesimally young shell has no well-defined density, but means this
model is only meaningful evaluated at t appreciably greater than 0, not
at the exact origin).
"""
function shock_radius(p::ShockModelParams, t::Real)
    x = t / p.tau
    return (p.vf / 2) * p.tau * (x - 1 + exp(-x))
end

"""shock_power(p, t) -> erg/s. Eq (8): Lsh = (1/2) Mdot_w (v_w^2 - v_s^2)."""
function shock_power(p::ShockModelParams, t::Real)
    vw = wind_velocity(p, t)
    vs = shell_velocity(p, t)
    return 0.5 * wind_mass_loss_rate(p, t) * (vw^2 - vs^2)
end

"""shock_temperature(p, t) -> K. Eq (10): Tsh = (3 mp / 16 k) v_sh^2."""
function shock_temperature(p::ShockModelParams, t::Real)
    vsh = shock_velocity(p, t)
    return (3 * PROTON_MASS_G / (16 * BOLTZMANN_ERG_K)) * vsh^2
end

"""postshock_density(p, t) -> g/cm^3. Eq (11): rho_sh = 4*rho_w(Rs) = Mdot_w/(pi Rs^2 v_w)."""
function postshock_density(p::ShockModelParams, t::Real)
    Rs = shock_radius(p, t)
    return wind_mass_loss_rate(p, t) / (pi * Rs^2 * wind_velocity(p, t))
end

"""
shell_temperature(p, t) -> K. Eq (12): Ts = [(Lsh + Lwd)/(4 pi sigma_SB Rs^2)]^(1/4),
with Lwd ~ L_Edd ~ 1.4e38 (Mwd/Msun) erg/s (white dwarf's own persistent luminosity).
"""
function shell_temperature(p::ShockModelParams, t::Real)
    Lwd = 1.4e38 * (p.Mwd / MSUN_G)
    Rs = shock_radius(p, t)
    return ((shock_power(p, t) + Lwd) / (4 * pi * SIGMA_SB * Rs^2))^0.25
end

"""column_density(p, t) -> cm^-2. Eq (13): N_H = Ms/(4 pi mp Rs^2)."""
function column_density(p::ShockModelParams, t::Real)
    Rs = shock_radius(p, t)
    return shell_mass(p, t) / (4 * pi * PROTON_MASS_G * Rs^2)
end

"""shell_density(p, t) -> g/cm^3. Eq (14): rho_s = Ms/(4 pi Delta_s Rs^2), Delta_s = delta_s_over_Rs * Rs."""
function shell_density(p::ShockModelParams, t::Real)
    Rs = shock_radius(p, t)
    delta_s = p.delta_s_over_Rs * Rs
    return shell_mass(p, t) / (4 * pi * delta_s * Rs^2)
end

# --- Cosmic rays / gamma-ray luminosity (Eq 19-25) --------------------------

"""
pion_creation_timescale(p, t) -> s. t_pi = (n_s sigma_pp c)^-1, with
n_s = rho_s/mp -- derived from Eq (14) and the definition below Eq (19);
matches the closed form in Eq (21)'s numerator/denominator.
"""
function pion_creation_timescale(p::ShockModelParams, t::Real)
    Rs = shock_radius(p, t)
    delta_s = p.delta_s_over_Rs * Rs
    return 4 * pi * PROTON_MASS_G * delta_s * Rs^2 / (shell_mass(p, t) * SIGMA_PP_CM2 * C_CM_S)
end

"""cosmic_ray_luminosity(p, t) -> erg/s. L_CR = xiCR * Lsh (text preceding Eq 19)."""
cosmic_ray_luminosity(p::ShockModelParams, t::Real) = p.xiCR * shock_power(p, t)

"""
    cosmic_ray_energy(p::ShockModelParams, t_grid::AbstractVector) -> Vector{Float64}

Solve Eq (19), dECR/dt = LCR - ECR/t_pi - ECR/t, for the total cosmic-ray
energy confined in the shell (erg), at each time in `t_grid` (s, strictly
increasing). `t_grid[1]` should be small compared to `p.tau` (e.g.
1e-3*tau or smaller) -- the 1/t adiabatic-loss term is singular at t=0,
so integration effectively starts from ECR=0 at `t_grid[1]` rather than
from t=0 itself, which is an excellent approximation when t_grid[1] << tau
since LCR ~ t^2 near t=0 (negligible CR energy has accumulated yet).

This ODE is a *fast-relaxation* problem, not a smooth/non-stiff one: the
whole point of the "calorimetric limit" (Eq 21-23) is t_pi << t through
most of the physically relevant range (the paper finds t_pi/t ~ 2e-4 at
t_pk for fiducial parameters), so the loss rate 1/t_pi is enormous
compared to 1/t and to any reasonable step size. A fixed-step explicit
method (e.g. RK4) is numerically unstable against a loss term that
fast -- it does not converge by shrinking the step, it diverges. Instead
this uses an exponential (integrating-factor) step: treating the loss
rate gamma(t) = 1/t_pi(t) + 1/t and source L(t) = LCR(t) as constant at
their midpoint value over each small interval, the *exact* local solution
is E(t1) = E(t0)*exp(-gamma*h) + (L/gamma)*(1-exp(-gamma*h)) -- stable for
any gamma*h, and it correctly reduces to the calorimetric quasi-equilibrium
ECR -> L/gamma = LCR*t_pi (Eq 22's regime) when gamma*h is large, rather
than blowing up.
"""
function cosmic_ray_energy(p::ShockModelParams, t_grid::AbstractVector{<:Real})
    n = length(t_grid)
    n >= 2 || error("t_grid must have at least 2 points")
    ECR = zeros(Float64, n)
    for i in 2:n
        t0, t1 = t_grid[i - 1], t_grid[i]
        tm = 0.5 * (t0 + t1)
        h = t1 - t0
        E0 = ECR[i - 1]
        L = cosmic_ray_luminosity(p, tm)
        gamma = 1 / pion_creation_timescale(p, tm) + 1 / tm
        gh = gamma * h
        if gh < 1e-8
            # avoid catastrophic cancellation in (1-exp(-x))/x for tiny x
            ECR[i] = E0 + h * (L - gamma * E0)
        else
            decay = exp(-gh)
            ECR[i] = E0 * decay + (L / gamma) * (1 - decay)
        end
    end
    return ECR
end

"""
gamma_ray_luminosity(p, t, ECR) -> erg/s. Eq (20): Lgamma = fOmega*kappa*ECR/t_pi,
given an already-solved ECR(t) (see [`cosmic_ray_energy`](@ref)).
"""
gamma_ray_luminosity(p::ShockModelParams, t::Real, ECR::Real) =
    p.fOmega * p.kappa * ECR / pion_creation_timescale(p, t)

"""
gamma_ray_luminosity_calorimetric(p, t) -> erg/s. Eq (22): the calorimetric-limit
(t_pi << t) steady-state approximation, Lgamma ~= fOmega*xiCR*kappa*Lsh --
faithfully tracks the shock power directly, no ODE needed. Valid near peak
shock power for typical nova parameters (Eq 23); breaks down once cosmic
rays enter the adiabatic limit at late times (Eq 25), where the full
[`cosmic_ray_energy`](@ref) + [`gamma_ray_luminosity`](@ref) ODE solution
should be used instead.
"""
gamma_ray_luminosity_calorimetric(p::ShockModelParams, t::Real) =
    p.fOmega * p.xiCR * p.kappa * shock_power(p, t)

# --- Maximum particle energy (Eq 26-32) -------------------------------------

"""
radiative_cooling_length_ratio(p, t) -> dimensionless. Delta_rad/Rs, Eq (26):
the thickness (relative to Rs) that the immediate post-shock region would
have if it cooled radiatively in a laminar flow -- an upper bound on the
true (turbulence-suppressed) hot-layer thickness; see [`downstream_thickness`](@ref).

This implements the *middle* form of Eq (26), (3*pi/16)*k*Tsh*mp*vw^2*Rs/(Mdot_w*Lambda(Tsh))
-- verified symbol-for-symbol against both the published PDF and the arXiv
HTML/MathJax rendering (arXiv:2604.06310v1), including v_down = v_sh/4 = v_w/8
and Lambda_0 = 2e-27. Deliberately NOT the paper's own further-reduced
closed form (their quoted "5.29e-2 * Menv,-4^-1 * tau20^2 * vsh,8^4 * (t/tau) * exp(t/tau)"):
substituting this module's exact hydrodynamics into that middle form gives
a coefficient of ~1.95e-2 at v_sh fully saturated (large t), not 5.29e-2 --
confirmed independently three ways (direct symbolic reduction, and
back-solving the ratio this function's own output settles to at large
t/tau in examples/shock_model_validation.jl's Check 6), with no erratum
found for the paper (single-version, Apr 2026). The ~2.72 (~e) ratio
between the two candidate coefficients is unexplained; this function
follows the more primary, independently-reverifiable middle form rather
than the compressed closed form. If this ever gets resolved (e.g. the
paper issues a correction, or the missing factor is found), the fix
belongs here.
"""
function radiative_cooling_length_ratio(p::ShockModelParams, t::Real)
    Tsh = shock_temperature(p, t)
    vw = wind_velocity(p, t)
    Rs = shock_radius(p, t)
    Lambda_Tsh = LAMBDA0 * sqrt(Tsh)
    return (3 * pi / 16) * (BOLTZMANN_ERG_K * Tsh * PROTON_MASS_G * vw^2 * Rs) /
           (wind_mass_loss_rate(p, t) * Lambda_Tsh)
end

"""
downstream_thickness(p, t) -> cm. Eq (27): Delta_down = (32/9) fX fOmega^-1 Delta_rad
-- the true (empirically X-ray-efficiency-calibrated) thickness of the hot
post-shock acceleration region, much smaller than the naive radiative
cooling length ([`radiative_cooling_length_ratio`](@ref)) once turbulent
mixing with the cool shell is accounted for (Metzger et al. 2025).
"""
function downstream_thickness(p::ShockModelParams, t::Real)
    Rs = shock_radius(p, t)
    return (32 / 9) * p.fX / p.fOmega * radiative_cooling_length_ratio(p, t) * Rs
end

"""
postshock_field(p, t) -> Gauss. Eq (31): B_down = sqrt(xiB * Mdot_w * v_sh / Rs^2),
from assuming a fraction xiB of the shock's ram pressure goes into the
post-shock magnetic field's energy density.
"""
function postshock_field(p::ShockModelParams, t::Real)
    Rs = shock_radius(p, t)
    return sqrt(p.xiB * wind_mass_loss_rate(p, t) * shock_velocity(p, t) / Rs^2)
end

"""
    max_proton_energy_gev(p, t) -> GeV

Maximum proton energy from the Hillas-like confinement argument, Eq (29)+(30):
assuming Bohm diffusion, the downstream diffusion length
`4*c*rL(Emax)/(3*v_sh)` is set equal to the acceleration region's thickness
[`downstream_thickness`](@ref) and solved for `Emax`, using the Larmor
radius `rL[cm] = LARMOR_COEFF_CM * E[GeV] / B[G]` for an ultra-relativistic
proton in the field from [`postshock_field`](@ref).

DOMAIN OF VALIDITY WARNING: this grows unphysically large (far past the
Planck energy) if evaluated at `t` more than a handful of `tau` past
`t_pk` -- [`downstream_thickness`](@ref) inherits an `(t/tau)*exp(t/tau)`
growth factor (Eq 26) that the paper itself never plots or claims
validity for beyond a "couple months" (~3-5 tau; their own Figure 2 stops
there). Cross-checked against the paper's own worked examples: matches
its `t~t_pk` "~12 GeV" estimate only qualitatively (this function uses
the *exact* time-dependent v_sh(t), not the paper's `t>>tau` asymptotic
v_sh->vf/2 substitution, so the two aren't expected to agree precisely
at `t_pk` itself -- see `examples/shock_model_validation.jl` Check 6 for
the full diagnostic), and its `t~2.5*tau` "~10 TeV" checkpoint to within
a factor of ~2 (Check 7) -- but do not trust this function extrapolated
to `t >> a few tau` without further validation; the source of the
continued raw-vs-asymptotic disagreement even after v_sh fully saturates
(Check 6 settles near a ratio of ~0.36, not 1, for `t >~ 10*tau`) has not
been fully isolated.
"""
function max_proton_energy_gev(p::ShockModelParams, t::Real)
    vsh = shock_velocity(p, t)
    Ddown = downstream_thickness(p, t)
    B = postshock_field(p, t)
    rL = 3 * vsh * Ddown / (4 * C_CM_S)
    return rL * B / LARMOR_COEFF_CM
end

"""
    max_proton_energy_gev_closedform(p, t) -> GeV

Same physical quantity as [`max_proton_energy_gev`](@ref), computed
directly from Eq (32)'s general form, Emax = (3/4)(v_sh/c) e B_down Delta_down,
rather than via the rL/B rearrangement. Algebraically identical to
`max_proton_energy_gev` -- agreement between the two is an internal
arithmetic self-consistency check of this module, not an independent
validation against the paper (for that, see the worked fiducial-parameter
example in `examples/shock_model_validation.jl`).
"""
function max_proton_energy_gev_closedform(p::ShockModelParams, t::Real)
    vsh = shock_velocity(p, t)
    Ddown = downstream_thickness(p, t)
    B = postshock_field(p, t)
    Emax_erg = 0.75 * (vsh / C_CM_S) * ELEMENTARY_CHARGE_ESU * B * Ddown
    return Emax_erg / GEV_TO_ERG
end

# --- Gamma-ray spectra (Eq 34-36) -------------------------------------------

"""
Integrate `f` from `a` to `b` (both > 0) over `n` log-spaced points via the
trapezoidal rule. The injection-spectrum normalization integral below is a
smooth power-law-times-exponential-cutoff integrand spanning many decades
in energy -- log-spaced quadrature handles that well without needing an
incomplete-gamma/exponential-integral special function (and hence without
adding SpecialFunctions.jl as a dependency).
"""
function _log_trapz(f, a::Real, b::Real, n::Integer=400)
    logs = range(log(a), log(b); length=n)
    xs = exp.(logs)
    total = 0.0
    fprev = f(xs[1])
    for i in 2:n
        fcur = f(xs[i])
        total += 0.5 * (fcur + fprev) * (xs[i] - xs[i - 1])
        fprev = fcur
    end
    return total
end

"""
Normalization `Ai` (protons / GeV) of an injected spectrum
`Ai*(E/E0)^-2*exp(-E/Emax)` (E0 = 1 GeV) such that its energy integral from
the proton rest energy to `50*Emax` (the exponential cutoff makes anything
beyond this negligible) equals `target_energy_gev` -- Eq (34)'s normalization
constraint, `integral(E*phi_i dE, mp*c^2, Inf) = LCR(ti)*dti`.
"""
function _injection_normalization(Emax_gev::Real, target_energy_gev::Real)
    integrand(E) = E * E^-2 * exp(-E / Emax_gev)  # E0 = 1 GeV folded in
    denom = _log_trapz(integrand, PROTON_REST_ENERGY_GEV, 50 * Emax_gev)
    return target_energy_gev / denom
end

"""
    proton_spectrum(p, E_grid_gev, t_now, injection_grid) -> Vector{Float64}

phi(E, t_now) [protons/GeV], Eq (34)-(35): the cumulative proton spectrum at
`t_now` (s) from superposing power-law-with-cutoff spectra injected during
each epoch of `injection_grid` (s, strictly increasing, all entries
<= t_now -- each consecutive pair `injection_grid[i], injection_grid[i+1]`
is one epoch, injected at its midpoint), each evolved from its injection
time to `t_now` via adiabatic dilation (`Lad = Rs(t_now)/Rs(t_i) >= 1`,
Eq 35) and proton-proton losses (`Lpi = exp(-(t_now-t_i)/t_pi_avg)`, with
`t_pi_avg` approximated by `t_pi` at the epoch-to-now midpoint time).
"""
function proton_spectrum(p::ShockModelParams, E_grid::AbstractVector{<:Real}, t_now::Real,
    injection_grid::AbstractVector{<:Real})
    phi = zeros(Float64, length(E_grid))
    Rs_now = shock_radius(p, t_now)
    n_epochs = length(injection_grid) - 1
    n_epochs >= 1 || error("injection_grid must have at least 2 points")

    for idx in 1:n_epochs
        ti, ti_next = injection_grid[idx], injection_grid[idx + 1]
        ti_next > t_now && break
        dt = ti_next - ti
        tmid = 0.5 * (ti + ti_next)

        Emax_i = max_proton_energy_gev(p, tmid)
        LCR_i_dt_gev = cosmic_ray_luminosity(p, tmid) * dt / GEV_TO_ERG
        Ai = _injection_normalization(Emax_i, LCR_i_dt_gev)

        Rs_i = shock_radius(p, tmid)
        Lad = Rs_now / Rs_i
        t_avg = 0.5 * (tmid + t_now)
        Lpi = exp(-(t_now - tmid) / pion_creation_timescale(p, t_avg))

        for (k, E) in enumerate(E_grid)
            Eprime = Lad * E
            phi[k] += Lad * Lpi * Ai * Eprime^-2 * exp(-Eprime / Emax_i)
        end
    end
    return phi
end

"""
    gamma_ray_spectrum(p, E_grid_gev, t_now, injection_grid) -> Vector{Float64}

dNgamma/dE [photons/GeV] at `t_now`, Eq (36): `fOmega * phi(E/kappa, t_now) /
(kappa * t_pi(t_now))`, built from [`proton_spectrum`](@ref).

CONSISTENCY CAVEAT: the paper states integrating `E*dNgamma/dE` should
recover [`gamma_ray_luminosity`](@ref)'s ODE-based bolometric `Lgamma(t_now)`
to "within 1%". This implementation lands within a factor of ~1.2-1.7x
instead (see `examples/shock_spectrum_validation.jl`), for two identified,
non-overlapping reasons rather than one bug: (1) Eq (19)'s ODE adiabatic
loss term `ECR/t` implicitly assumes `Rs ~ t` (`d ln Rs/d ln t -> 1`),
which only holds once `t >> tau`, whereas [`proton_spectrum`](@ref)'s
`Lad = Rs(t_now)/Rs(t_i)` uses the *exact* `Rs(t)` at any time -- this
alone explains most of the gap near `t_pk` and shrinks as `t_now` moves
deeper past `tau` (confirmed: the mismatch improves from 1.49x at
`t=1.1*tau` to 1.24x at `t=3*tau`, tracking `t*d(ln Rs)/dt` converging to
1). (2) it then gets *worse* again at `t_now >~ 5*tau` (1.29x at `5*tau`,
1.72x at `8*tau`) because each summed epoch's [`max_proton_energy_gev`](@ref)
inherits the already-documented `Delta_rad/Rs` coefficient discrepancy
(see that function's docstring) -- epochs at large `t/tau` feed an
increasingly unreliable cutoff into the spectrum. Treat this function's
bolometric normalization as good to an order of magnitude / factor of
~2, not the paper's claimed 1%, until the upstream Eq (26) question
resolves; the spectral *shape* (power law + cutoff location) is on firmer
footing than the absolute normalization.
"""
function gamma_ray_spectrum(p::ShockModelParams, E_grid::AbstractVector{<:Real}, t_now::Real,
    injection_grid::AbstractVector{<:Real})
    phi_vals = proton_spectrum(p, E_grid ./ p.kappa, t_now, injection_grid)
    tpi_now = pion_creation_timescale(p, t_now)
    return (p.fOmega / (p.kappa * tpi_now)) .* phi_vals
end

# --- Gamma-ray absorption (Eq 37-40) ----------------------------------------

"""
    bethe_heitler_cross_section(Eg_gev; Z=7) -> cm^2

Eq (37): Bethe-Heitler pair-production cross section for a gamma-ray of
energy `Eg_gev` interacting with a nucleus of charge `Z` (default 7,
nitrogen -- representative of CNO-dominated nova ejecta per the paper's
own worked example).
"""
function bethe_heitler_cross_section(Eg_gev::Real; Z::Real=7)
    me_c2_gev = ELECTRON_REST_ENERGY_MEV * 1e-3
    return (3 / (8 * pi)) * ALPHA_FS * SIGMA_T_CM2 * Z^2 *
           ((28 / 9) * log(2 * Eg_gev / me_c2_gev) - 218 / 27)
end

"""
    bethe_heitler_optical_depth(p, t, Eg_gev; A=14, Z=7) -> dimensionless

Eq (38): tau_BH = sigma_BH * (N_H/A), the pair-production optical depth
through the cool shell for a gamma-ray of energy `Eg_gev`, assuming ejecta
nuclei of mass number `A` and charge `Z` (defaults A=14, Z=7: nitrogen,
representative of CNO-dominated ejecta, matching the paper's own worked
example -- `N_H` from [`column_density`](@ref) is a hydrogen-mass-equivalent
column, so `N_H/A` converts it to an actual nuclei column).
"""
function bethe_heitler_optical_depth(p::ShockModelParams, t::Real, Eg_gev::Real; A::Real=14, Z::Real=7)
    return bethe_heitler_cross_section(Eg_gev; Z=Z) * column_density(p, t) / A
end

"""
optical_radiation_energy_density(p, t) -> erg/cm^3. Eq (39): u_gamma =
Lsh/(4*pi*c*Rs^2), the seed-photon energy density from shock power
reprocessed into optical/UV emission, used by [`gamma_gamma_optical_depth`](@ref).
"""
optical_radiation_energy_density(p::ShockModelParams, t::Real) =
    shock_power(p, t) / (4 * pi * C_CM_S * shock_radius(p, t)^2)

"""
    gamma_gamma_optical_depth(p, t, Eg_gev) -> dimensionless

Eq (40): pair-production optical depth for a gamma-ray of energy `Eg_gev`
against the nova's own reprocessed optical/UV light (seed photon energy
`Eopt = (me*c^2)^2/Eg`, cross section `sigma_gg ~= sigma_T/4`). Only
significant for `Eg` at the TeV scale interacting with ~0.25 eV optical
seed photons -- negligible for `Eg <~ 100 GeV` (see paper's discussion
following Eq 40).
"""
function gamma_gamma_optical_depth(p::ShockModelParams, t::Real, Eg_gev::Real)
    me_c2_gev = ELECTRON_REST_ENERGY_MEV * 1e-3
    Eopt_gev = me_c2_gev^2 / Eg_gev
    Eopt_erg = Eopt_gev * GEV_TO_ERG
    n_seed = optical_radiation_energy_density(p, t) / Eopt_erg
    sigma_gg = SIGMA_T_CM2 / 4
    return n_seed * shock_radius(p, t) * sigma_gg
end

# --- Thermal bremsstrahlung (shock-heated plasma continuum) -----------------

const H_ERG_S = 6.62607015e-27     # erg s (exact, 2019 SI) -- Planck constant
const EV_TO_ERG = 1.602176634e-12  # erg/eV (exact)
const GAUNT_FACTOR = 1.2           # order-unity thermal-averaged Gaunt factor, standard approximation

"""
    bremsstrahlung_emissivity_nu(n_e, n_i, T, nu; Z=1) -> erg/s/cm^3/Hz

Thermal (free-free) volume emissivity at frequency `nu` (Hz) from a fully
ionized plasma of electron density `n_e`, ion density `n_i` (cm^-3), and
temperature `T` (K):

    epsilon_ff(nu) = 6.8e-38 * Z^2 * n_e * n_i * T^(-1/2) * exp(-h*nu/(k T)) * g_ff

(Rybicki & Lightman's standard cgs thermal bremsstrahlung emissivity; the
`6.8e-38` coefficient already has the frequency-independent physical
constants folded in). `g_ff` is the frequency/temperature-averaged Gaunt
factor, taken as a fixed `GAUNT_FACTOR=1.2` (order-unity correction,
standard simplification -- not resolved as a function of `nu`/`T` here).
"""
function bremsstrahlung_emissivity_nu(n_e::Real, n_i::Real, T::Real, nu::Real; Z::Real=1)
    return 6.8e-38 * Z^2 * n_e * n_i * T^-0.5 * exp(-H_ERG_S * nu / (BOLTZMANN_ERG_K * T)) * GAUNT_FACTOR
end

"""
    shock_bremsstrahlung_luminosity_nu(p, t, nu) -> erg/s/Hz

Thermal bremsstrahlung spectral luminosity from the hot, shock-heated
downstream layer (NOT the cool pre-shock shell -- bremsstrahlung needs
the ~1e7-1e9 K post-shock gas, which sits in a thin layer of thickness
[`downstream_thickness`](@ref) just behind the shock, at
[`postshock_density`](@ref)/[`shock_temperature`](@ref)):

    n_e = n_i = postshock_density(p,t) / m_p            (fully ionized H, Z=1)
    V   = 4 pi Rs^2 * downstream_thickness(p,t)           (thin-shell volume)
    L_nu(nu) = bremsstrahlung_emissivity_nu(n_e, n_i, Tsh, nu) * V

This is the "hard X-ray" channel flagged in the project's own messenger
decoder table as `shock_temperature` being computed but not yet turned
into an actual luminosity -- this closes that gap.
"""
function shock_bremsstrahlung_luminosity_nu(p::ShockModelParams, t::Real, nu::Real)
    n = postshock_density(p, t) / PROTON_MASS_G
    V = 4 * pi * shock_radius(p, t)^2 * downstream_thickness(p, t)
    Tsh = shock_temperature(p, t)
    return bremsstrahlung_emissivity_nu(n, n, Tsh, nu) * V
end

"""
    shock_bremsstrahlung_luminosity_ev(p, t, E_ev) -> erg/s/eV

Same as [`shock_bremsstrahlung_luminosity_nu`](@ref), converted to
spectral luminosity per unit photon energy (eV) via `nu = E/h`:

    L_E(E) = L_nu(E/h) / h * EV_TO_ERG

the natural unit for compositing this channel into a spectrum alongside
[`QuiescentContinuum.spectral_luminosity_ev`](@ref).
"""
function shock_bremsstrahlung_luminosity_ev(p::ShockModelParams, t::Real, E_ev::Real)
    E_erg = E_ev * EV_TO_ERG
    nu = E_erg / H_ERG_S
    L_nu = shock_bremsstrahlung_luminosity_nu(p, t, nu)
    return L_nu / H_ERG_S * EV_TO_ERG
end

end # module ShockAcceleration
