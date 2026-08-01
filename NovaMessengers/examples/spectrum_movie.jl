using NovaMessengers
using CairoMakie
using Printf

# End-to-end multi-messenger spectrum movie: composes every channel this
# project can currently compute (neutrino lines, gamma decay lines,
# freeze-out lines if a ReacNetJl run is available, WD+companion
# continuum, shock bremsstrahlung + GeV continuum if shock params are
# calibrated) into one flux-vs-energy animation across star_age.
#
# Degrades gracefully: channels that need inputs you haven't generated
# yet (a ReacNetJl mass_fractions.csv, in particular) are simply left
# out rather than erroring, so this runs with just the MESA output too.

const WORK_DIR = joinpath(@__DIR__, "..", "..", "mesa_work", "wd_nova_burst_co")
const REACNETJL_OUT = joinpath(@__DIR__, "reacnetjl_out")
const OUTPUT_DIR = joinpath(@__DIR__, "plt_out")
const YEAR_S = 365.25 * 24 * 3600
const GEV_TO_EV = 1.0e9
const MEV_TO_EV = 1.0e6
const EV_TO_ERG = 1.602176634e-12

mkpath(OUTPUT_DIR)
run = MesaRun(WORK_DIR)

# --- Lines: the 7 MESA-tracked hot-CNO isotopes ------------------------

lines = LineChannel[]
for (iso, info) in MESSENGER_ISOTOPES
    r = decay_rate(run, iso)
    push!(lines, LineChannel("$(iso) neutrino", info.q_neu_mev * MEV_TO_EV, r.age, r.rate))
end

g511 = gamma_lightcurve(run, 0.511)
push!(lines, LineChannel("511 keV (e+e- annihilation)", 0.511 * MEV_TO_EV, g511.age, g511.rate))

# --- Lines: freeze-out channels, only if a ReacNetJl run is available ---

peak = peak_temperature_zone(run)
mfh_path = joinpath(REACNETJL_OUT, "mass_fractions.csv")
if isfile(mfh_path)
    println("Found ReacNetJl output at ", mfh_path, " -- adding na22/al26 freeze-out lines.")
    mfh = read_mass_fraction_history(mfh_path)
    pn = pristine_profile_number(run, peak.mass_coordinate)
    hist = history(run).data
    pidx = profiles_index(run)
    pmn = pidx[findfirst(==(pn), pidx.profile_number), :model_number]
    j = findfirst(==(pmn), hist.model_number)
    t0_star_age_s = hist.star_age[j] * YEAR_S

    # Explicit, documented approximation (see ExtendedMessengers.jl):
    # assume the envelope mass above the ignition zone shares this
    # trajectory's history.
    total_mass_g = max(hist.star_mass[end] - peak.mass_coordinate, 1e-6) * 1.9884098706980504e33

    for iso in (:na22, :al26)
        g = freezeout_gamma_lightcurve(run, mfh, iso, peak.mass_coordinate, total_mass_g; t0_star_age_s=t0_star_age_s)
        line_energy_ev = EXTENDED_ISOTOPES[iso].gamma_line_mev * MEV_TO_EV
        push!(lines, LineChannel("$(iso) $(round(EXTENDED_ISOTOPES[iso].gamma_line_mev, digits=3)) MeV", line_energy_ev, g.age, g.rate))
    end
else
    println("No ReacNetJl output found at ", mfh_path, " -- skipping freeze-out lines.")
    println("Run `julia --project=. examples/trajectory_postprocessing.jl` first to add them.")
end

# --- Continua: WD (+ placeholder companion) and shock channels ---------

# Time-dependent WD photosphere, not a single frozen snapshot: MESA's own
# log_Teff/radius_cm history already carries the WD's surface through
# quiescence, the TNR flash/envelope inflation, and the post-burst
# recontraction back to a hot, compact, supersoft-X-ray-like state (this
# run's own data reaches Teff > 6e5 K, Wien peak ~140 eV) -- see
# wd_photosphere_at's own docstring for the physical picture. The
# companion stays static: nothing in this run models binary irradiation.
companion = PlanckSource(4000.0, 3.0e10)  # generic cool-dwarf placeholder -- adjust to your system
continua = ContinuumChannel[
    ContinuumChannel("WD photosphere (flash/decline/SSS) + companion",
        (t, E) -> spectral_luminosity_ev([wd_photosphere_at(run, t), companion], E)),
]

shock_p = calibrate_shock_params(run; ignition_mass_coordinate=peak.mass_coordinate)

# ShockModelParams' own `t` is time SINCE SHOCK/WIND ONSET (scale: `tau`,
# ~days), not absolute star_age (scale: this run's ~1e11 s, dominated by
# the multi-millennium pre-outburst accretion tail) -- confirmed the hard
# way: passing star_age directly sent `max_proton_energy_gev` to
# `t/tau ~ 5e4` for most of the timeline, catastrophically past the
# "couple months" validity window documented on that function, and
# crashed proton_spectrum's log-spaced integration outright once
# `50*Emax_gev` went negative. Rebase to the burst-peak profile's own
# star_age (a reasonable proxy for "ejection begins") and treat anything
# before that as pre-shock (zero signal), matching how the freeze-out/
# gamma channels already only turn on after their own relevant epoch.
hist_r = history(run).data
pidx_r = profiles_index(run)
peak_model_number = pidx_r[findfirst(==(peak.profile_number), pidx_r.profile_number), :model_number]
t_shock_onset_s = hist_r.star_age[findfirst(==(peak_model_number), hist_r.model_number)] * YEAR_S

push!(continua, ContinuumChannel("Shock bremsstrahlung", (t, E) -> begin
    t_sh = t - t_shock_onset_s
    t_sh <= 0 ? 0.0 : shock_bremsstrahlung_luminosity_ev(shock_p, t_sh, E)
end))

function shock_gamma_continuum_ev(t::Real, E_ev::Real)
    t_sh = t - t_shock_onset_s
    t_sh <= 0 && return 0.0
    injection_grid = range(1.0e-3 * shock_p.tau, t_sh; length=50)
    E_gev = E_ev / GEV_TO_EV
    dNdE = gamma_ray_spectrum(shock_p, [E_gev], t_sh, injection_grid)[1]  # photons/s/GeV
    return (E_ev * EV_TO_ERG) * dNdE / GEV_TO_EV                          # erg/s/eV
end
push!(continua, ContinuumChannel("Shock GeV gamma-rays (pion decay)", shock_gamma_continuum_ev))

# --- Build the timeline --------------------------------------------------

# Log-spaced around t_shock_onset_s (burst peak), not uniform-by-index
# over the saved history rows. MESA's saved ages are extremely densely
# packed near the burst (the same clustering pristine_profile_number's
# tie-break has to account for) -- uniform-by-index sampling spent
# nearly the whole frame budget inside that dense cluster and gave the
# multi-millennium quiescent baseline almost no screen time. Composing
# the spectrum at an arbitrary time doesn't require a saved row to exist
# there (LineChannel rates are interpolated, ContinuumChannel closures
# are evaluated directly -- see SpectralEvolution.jl), so nothing stops
# sampling on whatever time grid best represents the story: like the
# real observational picture for a nova's multi-wavelength evolution
# (e.g. Chomiuk, Metzger & Shen 2021, "New Insights into Classical
# Novae", Fig 1), that's log(time since eruption), not linear absolute
# time. This grid is symmetric-log around t_shock_onset_s: dense near
# peak (captures the TNR flash and the day-to-week shock/outflow
# activity), sparse far from it on both sides (the multi-millennium
# quiescent past and the ~1.5 yr post-burst decline this run covers).
hist = history(run).data
age_min_s = minimum(hist.star_age) * YEAR_S
age_max_s = maximum(hist.star_age) * YEAR_S
n_frames = 200
n_pre = n_frames ÷ 2
n_post = n_frames - n_pre
pre_ages = t_shock_onset_s .- exp.(range(log(t_shock_onset_s - age_min_s), log(1.0); length=n_pre))
post_ages = t_shock_onset_s .+ exp.(range(log(1.0), log(age_max_s - t_shock_onset_s); length=n_post))
frame_ages = vcat(pre_ages, post_ages)

energy_grid_ev = photon_energy_grid_ev(1.0e-2, 1.0e10; n=150)  # 0.01 eV (radio-ish) to 10 GeV

println("Building spectrum timeline: ", n_frames, " frames, ", length(lines), " lines, ", length(continua), " continua...")
timeline = spectrum_timeline(frame_ages, lines, continua, energy_grid_ev)

# --- Render the movie ------------------------------------------------------

fig = Figure(size=(1000, 650))
ax = Axis(fig[1, 1], xscale=log10, yscale=log10,
    xlabel="Photon/neutrino energy [eV]", ylabel="Spectral luminosity [erg/s/eV]  (lines: rate x 1 eV placeholder height)",
    title="Nova multi-messenger spectrum")
xlims!(ax, 1.0e-2, 1.0e10)
# 1e46, not 1e40: the WD photosphere briefly exceeds the old ceiling
# during/just after the burst -- consistent with the super-Eddington
# wind episodes MESA's own log reports for this run (L/L_Edd up to
# ~3.45), not a bug, but the old fixed ceiling was clipping real signal.
ylims!(ax, 1.0e10, 1.0e46)

# Observable-backed plot data: this is the pattern `record` actually
# needs to redraw each frame (directly reassigning a plot object's [i]
# index, e.g. `plot_obj[2] = new_ys`, does NOT trigger a redraw and in
# some Makie versions errors outright -- confirmed by testing both
# before settling on this).
continuum_L_obs = [Observable(fill(1.0e10, length(energy_grid_ev))) for _ in continua]
for (i, cc) in enumerate(continua)
    lines!(ax, energy_grid_ev, continuum_L_obs[i]; label=cc.label)
end
line_E_obs = Observable(Float64[])
line_R_obs = Observable(Float64[])
scatter!(ax, line_E_obs, line_R_obs; markersize=12, color=:black, label="decay/neutrino lines")
axislegend(ax; position=:lt, labelsize=10)
time_label = Label(fig[0, 1], "")

"""
Human-scaled relative-time string: seconds near `dt_s=0`, growing to
hours/days/years as `|dt_s|` grows. A single fixed unit (e.g. always
years) either loses all resolution during the seconds-to-days burst
dynamics or is unreadably large during the multi-year quiescent
baseline -- this is the same "pick units matched to the scale you're
actually looking at" logic as `TrajectoryPostProcessing`'s log-spaced
time sampling, applied to display rather than solver stepping.
"""
function format_relative_time(dt_s::Real)
    adt = abs(dt_s)
    if adt < 3600
        return @sprintf("%+.4g s", dt_s)
    elseif adt < 86400
        return @sprintf("%+.4g hr", dt_s / 3600)
    elseif adt < YEAR_S
        return @sprintf("%+.4g days", dt_s / 86400)
    else
        return @sprintf("%+.4g yr", dt_s / YEAR_S)
    end
end

record(fig, joinpath(OUTPUT_DIR, "spectrum_movie.mp4"), timeline; framerate=12) do snap
    for i in eachindex(continua)
        continuum_L_obs[i][] = max.(snap.continuum_L_ev[i, :], 1e-30)
    end
    nonzero = snap.line_rates .> 0
    line_E_obs[] = snap.line_energies_ev[nonzero]
    line_R_obs[] = max.(snap.line_rates[nonzero], 1e10)
    dt = snap.age_s - t_shock_onset_s
    time_label.text[] = @sprintf("star_age = %.10g yr   (t - t_peak = %s)", snap.age_s / YEAR_S, format_relative_time(dt))
end

println("\nWrote ", joinpath(OUTPUT_DIR, "spectrum_movie.mp4"))
