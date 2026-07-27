using NovaMessengers
using CairoMakie
using Printf

# The combined observable gamma-ray spectrum from both messenger-production
# channels this package tracks:
#   1. Radioactive decay -> positron annihilation -> the 511 keV line,
#      computed from mesa_work/wd_nova_burst_co's REAL simulated burst.
#   2. Shock-accelerated protons -> pion decay -> a continuum with an
#      E^-2 power law and an exponential cutoff near Emax, computed from
#      Diesing & Metzger's STANDALONE, literature-parameterized model
#      (NOT tied to our specific run -- see ShockAcceleration.jl's
#      ShockModelParams docstring for why).
# These are two different scenarios (a specific simulated CO WD nova vs. a
# generic fiducial nova), shown together to compare channel character and
# rough scale, not as a single self-consistent prediction for one event.

const WORK_DIR = joinpath(@__DIR__, "..", "..", "mesa_work", "wd_nova_burst_co")
const MEV_TO_GEV = 1e-3

# --- Channel 1: decay line, from the real run -------------------------------

run = MesaRun(WORK_DIR)
line_energy_mev = 0.511
g = gamma_energy_lightcurve(run, line_energy_mev)  # (age, rate) erg/s vs star_age
peak_idx = argmax(g.rate)
line_power_erg_s = g.rate[peak_idx]
line_age_yr = g.age[peak_idx] / (365.25 * 24 * 3600)
@printf("Decay-line channel: peak 511 keV luminosity = %.3e erg/s at star_age = %.4f yr\n",
    line_power_erg_s, line_age_yr)
# CAVEAT (see examples/isotope_abundances.jl for the full derivation): at this
# peak, the effective escape fraction is ~5.7e-12 -- this run's window ends
# before the envelope expands enough to actually become transparent to 511
# keV photons, so this is the least-opaque moment available in the data, not
# a genuinely observable one. The isotope decay rates driving it are real;
# this specific line luminosity is illustrative of the method, not a
# physically meaningful predicted peak.

# --- Channel 2: shock continuum, standalone fiducial model ------------------

p = ShockModelParams()
t_now = 2.0 * p.tau  # where the bolometric consistency check was tightest (~1.24x)

injection_grid = collect(range(1e-3 * p.tau, t_now; length=1500))
Emax_now = max_proton_energy_gev(p, t_now)
# lower bound ~10 MeV: below ~kappa*(proton rest energy) = 0.1*0.938 GeV ~= 94 MeV,
# the underlying proton energies (E/kappa) go sub-relativistic and the E^-2
# injection spectrum is no longer physically meaningful -- keep some margin
# above that rather than extending the computed grid down to the line's own
# 0.511 MeV (the line is plotted separately and doesn't need the continuum's
# grid to reach it).
E_grid_gev = exp.(range(log(1e-2), log(50 * Emax_now); length=1500))

dNdE_intrinsic = gamma_ray_spectrum(p, E_grid_gev, t_now, injection_grid)  # photons/(GeV s)

# apply absorption: what actually escapes to a distant observer
transmission = [exp(-bethe_heitler_optical_depth(p, t_now, E) - gamma_gamma_optical_depth(p, t_now, E))
                for E in E_grid_gev]
dNdE_observed = dNdE_intrinsic .* transmission

E2dNdE_intrinsic = E_grid_gev .^ 2 .* dNdE_intrinsic .* ShockAcceleration.GEV_TO_ERG  # erg/s
E2dNdE_observed = E_grid_gev .^ 2 .* dNdE_observed .* ShockAcceleration.GEV_TO_ERG

@printf("Shock channel: Emax(t=2*tau) = %.3e GeV; peak E^2 dN/dE (intrinsic) = %.3e erg/s\n",
    Emax_now, maximum(E2dNdE_intrinsic))

# --- Combined plot -----------------------------------------------------------

E_grid_mev = E_grid_gev ./ MEV_TO_GEV

set_theme!(fontsize=11)
fig = Figure(size=(750, 500))
ax = Axis(fig[1, 1]; xscale=log10, yscale=log10,
    xlabel="Photon Energy [MeV]", ylabel="E² dN/dE [erg/s]",
    title="Predicted gamma-ray spectrum: decay line (this simulated nova) + shock continuum (fiducial model)")

valid = E2dNdE_intrinsic .> 0
lines!(ax, E_grid_mev[valid], E2dNdE_intrinsic[valid]; linestyle=:dash, color=:steelblue,
    label="shock continuum, intrinsic (t=2τ)")
valid_obs = E2dNdE_observed .> 0
lines!(ax, E_grid_mev[valid_obs], E2dNdE_observed[valid_obs]; color=:steelblue, linewidth=2,
    label="shock continuum, observed (absorption applied)")

# the line is a delta function in energy -- represent its own luminosity directly as a marker,
# not as an E^2 dN/dE density (which isn't well-posed for a discrete line the same way)
scatter!(ax, [line_energy_mev], [line_power_erg_s]; color=:firebrick3, markersize=14,
    marker=:diamond, label="511 keV decay line (peak, this run)")

axislegend(ax; position=:lb, labelsize=9)

outdir = joinpath(@__DIR__, "plt_out")
mkpath(outdir)
outpath = joinpath(outdir, "final_spectrum.png")
save(outpath, fig)
println("Wrote $outpath")
