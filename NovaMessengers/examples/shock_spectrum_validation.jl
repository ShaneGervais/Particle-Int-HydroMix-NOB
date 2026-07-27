using NovaMessengers
using Printf

# Validates the Eq (34)-(36) spectral synthesis against the paper's own
# stated consistency claim: "integrating this [dNgamma/dE] spectrum gives
# a bolometric Lgamma that is within 1% of that calculated by solving
# Equation (19)" (i.e. the ODE-based gamma_ray_luminosity).

p = ShockModelParams()  # fiducial
t_now = 1.1 * p.tau      # t_pk

println("== Check: bolometric Lgamma from integrating the spectrum vs. the ODE solution ==")

# injection epochs: fine grid from well before t_pk up to t_now
injection_grid = collect(range(1e-3 * p.tau, t_now; length=2000))

# energy grid spanning proton rest energy up past the cutoff -- log-spaced,
# since this is a power-law-times-exponential spanning many decades
Emax_now = max_proton_energy_gev(p, t_now)
E_grid = exp.(range(log(1e-2), log(50 * Emax_now); length=2000))  # GeV

dNdE = gamma_ray_spectrum(p, E_grid, t_now, injection_grid)  # photons/GeV

# integrate E*dN/dE dE (log-spaced trapezoidal, same method the module uses internally)
L_from_spectrum_gev_s = 0.0
f_prev = E_grid[1] * dNdE[1]
E_prev = E_grid[1]
for (E, dN) in zip(Iterators.drop(E_grid, 1), Iterators.drop(dNdE, 1))
    f_cur = E * dN
    global L_from_spectrum_gev_s += 0.5 * (f_cur + f_prev) * (E - E_prev)
    global f_prev, E_prev = f_cur, E
end
L_from_spectrum = L_from_spectrum_gev_s * ShockAcceleration.GEV_TO_ERG  # erg/s

# ODE-based bolometric Lgamma at the same time, for comparison
t_grid = collect(range(1e-3 * p.tau, 4 * p.tau; length=20000))
ECR = cosmic_ray_energy(p, t_grid)
k_now = argmin(abs.(t_grid .- t_now))
L_from_ode = gamma_ray_luminosity(p, t_grid[k_now], ECR[k_now])

@printf("Lgamma from integrating dN/dE spectrum: %.4e erg/s\n", L_from_spectrum)
@printf("Lgamma from ODE (Eq 19):                %.4e erg/s\n", L_from_ode)
@printf("ratio (paper claims within 1%%):         %.4f\n\n", L_from_spectrum / L_from_ode)

println("== Spectrum shape: dN/dE at a few energies (should show E^-2 power law, cutoff near Emax) ==")
@printf("Emax(t_now) = %.3f GeV\n", Emax_now)
for E in (0.1, 1.0, 10.0, Emax_now / 2, Emax_now, 5 * Emax_now)
    idx = argmin(abs.(E_grid .- E))
    @printf("E = %8.3f GeV:  dN/dE = %.4e photons/GeV\n", E_grid[idx], dNdE[idx])
end

println("\n== Absorption at t_now: Bethe-Heitler and gamma-gamma optical depths ==")
for Eg in (1.0, 10.0, 100.0, 1000.0, 10000.0)
    tau_bh = bethe_heitler_optical_depth(p, t_now, Eg)
    tau_gg = gamma_gamma_optical_depth(p, t_now, Eg)
    @printf("Eg = %8.1f GeV:  tau_BH = %.4e   tau_gg = %.4e\n", Eg, tau_bh, tau_gg)
end

println("\n== Diagnostic: does the 1.49x mismatch shrink as t_now moves deeper into t>>tau? ==")
println("(Hypothesis: Eq 19's adiabatic-loss term ECR/t implicitly assumes d(ln Rs)/dt = 1/t,")
println(" i.e. Rs ~ t, which only holds once t >> tau -- but this module's spectral synthesis")
println(" uses the *exact* Rs(t) ratio (Lad = Rs(t_now)/Rs(t_i)) at any t. At t_now = t_pk =")
println(" 1.1*tau (not deep in that regime), the two would be expected to disagree; if the")
println(" ratio below converges toward 1 together with t*d(ln Rs)/dt converging toward 1,")
println(" that confirms the ODE's adiabatic term -- not the spectral synthesis -- is the")
println(" approximate one here, same root cause as the earlier Emax/Delta_rad finding.)\n")

function bolometric_ratio_at(p, t_now; n_epoch=1500, n_E=1500)
    injg = collect(range(1e-3 * p.tau, t_now; length=n_epoch))
    Emax_now = max_proton_energy_gev(p, t_now)
    Eg = exp.(range(log(1e-2), log(50 * Emax_now); length=n_E))
    dN = gamma_ray_spectrum(p, Eg, t_now, injg)
    L_spec = 0.0
    fp, Ep = Eg[1] * dN[1], Eg[1]
    for (E, d) in zip(Iterators.drop(Eg, 1), Iterators.drop(dN, 1))
        fc = E * d
        L_spec += 0.5 * (fc + fp) * (E - Ep)
        fp, Ep = fc, E
    end
    L_spec *= ShockAcceleration.GEV_TO_ERG

    tg = collect(range(1e-3 * p.tau, max(4 * p.tau, 1.5 * t_now); length=20000))
    E_cr = cosmic_ray_energy(p, tg)
    k = argmin(abs.(tg .- t_now))
    L_ode = gamma_ray_luminosity(p, tg[k], E_cr[k])
    return L_spec / L_ode
end

function d_ln_Rs_dt_times_t(p, t; eps_frac=1e-4)
    eps = eps_frac * t
    Rs_plus = shock_radius(p, t + eps)
    Rs_minus = shock_radius(p, t - eps)
    dlnRs_dt = (log(Rs_plus) - log(Rs_minus)) / (2 * eps)
    return t * dlnRs_dt
end

@printf("%10s  %18s  %20s\n", "t/tau", "t*d(lnRs)/dt", "spectrum/ODE ratio")
for t_over_tau in (1.1, 2.0, 3.0, 5.0, 8.0)
    tt = t_over_tau * p.tau
    @printf("%10.1f  %18.4f  %20.4f\n", t_over_tau, d_ln_Rs_dt_times_t(p, tt), bolometric_ratio_at(p, tt))
end
