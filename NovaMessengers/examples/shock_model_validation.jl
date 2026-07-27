using NovaMessengers
using Printf

# Validates ShockAcceleration against numbers explicitly quoted in
# Diesing & Metzger (2026), references/model_shock_classical_novae.pdf,
# for their fiducial nova (Table 1): Menv=1e-4 Msun, tau=20 d, vf=6000 km/s,
# fOmega=0.3, xiCR=0.03, xiB=0.01, kappa=0.1.

p = ShockModelParams()  # all defaults = the paper's fiducial values
tpk = 1.1 * p.tau        # Eq (9): shock power peaks at t_pk ~= 1.1*tau

println("== Check 1: peak shock power (Eq 9) ==")
println("paper's own normalized formula: 2.6e38 * Menv,-4 * (vf,8/2)^2 * tau20^-1 erg/s")
println("(NOTE: earlier version of this check wrongly used vf,8^2 instead of (vf,8/2)^2 --")
println(" the paper normalizes by v_sh,8 = vf,8/2, not vf,8 itself. Fixed here.)")
expected_9 = 2.6e38 * 1.0 * (6.0 / 2)^2 / 1.0
@printf("paper's formula, evaluated at fiducial params: %.3e erg/s\n", expected_9)
@printf("computed:                                       Lsh(t_pk) = %.3e erg/s\n\n", shock_power(p, tpk))

println("== Check 2: peak gamma-ray luminosity, calorimetric limit (Eq 24) ==")
fOmega_xiCR_kappa = p.fOmega * p.xiCR * p.kappa
expected_24 = 2.6e35 * (fOmega_xiCR_kappa / 1e-3) * 1.0 * (6.0 / 2)^2 / 1.0
@printf("paper's Eq (24) formula, evaluated at fiducial params: %.3e erg/s\n", expected_24)
@printf("computed (calorimetric approx):                        %.3e erg/s\n",
    gamma_ray_luminosity_calorimetric(p, tpk))
println()

println("== Check 3: full ODE-solved Lgamma agrees with the calorimetric approximation near peak ==")
println("(paper states the fiducial nova is safely in the calorimetric limit at t_pk, Eq 23)")
t_grid = collect(range(1e-3 * p.tau, 4 * p.tau; length=20000))
ECR = cosmic_ray_energy(p, t_grid)
k_pk = argmin(abs.(t_grid .- tpk))
Lgamma_ode = gamma_ray_luminosity(p, t_grid[k_pk], ECR[k_pk])
Lgamma_calorimetric = gamma_ray_luminosity_calorimetric(p, t_grid[k_pk])
@printf("ODE-solved Lgamma(t_pk):    %.3e erg/s\n", Lgamma_ode)
@printf("calorimetric approx:        %.3e erg/s\n", Lgamma_calorimetric)
@printf("ratio (should be close to 1): %.4f\n\n", Lgamma_ode / Lgamma_calorimetric)

println("== Check 4: max proton energy (Eq 32 worked example) ==")
println("NOTE: the paper's Table 1 lists fX=5e-5 as 'fiducial', but this specific")
println("worked numerical example (Eq 32, '~12 GeV') uses fX=1e-4 (fX,-4=1) --")
println("using fX=1e-4 here to match THIS quoted number, not Table 1's fX.")
p_emax = ShockModelParams(fX=1e-4, xiB=0.01, fOmega=0.3)
@printf("paper quotes: Emax ~= 12 GeV at t ~ t_pk\n")
@printf("computed (rL/B chain):    %.3f GeV\n", max_proton_energy_gev(p_emax, tpk))
@printf("computed (closed form):   %.3f GeV\n", max_proton_energy_gev_closedform(p_emax, tpk))
@printf("internal agreement (should match to float precision): %.10f\n\n",
    max_proton_energy_gev(p_emax, tpk) / max_proton_energy_gev_closedform(p_emax, tpk))

println("== Check 5: Emax grows with time (paper: ~100 GeV at t_pk -> >=10 TeV over a couple months) ==")
for t_over_tau in (1.1, 2.0, 3.0, 5.0)
    p_fig2 = ShockModelParams(fX=5e-5, xiB=0.01, fOmega=0.3)  # Figure 2's actual fX
    Emax = max_proton_energy_gev(p_fig2, t_over_tau * p.tau)
    @printf("t = %.1f*tau (%.0f d): Emax = %.3e GeV\n", t_over_tau, t_over_tau * 20, Emax)
end

println("\n== Diagnostic: intermediate quantities feeding Emax, at t_pk, fX=1e-4 case ==")
println("(pinpointing the still-unresolved Check 4 discrepancy against the paper's raw vs.")
println(" normalized-closed-form Delta_rad/Rs -- printed precisely here rather than by hand)")
vsh = shock_velocity(p_emax, tpk)
vw = wind_velocity(p_emax, tpk)
Rs = shock_radius(p_emax, tpk)
Mdw = wind_mass_loss_rate(p_emax, tpk)
Tsh = shock_temperature(p_emax, tpk)
@printf("vsh(t_pk)      = %.4e cm/s  (vsh,8 = %.4f)\n", vsh, vsh / 1e8)
@printf("vw(t_pk)       = %.4e cm/s\n", vw)
@printf("Rs(t_pk)       = %.4e cm\n", Rs)
@printf("Mdot_w(t_pk)   = %.4e g/s\n", Mdw)
@printf("Tsh(t_pk)      = %.4e K   (paper's rounded form: 2e7*vsh,8^2 = %.4e K)\n",
    Tsh, 2e7 * (vsh / 1e8)^2)

delta_rad_raw = radiative_cooling_length_ratio(p_emax, tpk)
Menv_m4 = p_emax.Menv / (1e-4 * ShockAcceleration.MSUN_G)
tau20 = p_emax.tau / (20 * ShockAcceleration.DAY_S)
vsh8 = vsh / 1e8
x = tpk / p_emax.tau
delta_rad_normalized = 5.29e-2 * (1 / Menv_m4) * tau20^2 * vsh8^4 * x * exp(x)
@printf("\nDelta_rad/Rs (raw formula, Eq 26 middle form):        %.4e\n", delta_rad_raw)
@printf("Delta_rad/Rs (paper's normalized closed form, Eq 26): %.4e\n", delta_rad_normalized)
@printf("ratio (should be close to 1 if both readings agree):   %.4f\n\n",
    delta_rad_raw / delta_rad_normalized)

Ddown = downstream_thickness(p_emax, tpk)
B = postshock_field(p_emax, tpk)
@printf("Delta_down(t_pk) = %.4e cm  (Delta_down/Rs = %.4e)\n", Ddown, Ddown / Rs)
@printf("B_down(t_pk)     = %.4e G\n", B)
@printf("Emax(t_pk)       = %.4f GeV\n", max_proton_energy_gev(p_emax, tpk))

println("\n== Check 6: exact vs. asymptotic (t>>tau) formulas converge as t/tau grows ==")
println("(Eq 26's normalized closed form assumes v_sh has saturated to vf/2 -- true only")
println(" for t >> tau. This checks the exact and asymptotic Delta_rad/Rs agree once we're")
println(" actually deep in that regime, confirming t_pk's earlier mismatch was a regime-of-")
println(" validity issue, not a formula error.)")
vsh_inf = p_emax.vf / 2
@printf("%8s  %12s  %10s  %14s  %14s  %8s\n",
    "t/tau", "vsh/vsh_inf", "raw D/Rs", "asymptotic D/Rs", "Emax [GeV]", "")
for t_over_tau in (1.1, 3.0, 10.0, 30.0, 100.0)
    t = t_over_tau * p_emax.tau
    vsh_t = shock_velocity(p_emax, t)
    raw = radiative_cooling_length_ratio(p_emax, t)
    asym = 5.29e-2 * (1 / Menv_m4) * tau20^2 * (vsh_t / 1e8)^4 * t_over_tau * exp(t_over_tau)
    @printf("%8.1f  %12.4f  %10.4e  %14.4e  %14.4f  ratio=%.4f\n",
        t_over_tau, vsh_t / vsh_inf, raw, asym, max_proton_energy_gev(p_emax, t), raw / asym)
end

println("\n== Check 7: paper states Emax reaches 10 TeV at t ~ 2.5*tau (Figure 2 params, page 9) ==")
p_fig2b = ShockModelParams(fX=5e-5, xiB=0.01, fOmega=0.3)
Emax_25tau = max_proton_energy_gev(p_fig2b, 2.5 * p_fig2b.tau)
@printf("paper: Emax ~ 10 TeV (1e4 GeV) at t ~ 2.5*tau\n")
@printf("computed: Emax(2.5*tau) = %.3e GeV\n", Emax_25tau)
