using NovaMessengers
using Printf

# What isotope abundances actually produced the 511 keV signal in
# final_spectrum.jl, and how much of that signal's origin is actually
# recoverable from the line alone (spoiler: only the *total*, not the
# per-isotope breakdown -- all 7 tracked isotopes emit the same 511 keV
# photon, so the line by itself can't distinguish which one decayed).

const WORK_DIR = joinpath(@__DIR__, "..", "..", "mesa_work", "wd_nova_burst_co")
const YEAR_S = 365.25 * 24 * 3600

run = MesaRun(WORK_DIR)

g = gamma_energy_lightcurve(run, 0.511)
peak_idx = argmax(g.rate)
peak_age_yr = g.age[peak_idx] / YEAR_S
L511_observed = g.rate[peak_idx]

hist = history(run).data
j = argmin(abs.(hist.star_age .- peak_age_yr))
msun = history(run).header.msun

println("== Isotope inventory at the 511 keV emission peak (star_age = $(round(peak_age_yr, digits=4)) yr) ==\n")
@printf("%-6s %12s %12s %14s %16s %10s\n", "iso", "mass [g]", "mass frac", "N [nuclei]", "decay rate [/s]", "% of total")

total_decay_rate = 0.0
rates = Dict{Symbol,Float64}()
for iso in keys(MESSENGER_ISOTOPES)
    N = isotope_number(run, iso)[j]
    A = mass_number(iso)
    mass_g = N * A * MessengerProduction.AMU_G
    lambda = decay_constant(MESSENGER_ISOTOPES[iso])
    rate = lambda * N
    rates[iso] = rate
    global total_decay_rate += rate
end
for iso in sort(collect(keys(rates)); by=i -> -rates[i])
    N = isotope_number(run, iso)[j]
    A = mass_number(iso)
    mass_g = N * A * MessengerProduction.AMU_G
    mass_frac = mass_g / (msun * (hist.star_mass[j]))
    pct = 100 * rates[iso] / total_decay_rate
    @printf("%-6s %12.4e %12.4e %14.4e %16.4e %9.2f%%\n", iso, mass_g, mass_frac, N, rates[iso], pct)
end

if "total_mass_h1" in string.(propertynames(hist))
    @printf("\n%-6s %12.4e %12.4e\n", "h1", hist.total_mass_h1[j] * msun, hist.total_mass_h1[j] / hist.star_mass[j])
end
if "total_mass_he4" in string.(propertynames(hist))
    @printf("%-6s %12.4e %12.4e\n", "he4", hist.total_mass_he4[j] * msun, hist.total_mass_he4[j] / hist.star_mass[j])
end

println("\n== What the 511 keV line alone can and can't tell an observer ==")
@printf("Observed L(511 keV) at this snapshot:            %.4e erg/s\n", L511_observed)
implied_decay_rate = L511_observed / (2 * 0.511 * 1.602176634e-6)  # 2 photons/positron, MeV_to_erg
@printf("Implied decay rate IF fully transparent:         %.4e /s\n", implied_decay_rate)
@printf("Actual Sum(lambda_i * N_i) produced in the star: %.4e /s\n", total_decay_rate)
@printf("ratio = effective escape fraction at this snapshot: %.4f\n\n", implied_decay_rate / total_decay_rate)

println("That ratio is NOT expected to be 1 -- gamma_energy_lightcurve already applies Phase 3's")
println("Compton escape_probability_gamma (zone-by-zone) before summing, so naively inverting the")
println("observed line assuming full transparency recovers less than the true decay rate whenever")
println("photons are produced deep enough to be partially absorbed in transit. This IS the central")
println("point: the messenger's observed signal is not identical to what was produced -- transport")
println("reshapes it, and the size of this ratio is a direct, quantitative measure of how much.")
println()
println("Separately: even a perfectly transparent 511 keV line could only ever recover the TOTAL")
println("decay rate above, never the per-isotope breakdown printed above it -- every one of these")
println("7 isotopes emits the identical annihilation photon. Telling them apart needs different")
println("information: their distinct half-lives shape the light curve's TIMING (a fast burst-")
println("synchronous rise/fall points to short-lived o14/f17/ne18 dominating; a slower decay tail")
println("points to longer-lived f18/n13), or -- in principle -- the neutrino channel, since each")
println("isotope's neutrinos carry a distinct mean energy (q_neu_mev) even though all 7 also share")
println("the same free-streaming escape=1 transport.")
