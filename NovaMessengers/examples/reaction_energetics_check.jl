using NovaMessengers
using Printf

const WORK_DIR = joinpath(@__DIR__, "..", "..", "mesa_work", "wd_nova_burst_co")

run = MesaRun(WORK_DIR)
pidx = profiles_index(run)
last_profile_number = pidx.profile_number[argmax(pidx.model_number)]

println("== zone_energy_breakdown on the final profile (profile$(last_profile_number)) ==")
breakdown = zone_energy_breakdown(profile(run, last_profile_number).data)
peak_zone = argmax(breakdown.total)
@printf("peak zone: %d, total = %.3e erg/s\n", peak_zone, breakdown.total[peak_zone])
println("dominant reactions in that zone (top 5 by |rate|):")
row = breakdown[peak_zone, :]
ranked = sort([(r, row[r]) for r in TRACKED_REACTIONS]; by=p -> -abs(p[2]))
for (reaction, rate) in ranked[1:5]
    @printf("  %-16s %.3e erg/s\n", reaction, rate)
end

println("\n== reaction_energy_rate time series (whole star, profile cadence) ==")
for reaction in (:r_o14_ap_f17, :r_f17_pg_ne18, :r_n13_wk_c13)
    age, rate = reaction_energy_rate(run, reaction)
    peak = argmax(rate)
    @printf("%-16s peak %.3e erg/s at star_age = %.4f yr\n",
        reaction, rate[peak], age[peak] / (365.25 * 24 * 3600))
end
