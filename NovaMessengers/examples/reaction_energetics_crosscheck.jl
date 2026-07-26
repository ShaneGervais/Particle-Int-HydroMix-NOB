using NovaMessengers
using Printf

const WORK_DIR = joinpath(@__DIR__, "..", "..", "mesa_work", "wd_nova_burst_co")
run = MesaRun(WORK_DIR)

println("Cross-check: ReactionEnergetics-derived decay rate vs MessengerProduction.decay_rate")
println("(two independent data paths: screened_rate profile columns vs total_mass_<iso> history bookkeeping)\n")

for (iso, reaction) in DECAY_REACTIONS
    age_e, rate_e = reaction_energy_rate(run, reaction)  # erg/s, profile cadence
    q = REACTION_Q_VALUES_MEV[reaction]
    decays_from_energy = rate_e ./ (q * MEV_TO_ERG)        # back out decays/s

    age_h, rate_h = decay_rate(run, iso)  # decays/s, dense history cadence

    decays_from_history = similar(decays_from_energy)
    for (i, a) in enumerate(age_e)
        j = argmin(abs.(age_h .- a))
        decays_from_history[i] = rate_h[j]
    end

    k = argmax(decays_from_energy)
    ratio = decays_from_energy[k] / decays_from_history[k]
    @printf("%-6s peak: screened_rate-derived = %.3e /s, mass-bookkeeping = %.3e /s, ratio = %.3f\n",
        iso, decays_from_energy[k], decays_from_history[k], ratio)
end
