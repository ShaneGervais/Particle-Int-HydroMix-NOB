using NovaMessengers
using Printf

# Bridges the CO WD baseline into ReacNetJl: locates the peak-temperature
# zone, extracts its T9(t)/rho(t) trajectory and pristine (pre-burn)
# composition, and evolves it through ReacNetJl's full H-Ca network --
# recovering isotopes (22Na, 26Al, 7Be, ...) beyond the 7 tracked directly
# by MESA's own live net. See TrajectoryPostProcessing.jl's docstrings for
# the physics behind each non-obvious choice here (pristine-profile
# tie-break, dt_max rescaling, T9-floor clamping).

const WORK_DIR = joinpath(@__DIR__, "..", "..", "mesa_work", "wd_nova_burst_co")
const OUTPUT_DIR = joinpath(@__DIR__, "reacnetjl_out")

run = MesaRun(WORK_DIR)

peak = peak_temperature_zone(run)
@printf("Peak-temperature zone: profile=%d  mass=%.6f Msun  T=%.4e K\n",
    peak.profile_number, peak.mass_coordinate, peak.temperature_K)

pn = pristine_profile_number(run, peak.mass_coordinate)
println("Pristine (pre-burn) composition taken from profile ", pn)

result = postprocess_trajectory(run; output_dir=OUTPUT_DIR)

@printf("\nTrajectory handed to ReacNetJl: %d points, %.4e s -> %.4e s (%.2f days)\n",
    length(result.trajectory.time_s), result.trajectory.time_s[1], result.trajectory.time_s[end],
    result.trajectory.time_s[end] / 86400)
println("Output files written under: ", OUTPUT_DIR)
println("  ", result.result.output_files)

println("\n== Final mass fractions (top 20) ==")
fmf = result.result.final_mass_fractions
for (iso, X) in sort(collect(fmf); by=kv -> -kv[2])[1:min(20, length(fmf))]
    @printf("%-6s  %.6e\n", iso, X)
end

println("\n== Isotopes beyond MESA's own live net (the actual point of this bridge) ==")
for iso in ("be7", "na22", "al26", "al*6")
    @printf("%-6s  %.6e\n", iso, get(fmf, iso, 0.0))
end
