module NovaMessengers

include("MesaIO.jl")
using .MesaIO
export MesaIO
export read_header, read_history, read_profile, read_profiles_index
export MesaRun, history, profiles_index, profile

include("NuclearDecay.jl")
using .NuclearDecay
export NuclearDecay
export DecayIsotope, MESSENGER_ISOTOPES, EXTENDED_ISOTOPES, decay_constant, mass_number

include("MessengerProduction.jl")
using .MessengerProduction
export MessengerProduction
export isotope_number, decay_rate, positron_rate, neutrino_energy_loss_rate,
       annihilation_photon_rate, DECAY_REACTIONS, zone_decay_rate,
       zone_annihilation_photon_rate, reaction_decay_rate,
       FORMATION_REACTIONS, zone_formation_rate, formation_rate

include("Transport.jl")
using .Transport
export Transport
export ELECTRON_REST_ENERGY_MEV, escape_probability_neutrino, klein_nishina_factor,
       compton_opacity, optical_depth_gamma, escape_probability_gamma

include("ReactionEnergetics.jl")
using .ReactionEnergetics
export ReactionEnergetics
export REACTION_Q_VALUES_MEV, MEV_TO_ERG, TRACKED_REACTIONS,
       zone_reaction_energy_rate, reaction_energy_rate, zone_energy_breakdown,
       energy_breakdown_timeline

include("SignalSynthesis.jl")
using .SignalSynthesis
export SignalSynthesis
export neutrino_lightcurve, neutrino_energy_lightcurve, gamma_lightcurve,
       gamma_energy_lightcurve

include("ShockAcceleration.jl")
using .ShockAcceleration
export ShockAcceleration
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

include("TrajectoryPostProcessing.jl")
using .TrajectoryPostProcessing
export TrajectoryPostProcessing
export LIVE_NET_ISOTOPES, peak_temperature_zone, zone_trajectory,
       write_trajectory_file, zone_initial_abundances,
       write_initial_abundance_file, pristine_profile_number, postprocess_trajectory

include("ExtendedMessengers.jl")
using .ExtendedMessengers
export ExtendedMessengers
export read_mass_fraction_history, extended_decay_rate_per_gram,
       extended_positron_rate_per_gram, extended_gamma_rate_per_gram,
       freezeout_gamma_lightcurve, freezeout_gamma_energy_lightcurve

include("QuiescentContinuum.jl")
using .QuiescentContinuum
export QuiescentContinuum
export PlanckSource, wd_quiescent_source, spectral_luminosity_ev, photon_energy_grid_ev

include("SpectralEvolution.jl")
using .SpectralEvolution
export SpectralEvolution
export LineChannel, ContinuumChannel, SpectrumSnapshot, composite_spectrum, spectrum_timeline

end # module NovaMessengers
