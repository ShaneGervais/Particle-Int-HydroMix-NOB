module NovaMessengers

include("MesaIO.jl")
using .MesaIO
export MesaIO
export read_header, read_history, read_profile, read_profiles_index
export MesaRun, history, profiles_index, profile

include("NuclearDecay.jl")
using .NuclearDecay
export NuclearDecay
export DecayIsotope, MESSENGER_ISOTOPES, decay_constant, mass_number

include("MessengerProduction.jl")
using .MessengerProduction
export MessengerProduction
export isotope_number, production_rate, positron_rate, neutrino_energy_loss_rate,
       annihilation_photon_rate, ne18_production_rate

include("Transport.jl")
using .Transport
export Transport
export ELECTRON_REST_ENERGY_MEV, escape_probability_neutrino, klein_nishina_factor,
       compton_opacity, optical_depth_gamma, escape_probability_gamma

end # module NovaMessengers
