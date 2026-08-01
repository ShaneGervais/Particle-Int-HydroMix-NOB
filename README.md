# Particle-Int-HydroMix-NOB

**What can the particles coming off an exploding star tell us about what's happening inside it?**

This project simulates a classical nova — a thermonuclear explosion on the surface of a white dwarf star — and asks: if you could detect every kind of particle it sends out (light of every wavelength, neutrinos, cosmic rays), what would each one tell you, and when? Different particles come from different depths, different moments, and travel through the exploding star differently before they escape. The goal is to build the whole chain, from the explosion itself to the signal an observer would actually see, and watch how that signal changes over the seconds-to-years the event actually takes.

This README explains the physics in plain language, gives the equations the code actually uses, and shows how the pieces fit together and how to run them.

---

## 1. The physics story

A classical nova happens in a **binary star system**: a white dwarf (the compressed, burned-out core of a dead star) orbiting closely with a normal companion star. Gas slowly leaks off the companion and piles up on the white dwarf's surface. That gas — mostly hydrogen — gets squeezed by the white dwarf's enormous gravity until, after thousands of years of quiet accumulation, it ignites all at once. That's the nova.

Here's what happens, in order, and what each stage sends out:

1. **Quiescence** (thousands of years). Nothing dramatic — just the white dwarf and its companion, each glowing faintly at their own temperature. This is the baseline "before" picture.
2. **Thermonuclear runaway (TNR)** (minutes to hours). The hydrogen layer ignites. It's so deeply buried that light can't escape yet — but **neutrinos** can, instantly and without being altered. They are the one messenger that comes straight from the burning region with no distortion.
3. **The flash and envelope expansion** (hours to days). The energy release blows the envelope up to huge size. As it expands, its surface cools even as its total energy output rises — this is the classical "nova reaching optical maximum."
4. **Shocks** (days to weeks). The star doesn't eject material smoothly — it throws off a slower shell first, then a faster wind that catches up and rams into it. Where they collide, particles get accelerated to enormous energies, producing gamma-rays and X-rays that have nothing to do with the nuclear burning itself.
5. **Decline and freeze-out** (weeks to years). As the ejected gas thins out, it becomes transparent to gamma-ray light from the decay of longer-lived radioactive isotopes made during the burning. The white dwarf's own surface, now exposed, can be hot enough to shine in soft X-rays for a while — the "supersoft X-ray source" phase.
6. **Back to quiescence** (years later). The burning stops, the white dwarf cools, and the system returns to where it started, waiting to do it again.

The project's job is to model the source of each of these signals, how each one is reshaped (or not) on its way out, and stitch them into one evolving picture — literally a movie of the spectrum changing over the nova's whole lifetime.

---

## 2. Three tools, one pipeline

No single tool does all of this, so the project chains three together:

```
   MESA                    ReacNetJl                  NovaMessengers
(hydrodynamics +     -->  (nuclear network      -->   (turns all of the
 basic nuclear            post-processing:            above into particle
 network)                 recovers isotopes           messengers and an
                           MESA's own small             observable spectrum)
                           network doesn't track)
```

- **MESA** simulates the actual physics of the star: gravity, pressure, convection, and a *small* nuclear reaction network good enough to get the hydrodynamics right. It produces a time history of every zone of the star (temperature, density, composition, radius, luminosity...).
- **ReacNetJl** takes a single zone's temperature/density history from MESA and re-runs it through a *much larger* nuclear reaction network, recovering isotopes MESA's own small network doesn't bother tracking (like ²²Na, ²⁶Al, ⁷Be) — these matter for gamma-ray astronomy even though they're irrelevant to the hydrodynamics.
- **NovaMessengers** (this repository's Julia package, `NovaMessengers/`) takes the output of both and answers the actual science question: what particles come out, do they escape or get absorbed on the way, and what does the resulting signal look like at any given moment?

---

## 3. The physics and equations

Every equation below is implemented in the code (file names given) with the same notation, so you can go read the actual implementation once the idea makes sense.

### 3.1 Radioactive decay — the core "messenger" mechanism

Most of what NovaMessengers tracks comes from unstable isotopes made during the explosion decaying away afterward. Each decay releases:
- a **neutrino** (always escapes, see below),
- usually a **positron**, which immediately annihilates with an electron and makes two 511 keV gamma-ray photons,
- sometimes a **gamma-ray line** at a specific energy, if the decay leaves the daughter nucleus in an excited state that then relaxes by emitting a photon.

The number of undecayed nuclei falls off exponentially:

```
N(t) = N0 * exp(-lambda * t)          lambda = ln(2) / half-life
```

so the rate of decays (= the rate of messenger emission) at any instant is just:

```
decay_rate(t) = lambda * N(t)
```

Two isotopes matter differently here: ⁷Be decays purely by *electron capture* — no positron at all, so it contributes zero to the 511 keV line, only to the neutrino signal and a weak 478 keV gamma line. ²²Na and ²⁶Al are a mix of positron decay and electron capture. The code tracks each isotope's `positron_branching` and `gamma_branching` explicitly rather than assuming "one positron per decay" for everything.

*Code:* `NuclearDecay.jl` (decay constants, branching ratios), `MessengerProduction.jl` and `ExtendedMessengers.jl` (turning decay rates into positron/neutrino/gamma production rates).

### 3.2 Neutrinos — the one undistorted messenger

A neutrino's mean free path through the star is enormously larger than the star itself, so:

```
P_escape(neutrino) = 1, always
```

Whatever rate of neutrinos gets produced *is* exactly the rate an observer would see (if anyone could actually build a detector sensitive enough — in reality nobody can, for a nova this far away, but it's the perfect "ground truth" signal to compare everything else against). No other messenger in this project has that property.

*Code:* `Transport.jl`.

### 3.3 Gamma-rays — escape is a fight against absorption

Unlike neutrinos, gamma-ray photons *do* interact with the gas around them — mostly through Compton scattering off electrons. The probability a photon born at some depth actually escapes depends on how much material is between it and the surface:

```
tau(zone) = integral of (density * opacity) from that zone out to the surface
P_escape(zone) = exp(-tau(zone))
```

The opacity itself depends on photon energy through the Klein-Nishina formula (the quantum-mechanically correct version of Compton scattering, which becomes energy-dependent above ~511 keV):

```
kappa(E) = kappa_Thomson * [sigma_KN(E) / sigma_Thomson]
kappa_Thomson = 0.2 * (1 + X)      (X = hydrogen mass fraction)
```

This is why a 511 keV line and a 1.275 MeV line, produced at the same depth, don't escape with the same probability, and why gamma-ray lines characteristically switch on late — only once the ejecta has expanded and thinned out enough for `tau` to drop below 1.

*Code:* `Transport.jl` (`klein_nishina_factor`, `compton_opacity`, `optical_depth_gamma`, `escape_probability_gamma`).

### 3.4 The white dwarf's own light — a simple glowing sphere

Both the white dwarf and its companion star emit light because they're hot, the same way a hot piece of metal glows. That's blackbody (Planck) radiation:

```
B_nu(T) = (2 h nu^3 / c^2) / (exp(h nu / kT) - 1)      (energy per area per solid angle per frequency)
L_E(E)  = 8 pi^2 R^2 E^3 / (h^3 c^2 [exp(E/kT) - 1])   (converted to energy per second per unit photon energy)
```

The remarkable thing is that MESA's own simulation already tracks the white dwarf's temperature and radius through the *entire* event — quiescent (~30,000 K), then swelling and cooling during the explosion, then shrinking back down to a small, extremely hot state afterward (this run reaches over 600,000 K) that's actually the physical origin of the observed "supersoft X-ray" phase real novae show. One formula, sampled at different times, covers four different phases of the nova's life.

*Code:* `QuiescentContinuum.jl` (`spectral_luminosity_ev`, `wd_photosphere_at`).

### 3.5 Shocks — a second, independent messenger channel

A nova doesn't eject its material in one smooth puff. A slower shell leaves first; a faster wind follows and catches up to it. Where they collide, a shock forms, and shocks can accelerate a small fraction of particles to very high energies — completely unrelated to the nuclear burning that made the isotopes above. This project uses a published model (Diesing & Metzger 2026 — see references) for that shock physics. A few of the key relationships:

```
wind mass-loss rate:   Mdot_w(t) = (M_env / tau) * exp(-t/tau)
shock velocity:        v_sh = v_wind / 2   (momentum-conserving collision)
shock temperature:     T_sh = (3 * m_proton / 16 k_B) * v_sh^2
cosmic-ray luminosity: L_CR = xi_CR * L_shock
gamma-ray luminosity:  L_gamma ~= f_Omega * xi_CR * kappa * L_shock   (calorimetric limit)
```

The accelerated protons collide with other protons in the shocked gas and produce pions, which decay into the GeV gamma-rays real nova telescopes (like Fermi-LAT) actually detect. The same hot, shocked gas also glows via ordinary thermal bremsstrahlung ("braking radiation" from electrons deflected by ions):

```
bremsstrahlung emissivity:  epsilon_ff(nu) = 6.8e-38 * Z^2 * n_e * n_i * T^(-1/2) * exp(-h*nu/kT) * gaunt_factor
```

*Code:* `ShockAcceleration.jl` — every equation in this module cites its equation number from the source paper directly in the code comments, so cross-referencing is straightforward.

### 3.6 Reaction energetics — bookkeeping the nuclear energy budget

Each nuclear reaction releases (or, rarely, absorbs) a specific amount of energy, computed from the difference in nuclear mass between what goes in and what comes out:

```
Q = (mass of reactants) - (mass of products)     [converted to energy via E=mc^2]
```

This is used mainly as a cross-check: sum up the energy released by every tracked reaction, zone by zone, and it should account for a specific (traceable) fraction of the total nuclear energy generation MESA reports on its own.

*Code:* `ReactionEnergetics.jl`.

---

## 4. Code map

```
Particle-Int-HydroMix-NOB/
├── mesa_work/wd_nova_burst_co/     # the MESA simulation itself (1.1 Msun CO white dwarf)
└── NovaMessengers/                 # the Julia package that does everything downstream
    ├── src/
    │   ├── MesaIO.jl                    # reads MESA's output files
    │   ├── NuclearDecay.jl              # half-lives, branching ratios for every tracked isotope
    │   ├── MessengerProduction.jl       # decay/positron/gamma rates for MESA's own 7 tracked isotopes
    │   ├── Transport.jl                 # neutrino (trivial) and photon (Compton) escape
    │   ├── ReactionEnergetics.jl        # per-reaction energy bookkeeping
    │   ├── SignalSynthesis.jl           # combines production + transport into light curves
    │   ├── ShockAcceleration.jl         # the Diesing & Metzger shock/cosmic-ray/gamma-ray model
    │   ├── TrajectoryPostProcessing.jl  # bridges a MESA zone to ReacNetJl and back
    │   ├── ExtendedMessengers.jl        # decay/gamma rates for ReacNetJl's extra isotopes (22Na, 26Al, 7Be)
    │   ├── QuiescentContinuum.jl        # the white dwarf's and companion's own blackbody light
    │   └── SpectralEvolution.jl         # composites every channel above into one spectrum per moment
    ├── examples/
    │   ├── trajectory_postprocessing.jl # runs MESA -> ReacNetJl
    │   └── spectrum_movie.jl            # runs everything -> an animated spectrum movie
    └── test/                            # unit tests
```

---

## 5. Running it

```bash
# One-time setup
cd /home/sgervais/Documents/Particle-Int-HydroMix-NOB
julia --project=NovaMessengers -e 'using Pkg; Pkg.develop(path=expanduser("~/Documents/ReacNetJl")); Pkg.instantiate()'

# 1. MESA hydro simulation (optional to rerun -- LOGS/ from a prior run is already valid)
source mesa_work/env.sh
cd mesa_work/wd_nova_burst_co && ./rn && cd ../..

# 2. ReacNetJl post-processing: recover isotopes beyond MESA's own small network
julia --project=NovaMessengers NovaMessengers/examples/trajectory_postprocessing.jl

# 3. NovaMessengers: build the full multi-messenger spectrum movie
julia --project=NovaMessengers NovaMessengers/examples/spectrum_movie.jl
# -> NovaMessengers/examples/plt_out/spectrum_movie.mp4
```

---

## 6. What's not built yet

Two phases of a real nova's life are deliberately not modeled yet, because doing them correctly needs more reading first rather than guessing at the physics:

- **Dust formation.** Some novae go dim in visible light for a while a few months in, because dust condenses in the cooling ejecta and blocks the view, while glowing brightly in infrared instead.
- **Radio emission.** The expanding, ionized gas eventually becomes transparent at radio wavelengths too, and shocks can also produce non-thermal radio emission the same way they produce gamma-rays.

Both are flagged as needing dedicated literature review before implementation — see the reference list below for the specific papers already identified as starting points.

There's also one physics limitation baked into the current pipeline worth understanding: `TrajectoryPostProcessing`/`ExtendedMessengers` follow only *one* representative zone of the star through ReacNetJl, not the whole envelope at once, so isotope abundances beyond MESA's own network describe "one representative parcel of gas," not a full spatial picture. This is standard practice in the nova nucleosynthesis literature, but it does mean some effects (e.g. the ²²Na/²⁶Al gamma-ray lines only becoming visible once ejecta far outside the tracked zone thin out) aren't fully captured by the current single-zone run.

---

## 7. References

**Simulation tools**
- MESA (Modules for Experiments in Stellar Astrophysics) — the stellar evolution/hydrodynamics code. [docs.mesastar.org](https://docs.mesastar.org)
- ReacNetJl — the nuclear reaction network post-processor (private repository, `~/Documents/ReacNetJl`).

**Physics papers used directly**
- Diesing & Metzger (2026), "A Unified Model for Shock Interaction and gamma-Ray Emission in Classical Novae" — the source for every equation in `ShockAcceleration.jl`.
- Chomiuk, Metzger & Shen (2021), "New Insights into Classical Novae," *Annual Review of Astronomy and Astrophysics* 59, 391. [arXiv:2011.08751](https://arxiv.org/abs/2011.08751) — the review whose Figure 1 phase timeline this project's spectrum movie is structured around.

**Identified for the not-yet-built dust and radio phases** (see `reference_dust_radio_nova_literature` project notes for the full annotated list):
- Derdzinski, Metzger & Lazzati (2017), MNRAS 469, 1314 — dust formation in the dense shell behind an internal shock.
- Hachisu, Kato & Matsumoto (2024), [arXiv:2402.08287](https://arxiv.org/abs/2402.08287) — a multiwavelength light-curve model reproducing simultaneous dust dips and supersoft X-rays.
- Metzger, Hascoët, Vurm, Beloborodov, Chomiuk, Sokoloski & Nelson (2014), MNRAS 442, 713, "Shocks in nova outflows — I. Thermal emission" — thermal X-ray/optical/radio emission from the same shock geometry as the Diesing & Metzger model already implemented.
- Weston et al. (2016), MNRAS 457, 887 — non-thermal (synchrotron) radio emission from the same colliding-flow shocks.
- Chomiuk, Linford, Aydi et al. (2021), *ApJS* — a survey of 36 novae's radio light curves over five decades.
- Seaquist & Palimaka (1977); Hjellming et al. (1979) — the foundational thermal free-free "expanding photosphere" model for nova radio emission.

**Nuclear/atomic data**
- Isotope half-lives, branching ratios, and neutrino energies: MESA's own `weak_info.list` for the 7 MESA-tracked isotopes; standard nuclear data compilations (ENSDF/NNDC-style) for ²²Na, ²⁶Al, ⁷Be beyond MESA's network (see `NuclearDecay.jl`'s own documentation for exact values and caveats).

---

## 8. Things left to implement or fix

**New physics, needs literature research first (do not guess)**
- [ ] **Dust/IR channel.** Candidate mechanism: Derdzinski, Metzger & Lazzati (2017) — dust nucleating in the dense shell behind an internal shock, same shock geometry as `ShockAcceleration.jl` already computes. Read at equation level before implementing.
- [ ] **Radio channel**, two components: thermal free-free (classic Seaquist/Hjellming expanding-photosphere model) and non-thermal synchrotron (Weston et al. 2016, same colliding-flow shocks as the GeV channel). Metzger et al. (2014) is the strongest candidate for extending `ShockAcceleration.jl` to predict both self-consistently alongside the GeV/bremsstrahlung channels already there.
- [ ] **Optical/UV line spectrum** (P Cygni profiles, forbidden emission lines). This is the single biggest realism gap relative to how real nova abundances are actually measured observationally — most published nova abundance results come from optical spectroscopy, not gamma-ray lines. Needs at least a basic line-formation/radiative-transfer approach, not just the blackbody continuum `QuiescentContinuum.jl` already provides.

**Known approximations worth revisiting**
- [ ] **Single-zone limitation.** `TrajectoryPostProcessing`/`ExtendedMessengers` track one representative mass coordinate through ReacNetJl, not the whole envelope — this is why the ²²Na/²⁶Al freeze-out gamma-ray lines currently integrate to ~0 in this run's saved window (the tracked zone never becomes locally transparent, even though the outer envelope does). Fixing this for real needs either a MESA run that resolves full homologous ejection, or an explicit model for how freeze-out isotopes get transported to the transparent outer layers.
- [ ] **Shock parameters are only partially calibrated from this run.** `calibrate_shock_params` derives the white dwarf mass, envelope mass, and ejecta velocity from MESA's own data, but the ejection timescale (`tau_days`) and the shock microphysics efficiencies (`xiCR`, `xiB`, `fX`, `kappa`, `fOmega`) are still the paper's fiducial literature values, since `wd_nova_burst_co` doesn't resolve the actual mass-ejection event.
- [ ] **Companion is a generic placeholder** (`PlanckSource(4000K, 3e10cm)` in `spectrum_movie.jl`) — swap in real parameters once a specific system is chosen to model.
- [ ] **²⁶Al isomer (`al*6`) not tracked.** Negligible in this run's own ReacNetJl output (~1e-16 vs ~1e-6 for the ground state) and its name doesn't fit the isotope-name parsing convention (`mass_number`) without a dedicated exception — low priority unless a future run makes it non-negligible.

**Housekeeping**
- [ ] `examples/reaction_energetics_check.jl` vs. `examples/reaction_energetics_crosscheck.jl`, and `examples/shock_model_validation.jl` vs. `examples/shock_spectrum_validation.jl` — similarly-named example files that were never confirmed to be intentionally distinct (model-level vs. spectrum-level checks look likely) rather than leftover duplication. Quick review, not urgent.
- [ ] `examples/spectrum_movie.jl`'s companion/shock/dust choices are all currently hardcoded constants at the top of the script — fine for one system, but worth turning into named presets (or CLI arguments) once more than one nova scenario is being modeled.
