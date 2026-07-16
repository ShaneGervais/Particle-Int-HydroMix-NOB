# Nova Outburst Messenger-Physics Pipeline

## Context

The goal is to understand what information particle "messengers" (neutrinos,
positrons/annihilation photons, radioactive-decay gamma-ray lines, and
eventually the optical/UV light curve) carry away from a classical nova
thermonuclear runaway (TNR) about the underlying nuclear burning, convective
mixing, and white dwarf (WD) properties that produced them — and how their
interactions during transport (Compton scattering, absorption, annihilation,
free-streaming) reshape that information before it's observable at Earth.

MESA (already installed at `~/mesa-25.12.1`, SDK at `~/mesasdk`, confirmed
working) will do the hydrodynamics, convective mixing, and nuclear-network
integration — reusing decades of validated stellar-evolution physics rather
than reimplementing it. MESA ships a ready-to-run `wd_nova_burst` test case
(1.1 Msun CO WD, accreting H/He at 1e-9 Msun/yr through one full TNR cycle,
~586s runtime) that already has a populated `LOGS/` output from a prior run
and a small existing edit in `inlist_pgstar`, so this is a live, in-progress
starting point, not a cold install.

Julia's job is everything downstream and novel: reading MESA's output,
computing what particles those nuclear reactions actually produce, modeling
whether those particles escape freely or interact first, and synthesizing the
resulting "observed" messenger signals. This is a fresh, standalone Julia
package — no other project directories are referenced or reused for
conventions.

Confirmed decisions: (1) NuPPN (NuGrid) will be used as a **trajectory
post-processor** — MESA keeps its small live network for the hydro solve,
and representative T/rho/time histories of burning zones get fed through
NuPPN's larger network afterward to recover isotopes MESA's live net doesn't
track (22Na, 26Al, 7Be, etc.), which is the standard approach in the nova
nucleosynthesis literature. (2) Phase 0 targets the CO WD baseline first,
since it's runnable today and lets the whole pipeline (Phase 1-4) get built
and validated against fast beta+ messengers before the ONe/22Na/26Al work.

## Directory layout

```
Particle-Int-HydroMix-NOB/
├── README.md                      # science motivation, how the two halves relate
├── .gitignore
├── mesa_work/                     # copied MESA work dirs (never edit the install itself)
│   └── wd_nova_burst_co/          # Phase 0: copy of the CO/cno_extras test case
└── NovaMessengers/                # new Julia package (name is a placeholder, easy to rename)
    ├── Project.toml
    ├── src/
    │   ├── NovaMessengers.jl
    │   ├── MesaIO.jl          # Phase 1
    │   ├── NuclearDecay.jl    # Phase 2a
    │   ├── MessengerProduction.jl  # Phase 2b
    │   ├── Transport.jl       # Phase 3
    │   └── SignalSynthesis.jl # Phase 4
    ├── test/
    │   ├── Project.toml
    │   ├── runtests.jl
    │   └── data/               # small truncated real MESA output samples
    └── examples/
```

git-init at the repo root. `.gitignore` for `mesa_work/*`: `LOGS/`, `photos/`,
`make/*.o`, `make/*.mod`, `make/*.smod`, `restart_photo`, the compiled `star`
binary, and regenerable `.mod` snapshots — but track all inlists,
`src/run_star_extras.f90`, `history_columns.list`, `profile_columns.list`,
the seed model file, and the `mk`/`rn`/`clean`/`re`/`ck` scripts.

## Phase 0: working baseline (MESA side)

1. `cp -r ~/mesa-25.12.1/star/test_suite/wd_nova_burst mesa_work/wd_nova_burst_co`
   — never edit the MESA install directly.
2. `./clean && ./mk` inside the copy to confirm it rebuilds standalone.
3. Two inlist edits before the first fresh run (both in `&controls`):
   - `history_interval = 1` (currently 5) — near burst peak, individual
     timesteps shrink to ~5-60s, which is already comparable to o14/f17
     half-lives, so interval-5 sampling is too coarse to resolve them.
   - Add `total_mass n13`, `total_mass o14`, `total_mass o15`,
     `total_mass f17`, `total_mass f18`, `total_mass ne18`,
     `total_mass ne19` to `history_columns.list` — this file already has
     `total_mass h1`/`total_mass he4` enabled and documents "as many as
     desired," so this needs zero Fortran changes and turns per-isotope
     abundance tracking (needed for Phase 2) from "needs custom hooks" into
     "already in history.data."
4. Lower `profile_interval` from 50 to ~5-10 for Phase 3's spatial/optical-
   depth needs (cheap: profiles are a few MB each).
5. **ne18 exception**: its half-life (1.7s) is shorter than a single MESA
   hydro timestep near burst peak, so it's in secular equilibrium within
   each saved step — finite-differencing its abundance won't recover a real
   production rate. No custom Fortran needed for this after all: MESA
   already exposes named per-reaction rates as profile columns
   (`screened_rate <name>` / `raw_rate <name>`, in reactions/second,
   computed internally as `s% screened_rate(i,k) * s% dm(k)` — already
   integrated over each zone's mass). Added `screened_rate r_f17_pg_ne18`
   and `raw_rate r_f17_pg_ne18` to `profile_columns.list` directly, which
   gives the ne18 production rate (nuclei/second per zone) with zero
   Fortran-side work.
6. Rerun (~10 min), confirm the new columns show a clean rise/decay through
   the burst window.

### Build notes (this host)

This test case is compiled with `pgstar_flag = .true.` by default, which
links against the SDK's `libpgplot.so` regardless of the runtime flag value
(pgstar controls whether it's *invoked*, not whether it's *linked*). Two
environment-specific fixes were needed to get a working build here, both
already applied in this copy:

- `mesa_work/env.sh` exports `MESA_DIR`/`MESASDK_ROOT` and sources
  `mesasdk_init.sh` — `.bashrc` only sets `MESA_DIR`, not the SDK env.
  Source it before any MESA build/run command:
  `source mesa_work/env.sh`.
- `mesa_work/wd_nova_burst_co/make/makefile` sets
  `LOAD_EXTRAS = -L/usr/lib/x86_64-linux-gnu -lX11 -lxcb -lXau -lXdmcp -lrt -ldl -lpthread`.
  The SDK's bundled `ld` doesn't auto-resolve pgplot's transitive X11/xcb
  dependencies on this host even though the libraries are present
  system-wide, so they need to be listed explicitly.
- `pgstar_flag` is set to `.false.` in `inlist_wd_nova_burst` for headless
  batch runs (no `DISPLAY` here). This doesn't affect `LOGS/` output at all
  — live plotting is separate from data output, and `plot.py` (via
  `mesa_reader`) or `NovaMessengers` can still make static plots from
  `history.data`/`profile*.data` afterward.

## Phase 1: Julia MESA-output reader

MESA's `history.data`/`profile*.data` format: line 1 = numeric column-index
row (ignorable), line 2 = header field names, line 3 = header field values
(mixed quoted strings/floats), line 4 = blank, line 5 = numeric column-index
row (ignorable), line 6 = main-table column names, line 7+ = whitespace-
delimited data rows. Files are MB-scale (hundreds of rows x ~70-100
columns) — no need for a streaming parser.

- Hand-write the small 6-line header parser, then use `CSV.jl`/`DataFrames.jl`
  (`delim=' ', ignorerepeated=true`) for the bulk row parsing.
- `read_history(path) -> (header::NamedTuple, data::DataFrame)`
- `read_profile(path) -> (header::NamedTuple, data::DataFrame)` — profile
  zones are numbered from 1 at the **surface**, not the center.
- `read_profiles_index(dir)` — parses `LOGS/profiles.index` to map profile
  files to model number/age via a join against `history.data`.
- A `MesaRun` struct bundling a work-dir path + lazily-loaded history/
  profiles, designed to support multiple runs from day one (needed later for
  Phase 4 parameter sweeps).
- Column indices always resolved from the file's own header row into a
  name→index map — never hardcoded by position, since enabled columns (and
  therefore column count) change with `history_columns.list`/
  `profile_columns.list` and the active net.

## Phase 2: particle production

MESA ships `~/mesa-25.12.1/data/rates_data/weak_info.list`, which has
authoritative half-life and average-neutrino-energy (`Qneu`, MeV) values for
exactly the 7 candidate messenger isotopes (n13, o14, o15, f17, f18, ne18,
ne19) — the same table MESA itself uses internally. Freeze a copy into
`data/weak_info_subset.csv` in the repo so the package doesn't require a
MESA install at analysis time, and so the Julia physics stays numerically
consistent with what this exact MESA version used (avoids two independent
half-life tables silently disagreeing).

- **`NuclearDecay`**: `DecayIsotope(name, daughter, halflife_s, decay_const,
  q_neu_mev, mode)` built from the frozen table. No MESA-run dependency —
  independently unit-testable against literature half-lives.
- **`MessengerProduction`** (depends on `MesaIO` + `NuclearDecay`):
  - Neutrino/positron production rate per isotope from finite-differencing
    `total_mass <iso>` (now present per Phase 0) — directly viable for n13
    (598s) and f18 (7109s), noisier but usable for o14/o15/f17 (comparable
    to timestep spacing), and **not used for ne18** (uses the
    `screened_rate r_f17_pg_ne18` profile column from Phase 0 step 5 instead,
    read per-zone from `profile*.data` rather than finite-differenced from
    `history.data`).
  - `annihilation_photon_rate` = 2x positron rate (site-of-annihilation
    question deferred to Phase 3).
- **NuPPN trajectory extension**: for isotopes beyond the live net (22Na,
  26Al, 7Be), extract T(t)/rho(t) histories for the burning zones from
  Phase-1 profile data and hand them to NuPPN as fixed trajectories, rather
  than expanding MESA's live network. This works independent of the Phase 0
  CO-WD choice, but meaningful 22Na/26Al yields likely still require an ONe
  WD trajectory (Ne-rich core dredge-up supplies the seed nuclei) — flagged
  as a later phase once the CO-WD pipeline is validated end-to-end. Exact
  NuPPN invocation (trajectory file format, network choice) to be worked out
  when this phase starts, since it depends on the installed NuPPN version.

## Phase 3: transport/interaction

- **Neutrinos**: escape probability = 1, always — the Phase 2 production
  rate *is* the observable signal, unmodified. This is a physics finding to
  state explicitly (neutrinos are a direct real-time probe of the burning
  zone), not just an implementation shortcut.
- **Photons**: escape probability from local optical depth to a discrete
  MeV-scale line photon, whose dominant interaction is Compton scattering
  (roughly energy-independent above ~200 keV) — not the Rosseland-mean
  opacity MESA's own `tau`/`log_opacity` columns represent (that's for the
  thermal radiation field). Compute a separate
  `tau_gamma(zone) = integral(rho * kappa_Compton, dr)` from MESA's
  `logRho`/`radius`/`mass` structure with a fixed Compton opacity
  (~0.2-0.4 cm^2/g), then `P_escape(zone) = exp(-tau_gamma(zone))`.
- **Positrons**: stopping range is orders of magnitude shorter than a photon
  mean free path, so default to "annihilate in-situ, then transport the
  511 keV photons through the same Compton-escape machinery" — true
  in-flight annihilation flagged as a later refinement.
- Module layout: `Transport.jl` with clearly separate neutrino
  (trivial pass-through) vs. photon/positron (Compton escape) code paths, so
  the "neutrinos are physically different" decision stays visible in the
  code structure rather than being hidden inside one generic function.

## Phase 4: signal synthesis

- `neutrino_lightcurve(run)`: sum of Phase 2 production over all tracked
  isotopes vs. `star_age`.
- `gamma_lightcurve(run; line_energy)`: Phase 2 zone-resolved production x
  Phase 3 escape probability, summed over zones vs. time — this is the plot
  that directly answers "does the signal preserve or reshape the underlying
  nuclear information."
- Coarse time-resolved line spectra (511 keV now; 1.275/1.809 MeV once the
  NuPPN extension lands) — full Doppler/line-broadening from ejecta velocity
  explicitly deferred.
- Parameter sweeps (mixing efficiency, WD mass, CO vs ONe) require
  additional MESA runs, not just more Julia analysis — `MesaRun`/a
  `MesaRunSet` batch loader (Phase 1) is designed from the start so Phase 4
  just adds sibling `mesa_work/` directories and re-runs Phase 1-4 over each.

## Julia package tooling

- One package: `Project.toml` + `src/` + `test/` + `examples/`.
- Core deps: `DataFrames.jl`, `CSV.jl` (MESA output parsing), `Statistics.jl`,
  a plotting package (`CairoMakie.jl` recommended — headless-friendly).
  `Unitful.jl` is worth considering for the `NuclearDecay`/
  `MessengerProduction` constants layer, since the pipeline spans
  years-to-sub-second timescales and MeV-to-erg conversions.
- `DifferentialEquations.jl`/SciML is **not needed**: MESA already performs
  the hydro+network ODE integration, and Julia's role here is post-
  processing MESA's output, not re-solving the nucleosynthesis ODEs.
- Testing:
  - `NuclearDecay` unit tests: parsed half-lives vs. literature/NNDC values.
  - `MesaIO` integration tests against a small truncated real sample (can be
    built today from the existing `wd_nova_burst/LOGS/history.data` and
    `profile1.data`, no fresh run required).
  - `MessengerProduction`/`Transport` unit tests: analytic checks (pure
    exponential decay matches `N0*exp(-t/tau)`; escape probability -> 1 as
    tau -> 0, -> 0 as tau -> infinity).

## Verification

1. **Energetics check**: integrate Phase 2 neutrino energy loss
   (sum of decay_rate x Qneu) over the burst, compare to total nuclear
   energy release (`10^log_Lnuc` integrated, already in history.data) —
   should be a small but non-negligible fraction, consistent with hot-CNO
   literature.
2. **Internal MESA self-consistency check**: `eps_nuc` already has
   reaction-neutrino losses folded in; the gap between a naive Q-value
   energy-generation estimate (from the already-enabled `pp`/`cno`/
   `tri_alpha` profile columns) and `eps_nuc` should match the independently
   reconstructed neutrino-loss rate from Phase 2 — checkable from this exact
   run with no external literature number needed, and should gate trust in
   Phase 3/4 before going further.
3. **Reproducibility**: since a valid `LOGS/` already exists from a prior
   run, do the Phase 0/1 read-and-parse check against it first, and only
   spend the ~10 minute rerun once the `history_interval`/`total_mass`
   inlist edits are in place.
4. **Gamma-line check (NuPPN phase)**: once 22Na/26Al are available, compare
   predicted fluxes at a literature-typical distance against published nova
   gamma-line predictions — order-of-magnitude agreement is the bar, given
   Phase 3's first-order `exp(-tau)` transport.

## Critical files

- `~/mesa-25.12.1/data/rates_data/weak_info.list` — authoritative half-life/
  Qneu source for `NuclearDecay`
- `~/mesa-25.12.1/star/test_suite/wd_nova_burst/history_columns.list` and
  `profile_columns.list` — files to copy and edit for Phase 0
- `mesa_work/wd_nova_burst_co/profile_columns.list` — where the ne18
  `screened_rate r_f17_pg_ne18` column is enabled
- `mesa_work/wd_nova_burst_co/make/makefile` — `LOAD_EXTRAS` link workaround
  for this host's SDK/X11 linking (see below)
- `NovaMessengers/src/MesaIO.jl` — Phase 1 reader everything downstream
  depends on
- `NovaMessengers/src/NuclearDecay.jl` — Phase 2 decay-data anchor
