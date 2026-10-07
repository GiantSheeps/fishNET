# fishNET

Hybrid Lagrangian-Eulerian marine ecosystem model. Fish are tagged super-individuals in IBM
cells and stage-structured biomass (with four moments of every heritable trait) elsewhere;
both representations feed, respire, die, mature and spawn through the same equations and
couple two-way to an NPZD ocean driven by prescribed physics. Nitrogen is conserved to round-off.

    python fishnet.py namelist.toml              # run (output in runs/<name>/)
    python fishnet.py namelist.toml --resume     # continue it from its newest restart file (see Restarts)
    python fishnet.py namelist.toml --ensemble        # parameter sweep (output in runs/<name>_ensemble/)
    python fishnet.py namelist.toml --ensemble 40     # ...with 40 members, overriding the namelist
    python plot.py runs/<name>                   # figures + animations (add --no-anim to skip movies)
    python plot.py runs/<name>_ensemble          # ensemble spread, sensitivity heatmap, scatter
    python validate.py                   # 60 checks, ~60 s (add --long for a 90-day consistency test)
    python dashboard.py                  # browser GUI: launch runs and watch them as they go

## Files
- `fishnet.py` model: grid, ocean forcing (synthetic basin or netCDF), NPZD, agents, biomass fields,
  behaviour, genetics, hybrid regulator, output
- `plot.py` time series, stages, trait evolution, losses by cause, maps, trait maps, Hovmoller,
  life histories, agent tracks, depth use, and two MP4 animations
- `dashboard.py` + `webui/` interactive browser dashboard: start a run, animate it while it runs
- `get_forcing.py` (GLORYS, needs an account), `get_forcing_hycom.py` (HYCOM, no account) real ocean forcing
- `validate.py` conservation, positivity, transport, genetics, aggregation round trips,
  determinism, netCDF forcing, IBM-vs-biomass consistency
- `namelist.toml` run conditions; `species/*.toml` species parameters and trait distributions

## Output (netCDF)
- `fields.nc` temperature, currents, NPZD, fish biomass/numbers by species and stage, agent biomass,
  local mean breeding values, IBM mask
- `series.nc` domain totals every step: biomass, numbers, agents, trait means/SDs, losses by cause,
  catch, eggs, plankton inventories, N budget error
- `agents.nc` snapshots of every agent: id, parents, generation, stage, sex, position, depth,
  n represented, structure/reserve/gonad mass, length, age, breeding values and phenotypes
- `lifehist.nc` one record per agent at death/merge/end: cause, birth and end place and time, traits

## Spin-up
Before day 0 the model steps for `run.spinup_days` (default 30; 0 turns it off) while the forcing replays
its first `run.spinup_cycle_days` (default 7). Nothing is written during the spin-up. The clock runs from
-spinup_days to 0, so there is no fishing and fish born then have negative birth days. At day 0 the nitrogen
budget, external fluxes and life-history records start from the settled state. This takes away the first
weeks' plankton bloom and the fish die-off from the initial conditions. It does not bring fish populations
to equilibrium, which takes years. A resumed run does not spin up again.

## Initial state
- `run.initial_biomass = "literature"` replaces each species' namelist `biomass` with the global estimate in
  `observations/literature_biomass.toml` (spread as tonnes x 1e6 / 3.6e14 m2, the namelists' convention).
  Assessed species take the sum of their RAM Legacy stocks (1993-2021 mean), a lower bound because
  unassessed stocks are missing. Lanternfish, bristlemouth, polar cod and silverfish take published
  estimates. Each entry cites its source. Starting from RAM and then validating against RAM is partly
  circular, so judge the end state, not the start. `"namelist"` (default) keeps the old behaviour.
- `run.warm_start = "runs/<name>"` (or a restart `.pkl`) starts from an earlier run's fish instead: biomass
  classes, agents with their genes, ages and generations, and (with `warm_start_plankton = true`, the
  default) its plankton and nutrients, and its oxygen unless `npzd.O2_init = "woa"`. The clock goes back to day 0 and birth days shift with it.
  Species are matched by name, so you can retune, add or drop species; new ones start fresh. The grid must
  match. Use it to carry a years-long spin-up into many shorter experiments.
- `npzd.O2_init = "woa"` starts oxygen from the World Ocean Atlas 2023 annual mean
  (`npzd.O2_file`, downloaded from NCAR's RDA mirror, dataset d285000), averaged onto the grid, and makes
  the sea bed relax toward it instead of the constant `O2_deep`. This gives real low-oxygen zones.

## Fishing
`[fishing] mode = "ram"` (used by `global_species.toml` since 2026-09-26) fishes each species at its RAM
Legacy fishing mortality **for the calendar year** of the model date, **only in the FAO major areas where it
has assessed stocks**, times `scale` (1.0 = the assessed rate). `year = 1993` fixes the year (for spin-ups);
`unassessed` applies that fraction of the species' mean F elsewhere (default 0). The table,
`observations/ram_fishing_series.json`, holds F per species, FAO area and year (1993-2021), from
-ln(1 - catch/biomass) per stock averaged over the stocks assessed in each area, weighted by their biomass;
years outside it take the nearest year. `mode = "constant"` (the default) uses each species' `fishing_F`
everywhere, all the time; in `global_species.toml` those values are RAM F x 0.5, kept for comparison with
the 2026-09-24 runs. Regenerate all the tables with `python obsval.py --tables`.

## Hindcast (real ocean years)
`get_forcing_global.py --start 1993-01 --end 2012-12 --out forcing_glorys_1993_2012` streams GLORYS monthly
means (about a minute a month), and `make_forcing_sequence.py --in forcing_glorys_1993_2012/monthly --out
forcing_glorys_1993_2012/forcing.nc` joins them into a dated file (mid-month stamps). With `run.start` at the
file's first month the model follows the real calendar: interannual variability, El Nino, warming. The
prepared pair of runs is `runs/_configs/spinup_1993.toml` (5 years on the looped climatology at 1993 fishing,
from literature biomass) and `runs/_configs/hindcast_1993_2012.toml` (warm-started from it, year-by-year
fishing). The Validation tab's year-by-year chart compares the model with each species' assessed stocks.

## Restarts
Every `run.restart_every_days` (default 30; 0 turns it off), and again at the end of the run, the full
model state is written to `runs/<name>/restart/restart_dayDDDDD.D.pkl`. That covers fish and plankton,
agents with their genes and learned behaviour, the random stream and any shallow-water ocean state.
`run.restart_keep` newest files are kept (default 2; 0 keeps all). A file starts at about 25 MB on the
48x32 Atlantic grid and 250 MB on the 180x80 global grid. It grows during the run, because it also holds
the life-history records that `lifehist.nc` gets written from at the end.

    python fishnet.py namelist.toml --resume                                   # newest restart file
    python fishnet.py namelist.toml --resume runs/<name>/restart/restart_day00365.0.pkl

The run continues in the same folder and appends to its netCDF files. Anything written after the restart
point (for example by a run that crashed a few days later) is dropped first, including chunks the crash
may have corrupted, so the files read as one unbroken run. A stopped-and-resumed run writes exactly what
an uninterrupted one would (`validate.py` checks this). `run.days` in the namelist is the end day, so
raising it and resuming extends a finished run. Output and restart settings and the ocean come from the
namelist; species parameters come from the restart. `dt_hours`, `start` and `[grid]` must not change.

## Dashboard
`dashboard.py` serves a browser GUI on http://127.0.0.1:8765 (standard library only; nothing to install
beyond what the model already needs):

    python dashboard.py                  # opens a browser
    python dashboard.py --port 9000 --no-browser --host 127.0.0.1

Runs are launched with the same interpreter the dashboard itself runs under, so activate your
environment first if the model needs packages from it (`global.toml` needs `global-land-mask`).

Left panel launches runs: pick a namelist, set name, days, dt, seed, output intervals, grid size, warming
rate and each species' initial biomass and fishing mortality, or edit the namelist text directly. The
copy that actually runs is archived in `runs/_configs/<name>.toml`. Model stdout streams into the log
pane, with a progress bar, ETA and a per-species table of the latest step. **Resume** continues the run
selected at the top from its newest restart file (it asks for the day to run until). **Start** with an
existing run name wipes that run folder, restarts included.

The main panel reads the netCDF output **while the model is still writing it**, so everything fills in as
the run goes (HDF5 file locking is disabled in the server for this; the files are re-opened as they grow):
- *Map & animation*: any gridded field — fish biomass or numbers by species and stage, agent biomass,
  local trait means, temperature, NPZD, oxygen, PAR, current speed, bathymetry — animated with play/scrub
  /fps controls, agent super-individuals as dots coloured by species and sized by fish represented,
  current arrows, IBM cells, and a colour scale that can be fixed across the run so frames are comparable.
  "Follow live" pins playback to the newest frame; uncheck it to replay history while the run continues.
  Space plays/pauses, arrow keys step frames.
- *Populations*: biomass by species and stage, egg production, catch, mortality by cause.
- *Ecosystem*: NPZD inventories, oxygen minimum, carbon export, metabolic index, IBM cells, N budget error,
  behaviour budget.
- *Evolution*: trait means with ±1 SD bands, generations, F_ST.
- *Figures*: runs `plot.py` for the selected run and shows the PNGs (and MP4s if you ask for them).

One run at a time; "Stop" terminates it and the partial output stays readable. Runs already in `runs/`
can be browsed and animated the same way without starting anything.

## Skill assessment
`skill.py` scores a run against observations and writes a scorecard:

    python skill.py runs/<name> observations.toml
    python skill.py runs/<name> observations.toml --climatology runs/control --ensemble runs/x_ensemble/ensemble.nc

`observations.toml` holds `[[target]]` entries (single numbers, same format calibration uses), `[[series]]`
entries (observed time series), an optional `[size_spectrum]` slope, and `[[event]]` entries (an observed
heatwave or collapse with a date and tolerance). The report gives bias, RMSE, centred RMS, correlation,
sd ratio and Nash-Sutcliffe per series, a Taylor diagram, the community size-spectrum slope, detected
heatwaves and collapses matched against the observed ones, and CRPS with rank histograms if you pass an
ensemble. Note that anomalies need a climatology: on runs shorter than two years, pass `--climatology`
with a control run, or the fitted seasonal cycle will swallow the event you are looking for.

## Calibration
`calibrate.py` fits species parameters to observations:

    python calibrate.py namelist.toml targets.toml --method abc --members 60 --processes 4
    python calibrate.py namelist.toml targets.toml --method emulator --members 60 --chain 20000

Observations go in `targets.toml`, one `[[target]]` each (biomass, catch, chlorophyll, export or metabolic
index) with a species, day range, value and uncertainty; priors go in the namelist under
`[calibration.prior]`, using the same parameter paths as an ensemble sweep. Methods: `abc` (rejection on a
Latin-hypercube prior), `smc` (rounds resampling and perturbing survivors) and `emulator` (Gaussian process
fitted to the loss surface, then MCMC on it, with a leave-one-out R2 so you can tell whether to trust it).
Everything lands in `posterior.nc` with a printed report of observed against predicted.

Run members at the same length as the observations describe: a mismatch there biases the posterior.

## Where species live
Nothing prescribes a species' range. Distribution follows from the same machinery that drives everything
else: thermal preference (`T_opt`, `T_width`), food, predation risk, schooling, an evolvable preference for
shallow water (`w_shelf`, positive = shelf-associated, negative = open ocean), natal homing (ripe fish
return to their own birth place) and mate finding (a ripe female needs an adult male in her column, so
strays beyond the core range fail to reproduce). Spawning grounds are not listed in the species files;
biomass classes home up the gradient of ripe conspecifics, so aggregations form and persist on their own.

## Global configuration
`global.toml` runs the whole ocean: real coastlines (needs `pip install global-land-mask`), depth built from
distance offshore, seasons out of phase between hemispheres, and a wind-driven circulation with subtropical
and subpolar gyres in every basin plus a circumpolar current (`tau0`, `acc_u`, `wbc_cells`). At 2 degrees
with 8 levels that is 9,613 ocean columns at roughly 2.5-3.5 s per 12-hour step on one core.

    python fishnet.py global.toml

Stocks in `global.toml` are initialised from observed global biomass and use mortalities fitted by ABC, so
totals stay within 0.6-0.85x of observed over 90 days. The species files themselves are still tuned for a
North Atlantic basin; recalibrate with `calibrate.py` before trusting regional detail.

## Bathymetry
Columns have a bottom depth: cells hold only the water above the bed (partial bottom cells), so volumes,
light, mixing, sinking and deep nutrient supply all follow the topography. The synthetic ocean grows a
continental shelf and slope (`ocean.shelf_cells`, `ocean.shelf_depth`); the netCDF path reads a depth
variable (`ocean.bathymetry`) or infers the bed from the deepest valid level. Species with
`demersal = true` keep their juveniles and adults in the bottom cell. Horizontal transports are built from
a per-level streamfunction masked by the topography, so they are non-divergent by construction and a
uniform tracer stays uniform over a shelf.

See `ROADMAP.md` for the development plan toward a research model and what is done so far.

## Ocean forcing (`ocean.source`)
- `"shallow_water"`: an eddying layered ocean integrated with P. Connolly's Lax-Wendroff scheme
  on a grid `sw_refine` times finer than the ecosystem grid, forced by a double-gyre wind and diabatic
  heating (meridional gradient, seasonal and diurnal) that sets the upper-layer thickness. Spun up for
  `sw_spinup_days` once and cached in `sw_cache/`. What the ecosystem sees:
  - currents: the non-divergent part of each layer's flow (least-squares streamfunction fit, zero flow
    through coasts), mapped onto depth levels by layer thickness, so tracers and fish are stirred by
    gyres, meanders and eddies without spurious convergence;
  - temperature: SST anomalies from thermocline depth (warm deep rings, cold shallow ones), a seasonal
    mixed layer, the thermocline at the layer interface, and a daytime warm skin;
  - nutrients: upwelling where the thermocline shoals (subpolar gyre, cyclonic features), plus the coast.
  `sw_layers = 1` (1.5-layer reduced gravity) is the stable default. `sw_layers = 2` (stacked 2.5-layer)
  is experimental: the explicit cross-layer pressure coupling became unstable near the north-west corner
  in testing, even with a predictor-corrector.
- `"synthetic"` (default, the one species parameters are calibrated for): analytic double gyre with a
  seasonal thermocline and coastal upwelling (fast, no eddies).
- `"netcdf"`: reanalysis or your own fields. Variable names are auto-detected for CMEMS/GLORYS
  (`thetao`, `uo`, `vo`), NEMO and ROMS, or set `ocean.names`. Longitude convention, latitude order,
  packed integers, epoch and subsetting to the grid are handled; `ocean.bathymetry` reads a static
  `deptho` file. Two helpers:

      python get_forcing.py --preset north_atlantic --start 2010-01-01 --end 2019-12-31   # real GLORYS
      python get_forcing_hycom.py --preset north_atlantic --start 2010-01-01 --end 2011-12-31  # real HYCOM
      python make_test_forcing.py --out test_forcing                                      # offline lookalike

  `get_forcing.py` needs a free Copernicus Marine account and `pip install copernicusmarine`; it prints a
  namelist fragment matching what it downloaded. `get_forcing_hycom.py` needs **no account**: it pulls the
  HYCOM + NCODA GOFS 3.1 reanalysis (1994-2015, 1/12 degree, 3-hourly) over OPeNDAP, thinning in space
  (`--stride`), time (`--every`) and depth (`--levels`) as it goes, and writes the file incrementally so an
  interrupted download is still usable. HYCOM is the quicker route to a realistic run; GLORYS is the better
  product for published work. `make_test_forcing.py` writes invented data in real CMEMS packaging, for
  testing the reader without an account. Note the reader has not yet been run against real
  GLORYS files. Currents are projected onto their non-divergent part, which keeps about 77% of the mean
  speed; narrow boundary currents are smoothed unless the grid resolves them.

## Ensembles
`[ensemble.sweep]` accepts any namelist key (`npzd.mu_max`), species parameter (`species.cod.K_nursery`,
or `species.*.m0` for all species), trait component (`species.sardine.traits.T_opt.2` = heritability)
or list item (`ocean.heatwave.0.amplitude`). Values are lists, ranges `{min, max, n, log}` or factors on
each target's own value `{scale = [0.5, 1, 2]}`. Methods: `grid` (all combinations), `oat` (one at a
time), `lhs` (Latin hypercube, `members`), `random`; `replicates` repeats members with new seeds and
`processes` runs them in parallel. Each member writes a normal run folder (agent snapshots off unless
`keep_agents`), and `ensemble.nc` stacks every member's time series with its parameter values.

## Biogeochemistry
The NPZD carries a fifth tracer, dissolved oxygen: produced by photosynthesis, consumed by zooplankton,
remineralisation and fish respiration, exchanged with the atmosphere (Wanninkhof transfer velocity,
Garcia-Gordon solubility, `ocean.wind_speed`) and resupplied at the sea bed (`npzd.O2_deep`, or the World
Ocean Atlas field with `npzd.O2_init = "woa"`).
Remineralisation itself slows as oxygen runs out (`npzd.kO2_remin`).

Each species has a metabolic index, phi(O2, T, W): oxygen supply over resting demand. Below `phi_crit` the
aerobic scope for feeding shrinks, and below phi = 1 fish suffocate (`m_hypoxia`), a new cause of death.
Warming, growth and deoxygenation all push phi down, so the model reproduces metabolic habitat squeeze.
Species-level knobs: `A_o`, `phi_crit`, `m_hypoxia` (and optionally `E_o`, `eps_o`).

Carbon export past `npzd.export_depth` is diagnosed each step with Redfield C:N and split into sinking
detritus, fish faeces and carcasses released below that depth, and fish respiration there (the active
transport a migrating fish performs). These go to `series.nc` and the `biogeochem.png` figure.

## Key model choices
- Bioenergetics: C_max = cmax W^(2/3), respiration ~ W^0.8; respiration is set so fish at
  feeding level f_crit are at maintenance at L_inf. Surplus fills reserve, then structure or gonads.
- Feeding: multi-prey Holling II on plankton and fish classes, weighted by diet and a log-normal
  prey/predator length window; fish prey are gape-limited (ratio <= 0.8). `fish_diet` keys name species
  loosely (`warm_sardine`, `"warm sardine"` and `Warm-Sardine` are the same species), as do the
  `species.<name>.*` paths used by ensembles and calibration; a key matching no species in the run is
  listed in the startup banner rather than silently dropped.
- Feeding on plankton stops with size: above `plankton_L_max` (cm, smooth cut-off) a species eats only fish
  and the benthos proxy, so big predators depend on forage fish instead of grazing copepods. With
  `npzd.detritus_food = "bed_flux"` the benthos proxy (detritus) is available only on the sea bed and only as
  much as sinks onto it each step.
- Juveniles in biomass cells carry `hybrid.juvenile_bins` progress bins toward maturity. With one bin a fixed
  fraction matures each step whatever its age, so some fish mature almost at once (and generation counts
  inflate); K bins give an Erlang-K time to maturity with the same mean.
- Mortality: size-dependent natural, larval, thermal stress, fishing (logistic selectivity),
  predation, starvation, senescence; larval and juvenile mortality rises with same-stage density
  (`K_nursery`, g m-2), and that crowding penalty fades with length as exp(-L / `crowd_length`)
  (default a quarter of the maturity length), so it bites on larvae and small juveniles, not on large ones.
  Adults have their own crowding term, `K_adult` (g m-2 of same-species adults at which their natural
  mortality doubles; set at ~5-10x typical adult densities so it only stops runaway populations).
- Spawning grounds: `spawning_grounds = [[lon, lat], ...]` in a species file pulls ripe fish toward the
  nearest one (they may still spawn anywhere); `w_home` sets how strongly. This anchors stocks without walls.
- Starvation: agents die when their own reserve is exhausted. Biomass cells carry a reserve pool alongside
  structure and gonad, use the same allocation rule, and assume individual reserves spread uniformly within
  `hybrid.reserve_spread` (50%) of the class mean; exactly the fraction whose reserve would go negative
  starves, so the two representations agree (IBM-vs-biomass validation within ~5%) and mass still closes.
- Genetics: `genetics.loci` unlinked diploid loci per fish, with a pleiotropy matrix built from
  `[genetics.correlations]`, free recombination at meiosis and per-locus mutation. Genetic correlations
  therefore emerge from shared loci, and the realised G-matrix, F_ST between latitude bands and selection
  differentials by cause of death are reported each step. Set `loci = 0` for the older infinitesimal model (midparent breeding value + N(0, h2/2) segregation), phenotype =
  breeding value + environmental noise. Biomass cells carry offspring moments analytically.
- Mating and population structure: a ripe female looks for a mate in her own grid cell
  (`genetics.mate_search = "column"`, the default), so gene flow is limited by dispersal and regional races
  can form; `"global"` restores the old behaviour, where a female with no local male mated with any
  conspecific in the domain and divergence was continually erased. With `mate_choice_sd` she samples
  `mate_choice_sample` males and prefers those like herself in `mate_choice_traits` (assortative mating,
  width in standardised breeding-value units; 0 = random mating). `spawn_shift` is a heritable trait in
  days that shifts an individual's spawning window inside the species season, so populations can also
  diverge in *when* they breed, not only where. Note that within one biomass cell a trait distribution is
  carried as four moments, so two incipient forms sharing a cell are smeared into one: divergence that
  needs sympatry (rather than space or timing) requires `hybrid.mode = "all"`, which keeps individual
  genotypes everywhere.
- Behaviour (`behavior.engine`):
  - `"probabilistic"` (default, ported from the original fishNET): each decision, a fish computes an urgency
    for every behaviour from its stimuli (hunger from its reserve; local predation risk; thermal discomfort;
    isolation for schooling species; light for light-shy species; pre-spawning ripeness) and picks one with
    P ~ |drive weight| x urgency (plus a small random-move weight). It then heads for that behaviour's target:
    most food, least risk, best temperature, most company, darkest depth, or its spawning ground. Decisions
    persist for 1..`inertia_steps` steps. Drive weights are the heritable `w_*` traits (sign = toward/away),
    so drives can evolve. Biomass cells use the exact expectation of the same choice.
  - `"softmax"`: scores every candidate move by the weighted cues; `behavior.mode = "learn"` adds lifetime
    REINFORCE on top of the innate policy.
- Hybrid: IBM cells from static regions and/or SST-chlorophyll fronts; agents leaving IBM cells
  aggregate; biomass entering IBM cells is disaggregated at regulator updates, with trait mean and
  variance reproduced exactly and skew/kurtosis via Cornish-Fisher.

## Notes
- Species parameters are calibrated with the probabilistic engine for the synthetic North-Atlantic-like
  basin: all four species persist in bounded seasonal cycles over two years (0.3-3.7x initial biomass).
  The softmax engine runs with the same parameters but was calibrated earlier; re-tune for other forcing.
- Species parameters were calibrated on the synthetic ocean; the shallow-water ocean is more productive
  (thermocline upwelling), so expect different equilibria and re-tune if needed.
- `agents.nc` grows quickly (~200 bytes per agent per snapshot); set `agent_output_every_hours`.
- netCDF forcing needs `temp`, `u`, `v` on (time, depth, lat, lon); rename via `ocean.names`.
