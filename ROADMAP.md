# fishNET roadmap to a research model

Status legend: **[done]** implemented and validated | **[wip]** partly implemented | **[todo]** not started.
Updated every time the code is delivered. Validation counts refer to `python validate.py`.

---

## 1. Real forcing (reanalysis, bathymetry, shelves)

| Step | Status | Notes |
|---|---|---|
| 1.1 netCDF forcing path (temp, u, v) | [done] | `ocean.source = "netcdf"`, arbitrary variable names, land/below-bottom fill, time interpolation and cycling |
| 1.2 Bathymetry and shelves | [done] | partial bottom cells (`g.depth`, `g.wet`, `g.vol`, `g.bed`), shelf/slope in the synthetic ocean (`shelf_cells`, `shelf_depth`), depth read or inferred in the netCDF path, `demersal = true` pins species to the bed; transports rebuilt from a topography-masked per-level streamfunction so they stay exactly non-divergent |
| 1.3 CMEMS/GLORYS reader | [done] | `get_forcing.py` downloads GLORYS subsets via the copernicusmarine toolbox; the reader auto-detects CMEMS/NEMO/ROMS variable names, fixes longitude convention and latitude order, subsets to the grid, decodes packed integers and the 1950 epoch, and reads `deptho` from a static file. `make_test_forcing.py` writes a GLORYS-lookalike file for testing without an account. **Not yet tested against real GLORYS data** (no network access here) | depth/lat/lon/time renaming, `mlotst` mixed layer, `so` salinity, monthly-to-daily interpolation, on-the-fly subsetting |
| 1.4 Biogeochemical forcing | [todo] | also add `mlotst` mixed layer and `so` salinity, which the reader ignores today | optional prescribed nitrate/chlorophyll/oxygen from reanalysis instead of the internal NPZD |
| 1.5 Regional configs | [todo] | ready-made namelists for a shelf sea and an upwelling system |

## 2. Biogeochemistry (oxygen, carbon, export)

| Step | Status | Notes |
|---|---|---|
| 2.1 Oxygen tracer | [done] | fifth NPZD tracer: Garcia-Gordon solubility, Wanninkhof air-sea flux, production by photosynthesis, consumption by zooplankton, remineralisation (oxygen-limited) and fish respiration, deep resupply at the bed |
| 2.2 Metabolic index | [done] | Deutsch-style phi(O2, T, W) per species (`A_o`, `E_o`, `eps_o`, `phi_crit`, `m_hypoxia`); scales aerobic scope for feeding and adds suffocation mortality below phi = 1 |
| 2.3 Carbon tracers | [wip] | C:N stoichiometry and export diagnostics in carbon units are in; prognostic DIC and alkalinity (and therefore pCO2 and acidification) are not |
| 2.4 Fish-mediated export | [wip] | fish faeces, carcasses and respiration released below the reference depth are diagnosed separately from sinking detritus; faecal pellets do not yet have their own tracer and sinking speed |
| 2.5 Biological-pump diagnostics | [wip] | export by pathway is written to `series.nc` and plotted; sequestration depth and body-size sensitivity experiments remain |

## 3. Genetics (multilocus, mutation, G-matrix)

| Step | Status | Notes |
|---|---|---|
| 3.1 Infinitesimal model | [done] | midparent + segregation, exact offspring moments for biomass cells, validated against h2 |
| 3.2 Multilocus architecture | [done] | `genetics.loci` unlinked diploid loci; agents carry allele arrays, meiosis is free recombination. Verified to reproduce the infinitesimal model (midparent slope 1.0, segregation variance h2/2). Biomass classes keep the moment description and are resampled on disaggregation |
| 3.3 Mutation | [done] | `mutation_rate` per locus per gamete and `mutation_sd`; validated to add genetic variance across generations |
| 3.4 Trait covariance (G-matrix) | [done] | a pleiotropy matrix built from `[genetics.correlations]` reproduces requested genetic correlations within 0.02; realised variance reported each step. Note the correlations survive aggregation only approximately: biomass classes store one moment per trait and siblings are re-correlated on disaggregation |
| 3.5 Gene flow and adaptation diagnostics | [wip] | F_ST between latitude bands and selection differentials by cause of death are computed every step and written to `series.nc`. Breeder's equation vs realised response over generations is not yet automated |

## 4. Calibration (ABC, MCMC, emulator)

| Step | Status | Notes |
|---|---|---|
| 4.1 Ensemble sweep engine | [done] | grid / OAT / LHS / random, any parameter, parallel, `ensemble.nc` |
| 4.2 Observation interface | [done] | `targets.toml`: one `[[target]]` per observation (biomass, catch, chlorophyll, export, metabolic index) with species, day range, value and uncertainty. Readers for assessment/survey/FishMIP file formats still to write; values are entered by hand today |
| 4.3 Loss functions | [done] | weighted Gaussian log-likelihood, log-space for quantities spanning orders of magnitude; validated as 0 at the observation and 1.0 at one sigma |
| 4.4 ABC rejection and sequential Monte Carlo | [wip] | `--method abc` (rejection on a Latin-hypercube prior) run end to end; `--method smc` is implemented and unit-tested but has not yet been run on a full problem |
| 4.5 Gaussian-process emulator | [done] | RBF GP with marginal-likelihood hyperparameters and leave-one-out R2, then random-walk MCMC on the surface. Verified to reproduce a known surface (R2 0.999) and a known posterior. Active learning (proposing new members) not yet added |
| 4.6 Posterior ensembles | [wip] | `posterior.nc` holds the posterior sample, every member's parameters, losses and predictions; re-running the posterior as a predictive ensemble is not yet automated |

## 5. Validation (skill metrics, hindcasts)

| Step | Status | Notes |
|---|---|---|
| 5.1 Internal consistency suite | [done] | 25 checks: conservation, transport, genetics, aggregation, behaviour, shallow water, ensembles |
| 5.2 Skill metrics | [done] | `skill.py`: bias, RMSE, centred RMS, correlation, sd ratio, Nash-Sutcliffe for any series, plus a Taylor diagram; ensemble CRPS and rank histograms verified against analytic values |
| 5.3 Size spectra | [wip] | slope and intercept of the community abundance spectrum, from agent masses where available. The stage-structure fallback (used when agent snapshots are off) gives an implausible positive slope and should not be trusted |
| 5.4 Event hindcasts | [wip] | heatwave and collapse detectors validated on synthetic records and matched against observed event dates; real cases still need real forcing (area 1.3) |
| 5.5 Skill report | [done] | `python skill.py <run> observations.toml` writes skill.md, skill.nc and a Taylor diagram |

---

## Change log

- **This delivery:** prescribed ranges removed; biogeography left to the model. Species files no longer carry
  `native_range`, `bottom_max_m` or `spawning_grounds`. In their place: natal homing (ripe fish return to
  their birth place), emergent spawning aggregation, mate finding as a natural Allee effect, and a new
  evolvable cue for shallow water (`w_shelf`, a 13th trait with its own `seek_shelf` drive).
  Two real bugs surfaced. Biomass classes with no valid depth level over deep water produced no feeding unit
  and were silently zeroed when the fields were rebuilt (nitrogen error 3e-7); such classes are now carried
  over intact. And with fish free to spread globally, the feeding step built ~1.2M interacting units at once
  and was killed by memory; classes holding a negligible share of their species (`min_class_frac`) are now
  skipped and carried over instead, cutting that to 150k with conservation exact.
  A 45-day global run holds 117,727 agents over 2,870 individual-based cells (30% of the ocean) at ~7 s/step
  with nitrogen conserved to 4e-17.
  *Open:* 45 days is far too short to judge whether basin-scale stocks emerge. Occupancy is still broad
  (sardine 7,451 of 9,613 columns, flatfish 2,272), and `w_shelf` has not moved from its starting value in
  any species, as expected when no generation has turned over. Multi-year runs are needed before any claim
  that ranges self-organise.

- **This delivery:** realistic distributions and more individual cells.
  *Ranges.* Species now carry `native_range` boxes and a `bottom_max_m` depth limit, applied to habitat,
  movement, drift and agent dispersal. Cod occupy 154 columns of North Atlantic shelf, flatfish 68, sardine
  431 (temperate and upwelling seas), tuna 4,661 (Atlantic and Pacific), instead of every species forming
  continuous zonal bands. Enforcing this exposed a conservation leak: fish advected outside their range
  landed in cells with no valid habitat and were silently dropped when the fields were rebuilt from feeding
  units (nitrogen error 1e-6). Transport and agent drift now stop at the range edge, and the budget is back
  to 1e-16.
  *Coverage.* With `dynamic_fraction = 0.25`, `agents_per_class = 8` and a 200k agent budget, 27% of ocean
  columns are resolved as individuals (2,585 of 9,613) carrying ~33k super-individuals at ~2.5 s/step.
  Agent effort now tracks biomass: sardine 63% of agents for 70% of biomass, flatfish 4.5% for 4.3%
  (previously 55% for 5.5%). Tuna are over-sampled (26% for 6%) because their range spans 4,661 columns,
  and cod under-sampled (6% for 19%) because theirs spans 154.
  *Open:* range-restricted stocks decline faster (0.42-0.71x over 60 days) than the unrestricted global
  configuration, since the same biomass is packed into a fraction of the ocean. The mortality fit in
  `global.toml` predates the ranges and should be redone against them.

- **Patch:** `softmax` returned NaN for any row with no legal option (a cell whose every neighbour or depth
  level is blocked, which global land geometry produces); it now returns zeros, so those fish simply stay put.
  The symptom was a RuntimeWarning, but NaNs could reach the vertical distribution and the biomass fields.
  `--ensemble N` now sets the member count from the command line instead of being ignored.

- **This delivery:** the two global problems flagged last time.
  *Agent placement.* The front detector filled land with the ocean-mean temperature before differencing, so
  every coastline manufactured a fake front; land is now filled from the nearest ocean cell. The agent budget
  is also shared between species in proportion to their biomass (`budget_by_biomass`), filling each quota with
  the largest classes first. Flatfish fell from 55% of agents (5.5% of biomass) to 24%, sardine rose from 16%
  to 36% (70% of biomass), tuna now matches its share (7.9% vs 7.0%), and the coastal excess fell from 26% to
  22% of agents against 10% of biomass. Sardine remain under-sampled because their biomass sits in fewer,
  denser classes than the budget can subdivide.
  *Global biomass.* An ABC calibration on mortality (12 members, coarse global grid) converged but could not
  get within 10-30x of observed stocks, which was itself the finding: the initial condition was the problem,
  not the parameters. Species files carry biomass in g/m2 tuned for a productive basin, which spread over the
  world ocean gives ~7e8 t of sardine before the first step. Stocks are now initialised from observed global
  biomass (sardine 0.284 g/m2, tuna 0.023, cod 0.043, flatfish 0.017) with the calibrated mortalities, and a
  90-day run holds them at 0.60-0.85x of observed instead of 50x above.

- **This delivery:** global configuration (`global.toml`). Real coastlines from an offline land mask, depth
  from distance offshore, an equator-symmetric seasonal temperature field, and a wind-driven circulation:
  Sverdrup transport integrated westward from each basin's eastern boundary and closed by the least-squares
  projector (which holds the streamfunction at zero on every coast), plus a circumpolar term. Verified
  signs: westward subtropical interiors in all basins, northward Gulf Stream, southward Brazil Current,
  ACC 0.31 m/s, peak speed 0.66 m/s. Global runs conserve nitrogen to 1e-17 at 9,613 ocean columns
  (~2.5 s/step biomass-only, 3.5 s/step with 23,500 agents). Three attempts were needed on the streamfunction
  scaling and sign: depth-integrated Sverdrup transport is not the same quantity as a per-unit-depth
  streamfunction (100x error), and a hand-rolled western boundary layer produced 3 m/s spikes before the
  projector took over the closure. `fields.nc` now stores the land mask and bottom depth.
  **Still open:** agents remain shelf-heavy (flatfish 52% of the budget, down from 53%) because demersal
  biomass sits on shelves and the front detector finds its sharpest gradients at coastlines; and global
  biomass is still ~50x observed because species parameters were tuned for one basin and have not been
  recalibrated globally.

- **This delivery:** skill assessment (5.2-5.5) in `skill.py`, plus the area-4 documentation. Metrics, CRPS
  and rank histograms are checked against analytic values; the heatwave and collapse detectors are checked
  against synthetic records with known events. On a test hindcast with an imposed heatwave, the detector
  matched the observed event date (day 144 against 150). Two findings: a harmonic climatology fitted to a
  run shorter than two years absorbs much of the event it is meant to reveal (hence `--climatology`, which
  takes the seasonal cycle from a control run), and the size-spectrum fallback used when agent snapshots are
  off returns a positive slope, which is wrong. Validation: 52 checks.

- **This delivery:** calibration (4.2-4.6) in `calibrate.py`. A synthetic-truth test (targets generated from a
  run with known parameters) recovered all four parameters: the emulator posterior medians were 0.0059, 0.0065,
  1.599 and 0.095 against truths of 0.006, 0.005, 1.6 and 0.08, with every target matched inside its
  uncertainty. ABC recovered three of four, missing cod mortality because that test ran 60-day members against
  90-day targets. One bug found: relative uncertainties were divided by the observation twice in log space,
  giving losses of order 1e30. Validation: 48 checks.

- **This delivery:** multilocus genetics (3.2-3.4) and adaptation diagnostics (3.5). Traits are now built from
  unlinked diploid loci with a pleiotropy matrix, so genetic correlations are a property of the architecture
  rather than an imposed covariance, and mutation feeds variance back in. F_ST and per-cause selection
  differentials are reported each step: thermal deaths select on temperature optimum by up to 0.6 SD, which is
  the kind of signal the evolutionary-rescue experiments need. Two bugs found on the way: fish disaggregated
  from biomass initially had no genetic correlations at all, and selection differentials read zero because
  most mortality shrinks a super-individual rather than removing it. Validation: 43 checks.

- **This delivery:** oxygen and the metabolic index (2.1, 2.2) plus first carbon export diagnostics (2.3-2.5).
  Oxygen is a full tracer with realistic solubility (within 2% of published values), air-sea exchange, and
  biological sources and sinks including fish respiration; remineralisation slows as oxygen runs out. Each
  species has a metabolic index that scales its aerobic scope for feeding and kills it below phi = 1, giving
  hypoxia as a new cause of death. Carbon export past a reference depth is split into sinking detritus,
  fish faeces and carcasses, and fish respiration. Two sign/unit errors were caught by the new checks: the
  solubility polynomial returns umol/kg (not ml/l, a 50x error) and the index's temperature term had the
  wrong sign, making fish *more* tolerant when warmed. Validation: 38 checks.

- **This delivery:** CMEMS/GLORYS reader (1.3). Two scripts: `get_forcing.py` (real download, needs a
  Copernicus account) and `make_test_forcing.py` (offline lookalike). Testing against the lookalike found
  three real bugs: bathymetry files with their own coordinate names were not normalised, land was not being
  written as `_FillValue` in the generator, and the generator's velocity scaling was wrong. Currents lose
  about 23% of their mean speed to the non-divergent projection, and narrow boundary jets are smoothed by
  interpolation onto a coarser grid; domain-mean speed is resolution independent. Validation: 32 checks.

- **This delivery:** bathymetry and shelves (1.2). Partial bottom cells carry real volumes, sinking and deep
  nutrient restoring act at the sea bed, vertical mixing stops at the bottom, demersal species live on the bed,
  and fish cannot occupy dry cells. Horizontal transport was rebuilt: volume fluxes now come from a per-level
  streamfunction masked by the topography, which is non-divergent by construction, so a uniform tracer stays
  uniform over a shelf (previously it did not, and nitrogen drifted ~1%). Validation: 30 checks, all passing.

## Working notes

- Each numbered area gets its own module where it does not fit in `fishnet.py`: `calibrate.py` (area 4)
  and `skill.py` (area 5). Genetics and biogeochemistry stay in `fishnet.py` next to the code they touch.
- Every new mechanism needs a validation check before it is marked done.
- Nitrogen conservation to round-off is a hard invariant; carbon and oxygen budgets get the same treatment.
