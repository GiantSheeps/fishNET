fishNET validation data
=======================
Downloaded 2026-09-22 by fetch_observations.py (re-run it to refresh; --force re-downloads).
All of it is public and needed no account. Total on disk: 146 MB.
Nothing here is wired into the model yet - these are raw sources, not observations.toml entries.


ram_legacy/                                                              112 MB
  RAM Legacy Stock Assessment Database v4.65 (Zenodo doi:10.5281/zenodo.11995054,
  published 2024-06-17). ~500 assessed stocks worldwide with annual time series of
  spawning stock biomass, total biomass, recruitment, fishing mortality and catch.
  Archive verified intact. Useful members:
    RAMCore/*.csv                  ready-to-read series (Bio, Btrend, CPUE, BvB, CvMSY, ...)
    R Data/DBdata[asmt][v4.65].RData   the full relational database as R objects
    Excel/RAMLDB v4.65 (assessment data only).xlsx   same data, 63 MB spreadsheet
    Stock Summary Files (taxGroup)/  per-group PDF summaries
  Covers most fishNET species: cod, herring, capelin, pollock, hake, both sardines,
  anchoveta, halibut, Greenland halibut, flounders, skates, and the tunas.
  Caveat: these are stock-assessment model outputs, not raw observations.

ices_sag/                                                                3.1 MB
  ICES Stock Assessment Graphs, 2023 assessment round (sag.ices.dk).
    stocklist_2023.json            189 published stocks, with AssessmentKey per stock
    summary/<stock>_<key>.json     96 stocks matching fishNET species; 85 carry an SSB
                                   series. 4,279 stock-years spanning 1904-2024.
  Each summary file has a "Lines" array, one row per year, with Year, SSB, Low_SSB,
  High_SSB, Recruitment (+bounds), Catches, Landings, Discards. NOTE the capital L in
  "Lines". Low_SSB/High_SSB give a defensible sigma for a [[series]] entry instead of a
  guessed one. Units are in the file's "Units"/"StockSizeUnits" fields - check them, they
  are not uniform across stocks (tonnes vs thousand tonnes).
  Best represented: cod 17 stocks, herring 13, sole 12, plaice 9, haddock 6, whiting 5.

ices_datras/                                                              36 KB
  Index only, for the ICES bottom-trawl survey database (datras.ices.dk).
    survey_list.xml                all available surveys (NS-IBTS, BITS, EVHOE, ...)
    years_<survey>.xml             years available per survey for five key surveys
  The actual haul data is NOT here: it is served per survey/year/quarter and is large.
  Fetch it from the same SOAP service (getHHdata / getCAdata / getHLdata) once you know
  which survey-years you want. This is the closest thing to raw observation in this
  directory - haul-level, georeferenced, with length frequencies - and it is the right
  source for testing spatial distribution and range shifts on the shelf.

obis/                                                                    3.2 MB
  Ocean Biodiversity Information System (api.obis.org), per fishNET species.
    stats_<species>.json           record count, dataset count, year range
    grid3_<species>.geojson        occurrence counts binned into hexagonal cells
  6,422,801 records across the 23 species. The gridded files are the directly
  model-comparable product: counts per cell, straight against a gridded biomass field.
  Range: cod 3.1M records / 912 cells, herring 1.2M, yellowfin 247k / 2759 cells,
  down to bristlemouth 2,637 and anchoveta 141.
  NOTE: anchoveta has 141 records despite being the largest fishery on Earth - OBIS
  coverage reflects who submits survey data, not where fish are. Do not read low record
  counts as low abundance. Occurrence data is presence-only and effort-biased.

gbif/                                                                    140 KB
  Global Biodiversity Information Facility (api.gbif.org), per species.
    search_<species>.json          total occurrence count + year-frequency facet
  9,687,580 occurrences across the 23 species. Overlaps OBIS heavily (OBIS feeds GBIF).
  Kept for the year histograms, which show sampling effort over time - useful for
  deciding which decades an occurrence comparison is even meaningful in.

iccat/                                                                    28 MB
  International Commission for the Conservation of Atlantic Tunas (iccat.int/Data).
    t1nc_20260129.zip              Task 1 nominal catch: annual catch by species, gear,
                                   flag and area, 1950-present
    cdis5024_bySpecies.zip         catch distribution 1950-2024 GRIDDED IN SPACE, split
                                   by species - the spatially explicit product, and the
                                   best match for a spatial model
    cdis5024_all.zip               same, not split by species
    EFFDIS_LL2000-2024.csv         longline fishing effort distribution 2000-2024
    SCRScatalogues_MainSpecies_1995-2024.xlsx   data availability catalogue by stock
  Covers swordfish and the Atlantic tunas (bluefin, yellowfin, skipjack, albacore).
  Pacific tunas are NOT here - those live with IATTC and WCPFC.


Not downloaded, and why
-----------------------
  NOAA StockSMART   JavaScript app, CSV export by hand; no public REST endpoint found.
  FAO FishStat      portal is JS-rendered with no bulk link exposed; use FishStatJ.
  CCAMLR            bulletin is public, but toothfish fine-scale data needs a request.
  AquaMaps          download form, no API.
  IMARPE            anchoveta acoustic surveys are PDF reports only.
  FishBase          fishbase.ropensci.org does not resolve from here. Life-history
                    parameters (L_inf, L_mat, length-weight a/b) would validate the
                    species files directly, so this one is worth chasing via the
                    rfishbase R package or a FishBase dump.
  Mesopelagics      no database exists for lanternfish or bristlemouth biomass. Only
                    per-paper supplements. Those two stay order-of-magnitude checks.


Before using any of this
------------------------
  1. Scope. Assessments report a stock in a management area; fishNET reports biomass over
     its whole domain. "Atlantic cod" is ~17 ICES stocks plus others. Sum or subset to
     match, or the comparison is meaningless.
  2. Forcing. A looped climatology has no interannual variability, so it cannot reproduce
     a specific year's SSB or a dated collapse. Interannual validation needs the real
     GLORYS 1993-2021 sequence, not forcing_clim_025/.
  3. Units. Check StockSizeUnits per ICES stock; RAM units are documented per series in
     its DB Documents/ sheets.
