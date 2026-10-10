#!/usr/bin/env python3
"""Compare a fishNET run with the observations in observations/, on the model's own grid.

    python obsval.py runs/global_clim_5yr_gui_2 [--window 365] [--rebuild]

The raw sources measure different things from the model, so each is first converted to
something the model can be compared with:

  OBIS      occurrence records (presence only, effort biased). Binned onto the model grid and
            compared as presence: does the model put the species where it has been recorded?
            AUC of model biomass against record/no-record cells, rank correlation of biomass
            with record counts, and the latitude of the occupied range.
  RAM       stock assessments (tonnes per stock and year). Averaged over OBS_YEARS, summed over
            each species' current stocks, and set against model biomass summed over the FAO
            areas those stocks sit in: total biomass against juveniles + adults, SSB against
            adults. Reported catch against model catch.
  ICCAT     tuna and swordfish catch on 5 degree squares. Averaged over OBS_YEARS and compared
            by rank with model biomass aggregated to the same squares (catch is effort
            weighted, so ranks rather than magnitudes).

The model side is a time mean over the last `window` days of fields.nc (the spin-up is left
out). Observations are processed once and cached in observations/processed/.
Needs numpy and netCDF4 only.
"""
import argparse, csv, io, json, re, sys, time, zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OBS = ROOT / "observations"
CACHE = OBS / "processed" / "validation_products.json"
VERSION = 7
OBS_YEARS = (1993, 2021)                    # the GLORYS years the climatology forcing is built from
R_EARTH = 6.371e6

# fishNET species -> scientific names (first is the OBIS/GBIF name; the rest are synonyms RAM uses)
SPECIES = {
    "anchoveta": ["Engraulis ringens"], "sardine_eu": ["Sardina pilchardus"],
    "sardine_pacific": ["Sardinops sagax"], "herring": ["Clupea harengus"],
    "capelin": ["Mallotus villosus"], "cod_atlantic": ["Gadus morhua"],
    "pollock": ["Gadus chalcogrammus", "Theragra chalcogramma"], "hake": ["Merluccius merluccius"],
    "skipjack": ["Katsuwonus pelamis"], "yellowfin": ["Thunnus albacares"],
    "albacore": ["Thunnus alalunga"], "bluefin_atlantic": ["Thunnus thynnus"],
    "bluefin_pacific": ["Thunnus orientalis"], "swordfish": ["Xiphias gladius"],
    "halibut": ["Hippoglossus hippoglossus", "Hippoglossus stenolepis"], "greenland_halibut": ["Reinhardtius hippoglossoides"],
    "flounder": ["Pleuronectes platessa"], "skate": ["Raja clavata"],
    "lanternfish": ["Benthosema glaciale"], "bristlemouth": ["Cyclothone braueri"],
    "polar_cod": ["Boreogadus saida"], "silverfish": ["Pleuragramma antarcticum", "Pleuragramma antarctica"],
    "toothfish": ["Dissostichus mawsoni"],
    "sardinella": ["Sardinella longiceps"], "anchovy": ["Engraulis encrasicolus"],
    "jumbo_squid": ["Dosidicus gigas"], "flying_squid": ["Todarodes pacificus"],
    "shortfin_squid": ["Illex argentinus", "Illex illecebrosus"],
    "blue_shark": ["Prionace glauca"], "shortfin_mako": ["Isurus oxyrinchus"],
    "spiny_dogfish": ["Squalus acanthias", "Squalus suckleyi"],
    "blue_whale": ["Balaenoptera musculus"], "fin_whale": ["Balaenoptera physalus"],
    "humpback_whale": ["Megaptera novaeangliae"], "minke_whale": ["Balaenoptera acutorostrata", "Balaenoptera bonaerensis"],
    "sperm_whale": ["Physeter macrocephalus"], "killer_whale": ["Orcinus orca"],
    "common_dolphin": ["Delphinus delphis"], "bottlenose_dolphin": ["Tursiops truncatus"],
    # the old lumped categories, kept so archived configs still find a starting biomass; their scientific names
    # now belong to the representative species above
    "squid": [], "shark": [], "whale": [], "dolphin": [],
}
# Model species that stand for a whole genus or family. OBIS is read at that level (grid3_<taxon>.geojson,
# fetched by observations/fetch_observations.py), and RAM stocks of the genus/family count toward them,
# except stocks of species another model species already claims (e.g. halibut within Pleuronectidae).
TAXA = {"hake": ("genus", "Merluccius"), "flounder": ("family", "Pleuronectidae"), "skate": ("family", "Rajidae"),
        "lanternfish": ("family", "Myctophidae"), "bristlemouth": ("genus", "Cyclothone"),
        "sardinella": ("genus", "Sardinella"), "anchovy": ("genus", "Engraulis")}
REF_OCEAN_AREA = 3.6e14                    # m2: the namelists' biomass (g m-2) = global tonnes x 1e6 / this
FISHING_SCALE = 0.5                        # RAM Legacy F is deliberately halved for now (see ram_fishing.toml)

# Biomass estimates for species without stock assessments (tonnes). Assessed species take the sum of
# their RAM Legacy stocks instead (see tables()).
LITERATURE = {
    "lanternfish": dict(biomass_t=1.0e9, low_t=6e8, high_t=1e10, coverage="global",
        source="Myctophidae ~1 Gt from net surveys (Gjosaeter & Kawaguchi 1980, FAO Fish. Tech. Pap. 193). Acoustic "
               "surveys put all mesopelagic fish at 2-16 Gt (Irigoien et al. 2014, Nat. Commun. 5:3271, 11-15 Gt; "
               "Proud et al. 2019, ICES J. Mar. Sci. 76:718). The conservative net-based value is used."),
    "bristlemouth": dict(biomass_t=3.0e8, low_t=3e7, high_t=3e9, coverage="global, derived",
        source="No direct biomass estimate exists. Cyclothone numbers are put at 1e14-1e16 fish, over half of all "
               "vertebrates at 100-1000 m; at ~0.3 g each that is 0.03-3 Gt, and the geometric middle is used."),
    "polar_cod": dict(biomass_t=1.0e6, low_t=3.5e5, high_t=2.0e6, coverage="Barents Sea only (lower bound)",
        source="Barents Sea acoustic ecosystem survey: ~2 Mt in 2000, 0.36 Mt in 2017 (Institute of Marine "
               "Research / MOSJ indicator; Aune et al. 2021, Prog. Oceanogr.). Mid value used."),
    "sardinella": dict(biomass_t=3.0e6, low_t=1.5e6, high_t=6.0e6, coverage="global, partly derived",
        source="NW African round sardinella 1.17-1.43 Mt from R/V Dr Fridtjof Nansen acoustic surveys 2002-03 "
               "(FAO/CECAF); Indian oil sardine catches ~0.57 Mt/yr (Frontiers Mar. Sci. 2018, 5:443), about 1 Mt of "
               "fish at F ~ 0.5; SE Asian, Gulf of Guinea and Brazilian sardinellas taken as roughly as much again."),
    "silverfish": dict(biomass_t=5.92e5, low_t=3.26e5, high_t=8.66e5, coverage="Ross Sea shelf only (lower bound)",
        source="Acoustic survey of the western Ross Sea shelf, Feb-Mar 2008: 592 kt, 95% CI 326-866 kt "
               "(O'Driscoll et al. 2011, Deep-Sea Res. II 58:181)."),
    # Representative cephalopods, sharks and cetaceans (replacing the lumped squid/shark/whale/dolphin below).
    # Cetacean biomass is abundance x a population-mean body mass (calves included), so it carries both errors.
    "jumbo_squid": dict(biomass_t=2.0e6, low_t=5e5, high_t=5e6, coverage="eastern Pacific, derived",
        source="IMARPE acoustic surveys off Peru, 1999-2015: 0.1-1.7 Mt, usually 0.5-0.8 Mt in summer (Inf. Inst. Mar "
               "Peru 43(1)); the species also ranges off Chile, Mexico and offshore, and its fishery lands 0.5-1 Mt a "
               "year, so ~2 Mt over the whole range."),
    "shortfin_squid": dict(biomass_t=1.5e6, low_t=7.4e5, high_t=1.8e6, coverage="SW Atlantic Illex argentinus",
        source="Illex argentinus 1.32-1.80 Mt from an environmentally dependent surplus-production assessment; B_MSY "
               "~0.74 Mt in other assessments (fisheryprogress.org stock assessment review). I. illecebrosus is small "
               "beside it (NW Atlantic catches ~10-25 kt)."),
    "shortfin_mako": dict(biomass_t=1.5e5, low_t=5e4, high_t=3e5, coverage="global, derived",
        source="No assessed total biomass. N Atlantic catches 3.6-4.75 kt/yr while overfished (ICCAT SCRS 2017, "
               "B2015/BMSY 0.57-0.85); other oceans land about as much again; at F ~0.1-0.15 (RAM N Pacific) that "
               "implies ~0.1-0.2 Mt."),
    "blue_whale": dict(biomass_t=1.2e6, low_t=7e5, high_t=2.25e6, coverage="global",
        source="10,000-25,000 animals (IUCN Red List 2018; Antarctic >2,000, Branch et al. 2007) x ~80 t mean mass."),
    "fin_whale": dict(biomass_t=5.0e6, low_t=3.5e6, high_t=7e6, coverage="global",
        source="~100,000 animals (IUCN Red List 2018) x ~50 t mean mass."),
    "humpback_whale": dict(biomass_t=3.8e6, low_t=2.5e6, high_t=4.5e6, coverage="global",
        source="~135,000 animals (IUCN Red List 2018) x ~28 t mean mass."),
    "minke_whale": dict(biomass_t=4.3e6, low_t=2.5e6, high_t=6e6, coverage="global",
        source="Antarctic minke ~515,000 (IWC 2012 circumpolar surveys) plus common minke ~200,000, x ~6 t mean mass."),
    "sperm_whale": dict(biomass_t=7.2e6, low_t=3e6, high_t=1.5e7, coverage="global",
        source="~360,000 animals (Whitehead 2002, Mar. Ecol. Prog. Ser. 242:295) x ~20 t mean mass (females ~15 t, "
               "males ~45 t)."),
    "killer_whale": dict(biomass_t=1.75e5, low_t=1.5e5, high_t=5e5, coverage="global, minimum",
        source="At least 50,000 animals (Forney & Wade 2006, in Whales, Whaling and Ocean Ecosystems) x ~3.5 t."),
    "common_dolphin": dict(biomass_t=4.8e5, low_t=3e5, high_t=8e5, coverage="global",
        source="~6 million animals, the most abundant cetacean (NAMMCO; Hammond et al. 2008) x ~80 kg."),
    "bottlenose_dolphin": dict(biomass_t=1.4e5, low_t=1e5, high_t=3e5, coverage="global, minimum",
        source="At least 600,000 animals (NAMMCO; Wells & Scott 2009) x ~230 kg."),
    "squid": dict(biomass_t=1.5e7, low_t=1e7, high_t=3e7, coverage="global",
        source="Global cephalopod biomass estimates (Rodhouse et al. 2014; Hunsicker et al. 2010)."),
    "shark": dict(biomass_t=3.0e6, low_t=1.5e6, high_t=6e6, coverage="global",
        source="Pelagic and coastal shark biomass synthesis (Dulvy et al. 2014; Worm et al. 2013)."),
    "whale": dict(biomass_t=8.0e6, low_t=5e6, high_t=1.2e7, coverage="global",
        source="Baleen whale global biomass synthesis (IWC estimates; Christensen 2006)."),
    "dolphin": dict(biomass_t=4.0e6, low_t=2e6, high_t=7e6, coverage="global",
        source="Delphinidae and small odontocete population estimates (Hammond et al. 2013)."),
}
ICCAT_CODES = {"ALB": "albacore", "BFT": "bluefin_atlantic", "SKJ": "skipjack", "YFT": "yellowfin", "SWO": "swordfish"}

# Ocean-wide tuna stocks list one FAO area in RAM but span a basin; their area sets by stock-id suffix.
BASIN = {"NATL": [21, 27, 31, 34], "SATL": [41, 47], "EATL": [27, 34, 47], "WATL": [21, 31, 41],
         "ATL": [21, 27, 31, 34, 41, 47], "EIO": [51, 57], "IO": [51, 57], "CWPAC": [61, 71, 81],
         "EPAC": [77, 87], "NPAC": [61, 67, 77], "SPAC": [81, 87], "MED": [37]}
BASIN_STOCKS = {"ATBTUNAEATL": [27, 34, 37], "ATBTUNAWATL": [21, 31], "PACBTUNA": [61, 67, 71, 77]}


# ------------------------------------------------------------------ FAO major fishing areas
def fao_area(lon, lat):
    """Approximate FAO major fishing area for cell centres. The real boundaries follow
    meridians and parallels except along coasts, so boxes plus a few coast rules are close
    at model resolution (land cells are masked anyway). Returns 0 where nothing matches."""
    lon = (np.asarray(lon, float) + 180) % 360 - 180
    lat = np.asarray(lat, float)
    lon, lat = np.broadcast_arrays(lon, lat)
    out = np.zeros(lon.shape, int)
    b = lambda lo0, lo1, la0, la1: (lon >= lo0) & (lon < lo1) & (lat >= la0) & (lat < la1)
    wpac = (lon >= 100) | (lon < -175)
    # Atlantic vs Pacific across the Americas
    pac_am = np.where(lat < 9, lon < -77, np.where(lat < 18, lon < -84, lon < -98))
    pac_am = np.where(lat < -50, lon < -68, np.where(lat < -5, lon < -70, pac_am))
    rules = [
        (37, b(-5.6, 42, 30, 47.5) & ~((lon < 0) & (lat > 43.3))),
        (18, (lat >= 66) & ~b(-80, 68.5, 66, 90) | b(-80, -42, 78, 90)),
        (21, b(-80, -42, 35, 78)),
        (27, b(-42, 68.5, 36, 90)),
        (31, b(-98, -40, 5, 35) & ~pac_am),
        (34, b(-40, 20, -6, 36)),
        (41, b(-70, -20, -60, 5) & ~pac_am),
        (47, b(-20, 30, -50, -6)),
        (48, b(-70, 30, -90, -50)),
        (58, b(30, 80, -90, -45) | b(80, 150, -90, -55)),
        (88, ((lon >= 150) | (lon < -70)) & (lat < -60)),
        (51, b(30, 80, -45, 31)),
        (57, b(80, 100, -55, 25) | (b(100, 150, -55, -8) & ((lon < 129) | (lat < -25)))),
        (71, wpac & (lat >= -25) & (lat < 20)),
        (61, wpac & (lat >= 20) & (lat < 66)),
        (67, b(-175, -100, 40, 66)),
        (81, ((lon >= 150) | (lon < -120)) & (lat >= -60) & (lat < -25)),
        (87, b(-120, -67, -60, -5) & pac_am),
        (77, b(-175, -70, -25, 40) & pac_am),
    ]
    for code, m in rules:
        out = np.where((out == 0) & m, code, out)
    return out


def stock_areas(stockid, fao, region):
    if stockid in BASIN_STOCKS:
        return BASIN_STOCKS[stockid]
    if region in ("Atlantic Ocean", "Pacific Ocean", "Indian Ocean", "Mediterranean-Black Sea"):
        for suf in sorted(BASIN, key=len, reverse=True):
            if stockid.endswith(suf):
                return BASIN[suf]
    return [int(x) for x in re.findall(r"\d+", str(fao))]


# ------------------------------------------------------------------ raw sources -> products
XNS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def xlsx_sheet(z, name):
    """Rows of one sheet as dicts keyed by the header row (no openpyxl needed)."""
    wb, rels = ET.fromstring(z.read("xl/workbook.xml")), ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    rid = next(s.get(RNS + "id") for s in wb.iter(XNS + "sheet") if s.get("name") == name)
    target = next(r.get("Target") for r in rels if r.get("Id") == rid).lstrip("/").removeprefix("xl/")
    if not hasattr(xlsx_sheet, "ss"):
        xlsx_sheet.ss = ["".join(t.text or "" for t in si.iter(XNS + "t"))
                         for si in ET.fromstring(z.read("xl/sharedStrings.xml")).iter(XNS + "si")]
    ss, head, rows = xlsx_sheet.ss, None, []
    for _, el in ET.iterparse(z.open("xl/" + target)):
        if el.tag != XNS + "row":
            continue
        r = {}
        for c in el.iter(XNS + "c"):
            v = c.find(XNS + "v")
            if v is not None:
                r[re.match(r"[A-Z]+", c.get("r"))[0]] = ss[int(v.text)] if c.get("t") == "s" else v.text
        el.clear()
        if head is None:
            head = r
        else:
            rows.append({head[k]: v for k, v in r.items() if k in head})
    return rows


def build_ram():
    src = OBS / "ram_legacy" / "RAMLDB v4.65.zip"
    if not src.exists():
        return {}
    with zipfile.ZipFile(src) as outer:
        name = next(n for n in outer.namelist() if n.startswith("Excel/") and n.endswith(".xlsx"))
        z = zipfile.ZipFile(io.BytesIO(outer.read(name)))
    if hasattr(xlsx_sheet, "ss"):
        del xlsx_sheet.ss
    sci = {s: k for k, names in SPECIES.items() for s in names}
    taxo = {r["scientificname"]: r for r in xlsx_sheet(z, "taxonomy")}

    def owner(name):                                    # which model species a RAM stock counts toward
        if name in sci:
            return sci[name]
        t = taxo.get(name, {})
        return next((sp for sp, (rank, taxon) in TAXA.items() if t.get(rank) == taxon), None)
    stocks = {r["stockid"]: r for r in xlsx_sheet(z, "stock")
              if owner(r.get("scientificname")) and r.get("state", "Current") == "Current"}
    units = {r["stockid"]: r for r in xlsx_sheet(z, "timeseries_units_views") if r.get("stockid") in stocks}
    vals, yearly = {}, {}
    y0, y1 = OBS_YEARS
    need = {"TBbest": "MT", "SSB": "MT", "TCbest": "MT", "ERbest": "ratio", "F": "1/yr"}
    for r in xlsx_sheet(z, "timeseries_values_views"):
        sid = r.get("stockid")
        if sid not in stocks or not (y0 <= int(float(r["year"])) <= y1):
            continue
        yr = int(float(r["year"]))
        for k, unit in need.items():
            if r.get(k) not in (None, "", "NA") and units.get(sid, {}).get(k) == unit:
                vals.setdefault(sid, {}).setdefault(k, []).append(float(r[k]))
                yearly.setdefault(sid, {}).setdefault(k, {})[yr] = float(r[k])
    out = {}
    for sid, s in stocks.items():
        v = vals.get(sid, {})
        rec = dict(id=sid, name=s.get("stocklong", sid), region=s.get("region", ""), fao=s.get("primary_FAOarea", ""),
                   sci=s["scientificname"], areas=stock_areas(sid, s.get("primary_FAOarea", ""), s.get("region", "")))
        for k, key in (("TBbest", "tb"), ("SSB", "ssb"), ("TCbest", "catch")):
            if len(v.get(k, [])) >= 3:
                rec[key], rec[key + "_years"] = float(np.mean(v[k])), len(v[k])
        # fishing mortality (1/yr): from the exploitation rate (catch / biomass) where there is one, else F
        yv = yearly.get(sid, {})
        if len(v.get("ERbest", [])) >= 3:
            rec["F"] = float(np.mean(-np.log(1 - np.minimum(v["ERbest"], 0.95))))
            rec["F_years"] = {y: round(float(-np.log(1 - min(e, 0.95))), 4) for y, e in yv["ERbest"].items()}
        elif len(v.get("F", [])) >= 3:
            rec["F"] = float(np.mean(v["F"]))
            rec["F_years"] = {y: round(f, 4) for y, f in yv["F"].items()}
        for k, key in (("TBbest", "tb"), ("SSB", "ssb"), ("TCbest", "catch")):   # annual series for validation
            if len(yv.get(k, {})) >= 3:
                rec[key + "_series"] = {y: round(x, 1) for y, x in sorted(yv[k].items())}
        out.setdefault(owner(s["scientificname"]), []).append(rec)
    return out


def build_obis():
    out = {}
    for sp in SPECIES:
        p = OBS / "obis" / f"grid3_{TAXA[sp][1]}.geojson" if sp in TAXA else None
        p = p if p is not None and p.exists() else OBS / "obis" / f"grid3_{sp}.geojson"
        if not p.exists():
            continue
        feats = json.loads(p.read_text()).get("features", [])
        pts = []
        for f in feats:
            ring = np.asarray(f["geometry"]["coordinates"][0], float)
            pts.append([round(float(ring[:, 0].min() + ring[:, 0].max()) / 2, 3),
                        round(float(ring[:, 1].min() + ring[:, 1].max()) / 2, 3), int(f["properties"].get("n", 0))])
        if pts:
            out[sp] = pts
    return out


def build_iccat():
    src = OBS / "iccat" / "cdis5024_all.zip"
    if not src.exists():
        return {}
    y0, y1 = OBS_YEARS
    tot = {}
    with zipfile.ZipFile(src) as z:
        with z.open(z.namelist()[0]) as f:
            for r in csv.DictReader(io.TextIOWrapper(f, encoding="utf-8", errors="replace")):
                sp = ICCAT_CODES.get(r["SpeciesCode"])
                if not sp or not (y0 <= int(r["YearC"]) <= y1):
                    continue
                key = (sp, float(r["xLon5ctoid"]), float(r["yLat5ctoid"]))
                tot[key] = tot.get(key, 0.0) + float(r["Catch_t"] or 0)
    out = {}
    for (sp, lon, lat), c in tot.items():
        if c > 0:
            out.setdefault(sp, []).append([lon, lat, float(f"{c / (y1 - y0 + 1):.4g}")])
    return out


def products(rebuild=False):
    """The processed observations, rebuilt when missing, stale or asked for."""
    if CACHE.exists() and not rebuild:
        try:
            d = json.loads(CACHE.read_text())
            if d.get("version") == VERSION:
                return d
        except (OSError, ValueError):
            pass
    t0 = time.time()
    d = dict(version=VERSION, years=list(OBS_YEARS), built=time.strftime("%Y-%m-%d %H:%M"),
             species={k: (f"{TAXA[k][1]} ({TAXA[k][0]})" if k in TAXA else v[0]) for k, v in SPECIES.items() if v},
             obis=build_obis(), ram=build_ram(), iccat=build_iccat())
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(d))
    print(f"observation products built in {time.time() - t0:.0f} s -> {CACHE}", file=sys.stderr, flush=True)
    return d


# ------------------------------------------------------------------ statistics
def rankdata(a):
    a = np.asarray(a, float)
    order = np.argsort(a, kind="mergesort")
    s = a[order]
    new = np.r_[True, s[1:] != s[:-1]]
    dense = np.cumsum(new)
    cnt = np.r_[np.flatnonzero(new), len(a)]
    r = np.empty(len(a))
    r[order] = 0.5 * (cnt[dense] + cnt[dense - 1] + 1)
    return r


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def pearson(x, y):
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def auc(score, label):
    """Probability that a presence cell scores higher than a background cell (ties count half)."""
    label = np.asarray(label, bool)
    n1, n0 = label.sum(), (~label).sum()
    if n1 == 0 or n0 == 0:
        return None
    r = rankdata(score)
    return float((r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def lat_stats(lat, w):
    w = np.asarray(w, float)
    if w.sum() <= 0:
        return None
    o = np.argsort(lat)
    c = np.cumsum(w[o]) / w.sum()
    q = lambda p: float(lat[o][min(np.searchsorted(c, p), len(o) - 1)])
    return dict(mean=float(np.average(lat, weights=w)), p5=q(0.05), p95=q(0.95))


def core_range(B, frac=0.95):
    """Cells holding `frac` of the biomass: the model's occupied range."""
    f = B.ravel()
    o = np.argsort(-f)
    c = np.cumsum(f[o])
    occ = np.zeros(f.size, bool)
    if c[-1] > 0:
        occ[o[: np.searchsorted(c, frac * c[-1]) + 1]] = True
    return occ.reshape(B.shape)


# ------------------------------------------------------------------ model side
class Grid:
    def __init__(self, lon, lat, mask):
        self.lon, self.lat, self.mask = np.asarray(lon, float), np.asarray(lat, float), np.asarray(mask, bool)
        self.ny, self.nx = self.mask.shape
        dlon = abs(self.lon[1] - self.lon[0]) if self.nx > 1 else 1.0
        dlat = abs(self.lat[1] - self.lat[0]) if self.ny > 1 else 1.0
        self.dlon, self.dlat = dlon, dlat
        s = np.sin(np.radians(np.clip(self.lat + dlat / 2, -90, 90))) - np.sin(np.radians(np.clip(self.lat - dlat / 2, -90, 90)))
        self.area = (R_EARTH ** 2 * np.radians(dlon) * s)[:, None] * np.ones(self.nx) * self.mask
        self.LON, self.LAT = np.meshgrid(self.lon, self.lat)
        self.fao = fao_area(self.LON, self.LAT) * self.mask
        self.glob = dlon * self.nx > 300

    def index(self, lon, lat):
        """(j, i, inside) of points on the grid; longitude wraps on a global grid."""
        lon = np.asarray(lon, float)
        l0 = self.lon[0] - self.dlon / 2
        x = (lon - l0) % 360 if self.glob else lon - l0
        i = np.floor(x / self.dlon).astype(int)
        j = np.floor((np.asarray(lat, float) - (self.lat[0] - self.dlat / 2)) / self.dlat).astype(int)
        ok = (i >= 0) & (i < self.nx) & (j >= 0) & (j < self.ny)
        return np.clip(j, 0, self.ny - 1), np.clip(i, 0, self.nx - 1), ok

    def bin(self, pts):
        """Sum point values onto the grid, ocean cells only."""
        out = np.zeros((self.ny, self.nx))
        if not len(pts):
            return out
        p = np.asarray(pts, float)
        j, i, ok = self.index(p[:, 0], p[:, 1])
        np.add.at(out, (j[ok], i[ok]), p[ok, 2])
        return out * self.mask


def model_means(F, S, window):
    """Time means over the last `window` days: juvenile+adult and adult biomass (g m-2) per species
    from fields.nc, and catch (t yr-1) from series.nc."""
    t = np.asarray(F["time"][:], float)
    if not len(t):
        return None
    k0 = int(np.searchsorted(t, t[-1] - window + 1e-6)) if window > 0 else 0
    stages = [s for s in F.getncattr("stages").split(",")] if "stages" in F.ncattrs() else []
    ia = stages.index("adult") if "adult" in stages else len(F.dimensions["stage"]) - 1
    ij = stages.index("juvenile") if "juvenile" in stages else ia
    fb = F["fish_biomass"]
    B = np.zeros((len(F.dimensions["species"]), 2) + F["mask"].shape)
    for k in range(k0, len(t)):                         # one frame at a time keeps memory small
        B[:, 1] += np.nan_to_num(np.asarray(fb[k, :, ia], float))
        if ij != ia:
            B[:, 0] += np.nan_to_num(np.asarray(fb[k, :, ij], float))
    B /= len(t) - k0
    catch = None
    if S is not None and "catch" in S.variables and len(S.dimensions["time"]):
        ts = np.asarray(S["time"][:], float)
        s0 = int(np.searchsorted(ts, ts[-1] - window + 1e-6)) if window > 0 else 0
        catch = np.nan_to_num(np.asarray(S["catch"][s0:], float)).mean(0) * 365 / 1e6
    return dict(juvad=B[:, 0] + (B[:, 1] if ij != ia else 0), adult=B[:, 1], catch=catch,
                t0=float(t[k0]), t1=float(t[-1]), frames=len(t) - k0)


# ------------------------------------------------------------------ comparison
def compare(grid, names, M, P):
    """Every observation-model pairing for every species that has observations."""
    rows = []
    rowocean = grid.mask.sum(1)
    for q, sp in enumerate(names):
        r = dict(name=sp, q=q, sci=P["species"].get(sp))
        B = M["juvad"][q] * grid.mask
        occ = core_range(B) & grid.mask

        pts = P["obis"].get(sp)
        if pts:
            O = grid.bin(pts)
            pres = (O > 0) & grid.mask
            m = grid.mask
            both = pres & occ
            both_pos = pres & (B > 0)
            lat_m = lat_stats(grid.LAT[occ], np.ones(occ.sum())) if occ.any() else None
            lat_o = lat_stats(grid.LAT[pres], np.ones(pres.sum())) if pres.any() else None
            with np.errstate(invalid="ignore", divide="ignore"):
                prof_m = np.where(rowocean > 0, occ.sum(1) / rowocean, np.nan)
                prof_o = np.where(rowocean > 0, pres.sum(1) / rowocean, np.nan)
            r["obis"] = dict(
                cells=int(pres.sum()), records=int(O[pres].sum()), outside=int(sum(p[2] for p in pts) - O.sum()),
                auc=auc(B[m], pres[m]), rho=spearman(B[both_pos], O[both_pos]), n_rho=int(both_pos.sum()),
                sens=float(both.sum() / pres.sum()) if pres.any() else None,
                prec=float(both.sum() / occ.sum()) if occ.any() else None,
                lat_model=lat_m, lat_obs=lat_o, model_cells=int(occ.sum()),
                profile=dict(lat=grid.lat.tolist(), model=np.round(prof_m, 4).tolist(), obs=np.round(prof_o, 4).tolist()),
                scatter=dict(obs=O[both_pos].round(1).tolist(), model=B[both_pos].astype("f4").tolist()))

        stocks = P["ram"].get(sp)
        if stocks:
            # whichever of total or spawning biomass covers more stocks (assessments report one or both)
            use_tb = sum("tb" in s for s in stocks) >= sum("ssb" in s for s in stocks)
            key = "tb" if use_tb else "ssb"
            used = [s for s in stocks if key in s]
            areas = sorted({a for s in used for a in s["areas"]})
            inside = np.isin(grid.fao, areas) & grid.mask
            field = M["juvad"][q] if use_tb else M["adult"][q]
            model_t = float((field * grid.area * inside).sum() / 1e6)
            caught = [s for s in stocks if "catch" in s]
            r["ram"] = dict(
                metric="total biomass" if use_tb else "spawning biomass", key=key, areas=areas,
                obs_t=float(sum(s[key] for s in used)) if used else None,
                model_t=model_t if used else None, model_domain_t=float((field * grid.area).sum() / 1e6),
                area_cells=int(inside.sum()), n_used=len(used), n_stocks=len(stocks),
                catch_obs=float(sum(s["catch"] for s in caught)) if caught else None, n_catch=len(caught),
                catch_model=float(M["catch"][q]) if M["catch"] is not None else None,
                stocks=sorted(stocks, key=lambda s: -s.get(key, s.get("catch", 0))))

        boxes = P["iccat"].get(sp)
        if boxes:
            b = np.asarray(boxes, float)
            jj, ii, ok = grid.index(grid.LON.ravel(), grid.LAT.ravel())
            kx = np.floor(grid.LON.ravel() / 5).astype(int)
            ky = np.floor(grid.LAT.ravel() / 5).astype(int)
            wet = grid.mask.ravel()
            Bm, cnt = {}, {}
            for x, y, v in zip(kx[wet], ky[wet], B.ravel()[wet]):
                Bm[(x, y)] = Bm.get((x, y), 0.0) + v
                cnt[(x, y)] = cnt.get((x, y), 0) + 1
            obs, mod = [], []
            for lon, lat, c in b:
                k = (int(np.floor(lon / 5)), int(np.floor(lat / 5)))
                if k in cnt:
                    obs.append(c); mod.append(Bm[k] / cnt[k])
            obs, mod = np.asarray(obs), np.asarray(mod)
            pos = (mod > 0) & (obs > 0)
            r["iccat"] = dict(boxes=len(obs), boxes_total=len(b), rho=spearman(mod, obs),
                              r_log=pearson(np.log10(mod[pos]), np.log10(obs[pos])) if pos.sum() > 2 else None,
                              catch_t=float(b[:, 2].sum()),
                              scatter=dict(obs=obs.round(2).tolist(), model=mod.astype("f4").tolist()))
        if len(r) > 3:
            rows.append(r)
    return rows


def yearly(F, grid, names, P, stride_days=30):
    """Model biomass per calendar year in each species' assessed FAO areas (total biomass -> juveniles + adults,
    spawning biomass -> adults, as in compare()), beside the RAM Legacy total of the stocks reporting that year.
    Reads one fields.nc snapshot per ~stride_days."""
    import datetime as dt
    t = np.asarray(F["time"][:], float)
    if not len(t):
        return {}
    unit = F["time"].getncattr("units") if "units" in F["time"].ncattrs() else "days since 2000-01-01"
    t0 = dt.datetime.fromisoformat(unit.split("since")[1].strip().replace(" ", "T")[:19])
    step = max(1, int(round(stride_days / max(np.median(np.diff(t)) if len(t) > 1 else 1, 1e-9))))
    frames = list(range(0, len(t), step))
    stages = F.getncattr("stages").split(",") if "stages" in F.ncattrs() else []
    ia = stages.index("adult") if "adult" in stages else len(F.dimensions["stage"]) - 1
    ij = stages.index("juvenile") if "juvenile" in stages else ia
    use = {}
    for q, sp in enumerate(names):
        st = P["ram"].get(sp)
        if not st:
            continue
        key = "tb" if sum("tb" in x for x in st) >= sum("ssb" in x for x in st) else "ssb"
        stocks = [x for x in st if x.get(key + "_series")]
        if stocks:
            areas = sorted({a for x in stocks for a in x["areas"]})
            use[q] = (sp, key, stocks, np.isin(grid.fao, areas) & grid.mask)
    if not use:
        return {}
    sums = {q: {} for q in use}
    fb = F["fish_biomass"]
    for k in frames:
        yr = (t0 + dt.timedelta(days=float(t[k]))).year
        for q, (sp, key, stocks, inside) in use.items():
            a = np.nan_to_num(np.asarray(fb[k, q, ia], float))
            if key == "tb" and ij != ia:
                a = a + np.nan_to_num(np.asarray(fb[k, q, ij], float))
            v = float((a * grid.area * inside).sum() / 1e6)
            sums[q].setdefault(yr, []).append(v)
    out = {}
    for q, (sp, key, stocks, inside) in use.items():
        years = sorted(sums[q])
        obs, nobs = [], []
        for y in years:
            vals = [x[key + "_series"].get(str(y), x[key + "_series"].get(y)) for x in stocks]
            vals = [v for v in vals if v is not None]
            obs.append(float(sum(vals)) if vals else None)
            nobs.append(len(vals))
        mod = [float(np.mean(sums[q][y])) for y in years]
        pair = [(a, b) for a, b in zip(mod, obs) if a and b and a > 0 and b > 0]
        r = None
        if len(pair) >= 4:                                 # year-to-year changes, the part a hindcast can get right
            dm, do = np.diff(np.log([a for a, _ in pair])), np.diff(np.log([b for _, b in pair]))
            r = pearson(dm, do)
        out[sp] = dict(years=years, model=mod, obs=obs, n_stocks=nobs, of=len(stocks), metric=key, r_change=r)
    return out


def obs_field(grid, P, sp, source):
    """An observation layer on the model grid for the maps: OBIS records per cell, or ICCAT mean annual
    catch spread evenly over the model cells of each 5 degree square."""
    if source == "iccat":
        out = np.zeros((grid.ny, grid.nx))
        b = P["iccat"].get(sp) or []
        kx, ky = np.floor(grid.LON / 5).astype(int), np.floor(grid.LAT / 5).astype(int)
        for lon, lat, c in b:
            m = (kx == int(np.floor(lon / 5))) & (ky == int(np.floor(lat / 5))) & grid.mask
            if m.any():
                out[m] += c / m.sum()
        return out
    return grid.bin(P["obis"].get(sp) or [])


def summary(rows):
    lines = []
    for r in rows:
        o, s, c = r.get("obis", {}), r.get("ram", {}), r.get("iccat", {})
        f = lambda v, p=2: "—" if v is None else f"{v:.{p}f}"
        ratio = f(s["model_t"] / s["obs_t"], 2) if s.get("obs_t") and s.get("model_t") is not None else "—"
        lines.append(f"  {r['name']:18s} OBIS AUC {f(o.get('auc'))}  rho {f(o.get('rho'))}  "
                     f"lat {f((o.get('lat_model') or {}).get('mean'), 0)}/{f((o.get('lat_obs') or {}).get('mean'), 0)}  "
                     f"RAM model/obs {ratio}  ICCAT rho {f(c.get('rho'))}")
    return "\n".join(lines)


def tables(P, root=OBS):
    """Write literature_biomass.toml (starting biomass per species, for run.initial_biomass = "literature")
    and ram_fishing.toml (fishing mortality per species from RAM Legacy, before and after FISHING_SCALE)."""
    y0, y1 = P["years"]
    bio, fish = {}, {}
    for sp in SPECIES:
        st = P["ram"].get(sp, [])
        tb = [s for s in st if "tb" in s]
        ssb = [s for s in st if "tb" not in s and "ssb" in s]
        if tb or ssb:
            tot = sum(s["tb"] for s in tb) + sum(s["ssb"] for s in ssb)
            bio[sp] = dict(biomass_t=float(f"{tot:.4g}"), coverage="assessed stocks only (lower bound)",
                           source=f"RAM Legacy v4.65, {y0}-{y1} mean of {len(tb)} stocks' total biomass (TBbest)"
                                  + (f" and {len(ssb)} stocks' spawning biomass where no total is given" if ssb else "")
                                  + f", of {len(st)} current stocks")
        if sp in LITERATURE:
            bio[sp] = dict(LITERATURE[sp])
        f = [(s["F"], s.get("tb", s.get("ssb", 3 * s.get("catch", 0)))) for s in st if "F" in s]
        f = [(x, w) for x, w in f if w > 0]
        if f:
            F = float(np.average([x for x, _ in f], weights=[w for _, w in f]))
            fish[sp] = dict(F_ram=round(F, 4), F_model=round(FISHING_SCALE * F, 4), stocks=len(f), of=len(st),
                            min=round(min(x for x, _ in f), 3), max=round(max(x for x, _ in f), 3))
    q = lambda v: json.dumps(v) if isinstance(v, str) else repr(v)
    lines = ["# Starting biomass per fishNET species, used when [run] initial_biomass = \"literature\".",
             f"# Written by `python obsval.py --tables` on {time.strftime('%Y-%m-%d')}. biomass_t is the global total in",
             f"# tonnes; the model spreads it as biomass_t x 1e6 / {REF_OCEAN_AREA:.2g} m2, the namelists' convention.",
             "# Assessed species: the sum of their RAM Legacy stocks, which misses unassessed stocks (a lower bound).",
             "# Validating against RAM after starting from RAM is partly circular: compare the end state, not the start.", ""]
    for sp in SPECIES:
        if sp in bio:
            lines += [f"[{sp}]"] + [f"{k} = {q(v)}" for k, v in bio[sp].items()] + [""]
    (root / "literature_biomass.toml").write_text("\n".join(lines))
    lines = ["# Fishing mortality per fishNET species from RAM Legacy v4.65 stock assessments.",
             f"# Written by `python obsval.py --tables` on {time.strftime('%Y-%m-%d')}.",
             f"# F_ram: {y0}-{y1} mean of -ln(1 - catch/biomass) per stock (or the assessed F where there is no",
             "#        exploitation rate), averaged over the species' stocks weighted by their biomass. 1/yr.",
             f"# F_model = {FISHING_SCALE} x F_ram.  NOTE: RAM fishing is DELIBERATELY SCALED DOWN by {FISHING_SCALE} for now,",
             "#        so the model is fished at half the assessed rate. F_model is what the namelists use (fishing_F).", ""]
    for sp in SPECIES:
        if sp in fish:
            lines += [f"[{sp}]"] + [f"{k} = {q(v)}" for k, v in fish[sp].items()] + [""]
    (root / "ram_fishing.toml").write_text("\n".join(lines))
    # year-by-year fishing mortality per species and FAO area, for [fishing] mode = "ram"
    series = {}
    for sp in SPECIES:
        st = [x for x in P["ram"].get(sp, []) if x.get("F_years")]
        if not st:
            continue
        w = {x["id"]: x.get("tb", x.get("ssb", 3 * x.get("catch", 0))) or 1.0 for x in st}
        areas = {}
        for a in sorted({a for x in st for a in x["areas"]}):
            here = [x for x in st if a in x["areas"]]
            years = sorted({int(y) for x in here for y in x["F_years"]})
            yf = {}
            for y in years:
                have = [(x["F_years"][str(y)] if str(y) in x["F_years"] else x["F_years"].get(y), w[x["id"]]) for x in here]
                have = [(f, ww) for f, ww in have if f is not None]
                if have:
                    yf[str(y)] = round(float(np.average([f for f, _ in have], weights=[ww for _, ww in have])), 4)
            if yf:
                areas[str(a)] = dict(years=yf, mean=round(float(np.mean(list(yf.values()))), 4),
                                     stocks=[x["id"] for x in here])
        series[sp] = dict(areas=areas, mean=fish.get(sp, {}).get("F_ram"))
    (root / "ram_fishing_series.json").write_text(json.dumps(dict(
        note=f"Fishing mortality (1/yr) from RAM Legacy v4.65 per species, FAO major area and year ({y0}-{y1}): "
             "-ln(1 - catch/biomass) per stock (or the assessed F), averaged over the stocks assessed in that area, "
             "weighted by their mean biomass. Unscaled: the namelist's [fishing] scale multiplies it.",
        years=[y0, y1], species=series), indent=1))
    return bio, fish


def main():
    import netCDF4 as nc
    ap = argparse.ArgumentParser(description="compare a fishNET run with observations/")
    ap.add_argument("run", nargs="?")
    ap.add_argument("--tables", action="store_true", help="write observations/literature_biomass.toml and ram_fishing.toml")
    ap.add_argument("--window", type=float, default=365, help="days at the end of the run to average (0 = all)")
    ap.add_argument("--rebuild", action="store_true", help="re-process the raw observation files")
    a = ap.parse_args()
    P = products(a.rebuild)
    if a.tables:
        bio, fish = tables(P)
        print(f"wrote {OBS / 'literature_biomass.toml'} ({len(bio)} species) and {OBS / 'ram_fishing.toml'} ({len(fish)} species)")
    if not a.run:
        return
    run = Path(a.run)
    with nc.Dataset(run / "fields.nc") as F, nc.Dataset(run / "series.nc") as S:
        F.set_auto_mask(False); S.set_auto_mask(False)
        grid = Grid(F["lon"][:], F["lat"][:], np.asarray(F["mask"][:]) > 0.5)
        names = F.getncattr("species").split(",")
        M = model_means(F, S, a.window)
    print(f"{run.name}: model days {M['t0']:.0f}-{M['t1']:.0f} ({M['frames']} frames), observations {OBS_YEARS[0]}-{OBS_YEARS[1]}")
    print(summary(compare(grid, names, M, P)))


if __name__ == "__main__":
    main()
