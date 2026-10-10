#!/usr/bin/env python3
"""Fetch the open fisheries datasets fishNET can be validated against.

Everything here is public and needs no account. Re-run to refresh; existing files are kept
unless --force. Nothing in this script touches the model.

    python fetch_observations.py [--force] [--skip ram,obis,...]
"""
import argparse, json, sys, time, urllib.parse, urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
UA = {"User-Agent": "fishNET-validation/1.0 (research use)"}

# fishNET species -> the scientific name used by OBIS/GBIF
SPECIES = {
    "anchoveta": "Engraulis ringens",          "sardine_eu": "Sardina pilchardus",
    "sardine_pacific": "Sardinops sagax",      "herring": "Clupea harengus",
    "capelin": "Mallotus villosus",            "cod_atlantic": "Gadus morhua",
    "pollock": "Gadus chalcogrammus",          "hake": "Merluccius merluccius",
    "skipjack": "Katsuwonus pelamis",          "yellowfin": "Thunnus albacares",
    "albacore": "Thunnus alalunga",            "bluefin_atlantic": "Thunnus thynnus",
    "bluefin_pacific": "Thunnus orientalis",   "swordfish": "Xiphias gladius",
    "halibut": "Hippoglossus hippoglossus",    "greenland_halibut": "Reinhardtius hippoglossoides",
    "flounder": "Pleuronectes platessa",       "skate": "Raja clavata",
    "lanternfish": "Benthosema glaciale",      "bristlemouth": "Cyclothone braueri",
    "polar_cod": "Boreogadus saida",           "silverfish": "Pleuragramma antarcticum",
    "toothfish": "Dissostichus mawsoni",
    "sardinella": "Sardinella longiceps",      "anchovy": "Engraulis encrasicolus",
    "jumbo_squid": "Dosidicus gigas",          "flying_squid": "Todarodes pacificus",
    "shortfin_squid": "Illex argentinus",      "blue_shark": "Prionace glauca",
    "shortfin_mako": "Isurus oxyrinchus",      "spiny_dogfish": "Squalus acanthias",
    "blue_whale": "Balaenoptera musculus",     "fin_whale": "Balaenoptera physalus",
    "humpback_whale": "Megaptera novaeangliae", "minke_whale": "Balaenoptera acutorostrata",
    "sperm_whale": "Physeter macrocephalus",   "killer_whale": "Orcinus orca",
    "common_dolphin": "Delphinus delphis",     "bottlenose_dolphin": "Tursiops truncatus",
}
# Model species that stand for a whole genus or family: OBIS is also queried at that level, and
# obsval.py validates them against it (saved as grid3_<taxon>.geojson)
TAXA = {"hake": "Merluccius", "flounder": "Pleuronectidae", "skate": "Rajidae",
        "lanternfish": "Myctophidae", "bristlemouth": "Cyclothone", "sardinella": "Sardinella", "anchovy": "Engraulis"}
ICCAT = ["Data/t1nc_20260129.zip", "Data/Catdis/cdis5024_bySpecies.zip", "Data/Catdis/cdis5024_all.zip",
         "Data/EFFDIS_LL2000-2024.csv", "Data/Catalogues/SCRScatalogues_MainSpecies_1995-2024.xlsx"]
RAM = "https://zenodo.org/api/records/11995054"
SAG = "https://sag.ices.dk/SAG_API/api"
DATRAS = "https://datras.ices.dk/WebServices/DATRASWebService.asmx"
MATCH = ("cod", "herring", "capelin", "hake", "sole", "plaice", "halibut", "sardine",
         "anchovy", "sprat", "mackerel", "whiting", "haddock", "saithe", "ray", "skate")


def get(url, timeout=180):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
        return r.read()


def save(path, data, force=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not force:
        print(f"    kept {path.relative_to(HERE)} ({path.stat().st_size/1e6:.1f} MB)")
        return False
    path.write_bytes(data)
    print(f"    {path.relative_to(HERE)} ({len(data)/1e6:.2f} MB)")
    return True


def fetch_ram(force):
    print("RAM Legacy Stock Assessment Database v4.65")
    rec = json.loads(get(RAM))
    for f in rec["files"]:
        save(HERE / "ram_legacy" / f["key"], get(f["links"]["self"], timeout=900), force)
    save(HERE / "ram_legacy" / "zenodo_record.json", json.dumps(rec, indent=1).encode(), True)


def fetch_sag(force):
    print("ICES Stock Assessment Graphs (SSB / F / recruitment / catch)")
    d = HERE / "ices_sag"
    stocks = json.loads(get(f"{SAG}/StockList?year=2023"))
    save(d / "stocklist_2023.json", json.dumps(stocks, indent=1).encode(), True)
    hits = [s for s in stocks
            if any(k in (str(s.get("SpeciesName", "")) + str(s.get("StockDescription", ""))).lower()
                   for k in MATCH)]
    print(f"    {len(hits)} of {len(stocks)} stocks match fishNET species")
    ok = 0
    for s in hits:
        out = d / "summary" / f"{s['StockKeyLabel']}_{s['AssessmentKey']}.json"
        if out.exists() and not force:
            ok += 1
            continue
        try:
            body = get(f"{SAG}/SummaryTable?assessmentKey={s['AssessmentKey']}", timeout=90)
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_bytes(body)
            ok += 1
        except Exception as e:
            print(f"    skip {s['StockKeyLabel']}: {e}")
        time.sleep(0.3)
    print(f"    {ok} stock summary tables in ices_sag/summary/")


def fetch_datras(force):
    print("ICES DATRAS bottom-trawl survey index")
    d = HERE / "ices_datras"
    save(d / "survey_list.xml", get(f"{DATRAS}/getSurveyList"), force)
    for sv in ("NS-IBTS", "BITS", "EVHOE", "SWC-IBTS", "IE-IGFS"):
        try:
            save(d / f"years_{sv}.xml", get(f"{DATRAS}/getSurveyYearList?survey={sv}", timeout=90), force)
        except Exception as e:
            print(f"    skip {sv}: {e}")


def fetch_obis(force):
    print("OBIS occurrence statistics and gridded distributions")
    d = HERE / "obis"
    for key, sci in SPECIES.items():
        q = urllib.parse.quote(sci)
        try:
            save(d / f"stats_{key}.json", get(f"https://api.obis.org/v3/statistics?scientificname={q}"), force)
            save(d / f"grid3_{key}.geojson",
                 get(f"https://api.obis.org/v3/occurrence/grid/3?scientificname={q}", timeout=300), force)
        except Exception as e:
            print(f"    skip {key}: {e}")
        time.sleep(0.4)
    for key, taxon in TAXA.items():
        q = urllib.parse.quote(taxon)
        try:
            save(d / f"stats_{taxon}.json", get(f"https://api.obis.org/v3/statistics?scientificname={q}"), force)
            save(d / f"grid3_{taxon}.geojson",
                 get(f"https://api.obis.org/v3/occurrence/grid/3?scientificname={q}", timeout=300), force)
        except Exception as e:
            print(f"    skip {taxon}: {e}")
        time.sleep(0.4)


def fetch_gbif(force):
    print("GBIF occurrence counts and year histograms")
    d = HERE / "gbif"
    for key, sci in SPECIES.items():
        q = urllib.parse.quote(sci)
        try:
            save(d / f"search_{key}.json",
                 get("https://api.gbif.org/v1/occurrence/search"
                     f"?scientificName={q}&limit=0&facet=year&facetLimit=200"), force)
        except Exception as e:
            print(f"    skip {key}: {e}")
        time.sleep(0.3)


def fetch_iccat(force):
    print("ICCAT Atlantic tuna and swordfish catch")
    for rel in ICCAT:
        try:
            save(HERE / "iccat" / Path(rel).name, get(f"https://www.iccat.int/{rel}", timeout=600), force)
        except Exception as e:
            print(f"    skip {rel}: {e}")


JOBS = dict(ram=fetch_ram, sag=fetch_sag, datras=fetch_datras,
            obis=fetch_obis, gbif=fetch_gbif, iccat=fetch_iccat)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--force", action="store_true", help="re-download files that already exist")
    p.add_argument("--skip", default="", help="comma-separated: " + ",".join(JOBS))
    a = p.parse_args()
    skip = {x.strip() for x in a.skip.split(",") if x.strip()}
    for name, fn in JOBS.items():
        if name in skip:
            continue
        try:
            fn(a.force)
        except Exception as e:
            print(f"  FAILED {name}: {e}", file=sys.stderr)
    print("\ndone")
