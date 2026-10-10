#!/usr/bin/env python3
"""Aggregate a fishNET run into FEISTY's functional types, for comparison with FEISTY.

    python feisty_compare.py runs/<name>                       # setupBasic types, last 365 days
    python feisty_compare.py runs/<name> --setup vertical      # setupVertical types
    python feisty_compare.py runs/<name> --window 0 --map squid=large_pelagic

FEISTY (Petrik et al. 2019, Prog. Oceanogr. 176:102124; van Denderen et al. 2021, Glob. Ecol. Biogeogr.
30:1822; R/Fortran code at github.com/Kenhasteandersen/FEISTY) has no species. Its fish are size-structured
functional types eating two zooplankton pools and an unstructured benthic invertebrate pool:

  setupBasic     smallPel (forage fish, 0.001-250 g), largePel (large pelagics, 0.001-125,000 g) and demersals
                 (0.001-125,000 g), plus benthos
  setupVertical  adds mesoPel (mesopelagic fish, 0.001-250 g) and midwPred (midwater predators, 0.001-125,000 g)

Size classes are FEISTY's three in setupBasic: small 0.001-0.5 g (larvae), medium 0.5-250 g (juveniles, adult
forage fish) and large 250-125,000 g (adult large fish). Benthic production is 0.1 x the detrital flux reaching
the bottom (setupBasic: bprod = dfbot * 0.1).

How fishNET maps onto that:
  types    every species file carries `feisty = "forage" | "mesopelagic" | "large_pelagic" | "midwater_predator" |
           "demersal" | "none"` (none = not a fish FEISTY represents: cephalopods and marine mammals). In setupBasic
           the mesopelagic and midwater types are kept apart as "outside setupBasic" rather than folded in, since
           myctophids alone can outweigh every other type. --map name=type overrides a species.
  sizes    each species/stage/cell's mean individual mass (fish_biomass / fish_numbers in fields.nc, agents
           included) puts that class in a size class. Eggs and anything under 0.001 g are left out, as FEISTY's
           spectrum starts there; forage and mesopelagic biomass above 250 g counts as medium, and fish above
           125 kg count as large.
  SSB      adult-stage biomass (FEISTY's SSB is mature biomass).
  yield    catch from series.nc, by type.
  benthos  with npzd.detritus_food = "benthos", fishNET's benthic fauna pool (fields.nc `benthos`), the direct
           counterpart of FEISTY's; for older runs the proxy, the detritus in each column's bed cell, as g wet mass
           m-2 (1.8 mmol N per g). Benthic production is FEISTY's 0.1 x the sinking detritus flux onto the bed.
  forcing  the inputs FEISTY needs at each column, for running it against the same ocean: pelagic temperature
           (Tp, mean over the top 100 m), bottom temperature (Tb), bottom depth, detrital flux to the bed (dfbot),
           and zooplankton standing stocks (small = mesozooplankton, large = krill; biomass, not production).

Output in runs/<name>/feisty/: feisty.nc (series and time-mean maps), feisty_summary.csv, feisty.png.
Needs numpy and netCDF4 (matplotlib for the figure).
"""
import argparse, ast, csv
from pathlib import Path
import numpy as np
import netCDF4 as nc

RHO_N = 1.8                                     # mmol N per g wet mass (fishnet.py)
SIZE_EDGES = np.array([0.001, 0.5, 250.0])      # g: small / medium / large lower bounds (FEISTY setupBasic)
SIZES = ["small", "medium", "large"]
BENTHIC_EFF = 0.1                               # FEISTY setupBasic: bprod = dfbot * 0.1
SETUPS = {                                      # fishNET type -> FEISTY group, per setup
    "basic": {"forage": "smallPel", "large_pelagic": "largePel", "demersal": "demersals",
              "mesopelagic": "outside_setupBasic", "midwater_predator": "outside_setupBasic"},
    "vertical": {"forage": "smallPel", "mesopelagic": "mesoPel", "large_pelagic": "largePel",
                 "midwater_predator": "midwPred", "demersal": "demersals"},
}
SMALL_TYPES = {"forage", "mesopelagic"}          # FEISTY's small fish types stop at 250 g
# For runs whose species files predate the `feisty` key
DEFAULT = {**{n: "forage" for n in "anchoveta sardine_eu sardine_pacific sardinella anchovy herring capelin sardine "
                                   "polar_cod silverfish warm_sardine".split()},
           **{n: "mesopelagic" for n in ("lanternfish", "bristlemouth")},
           **{n: "large_pelagic" for n in "skipjack yellowfin albacore bluefin_atlantic bluefin_pacific swordfish tuna "
                                         "warm_tuna shark blue_shark shortfin_mako".split()},
           **{n: "demersal" for n in "cod_atlantic cod pollock hake halibut greenland_halibut flounder flatfish skate "
                                    "toothfish spiny_dogfish".split()}}
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7"]   # fixed order, never cycled


def species_types(cfg, names, overrides):
    spec = {s.get("name"): s for s in cfg.get("species", [])}
    out = []
    for n in names:
        t = overrides.get(n, spec.get(n, {}).get("feisty", DEFAULT.get(n, "none")))
        out.append(t)
    return out


def bed_cell(cfg, depth):
    """Index and water thickness of each column's deepest wet cell (partial bottom cells, as in fishnet.Grid)."""
    ze = np.asarray(cfg["grid"]["depth_edges"], float)
    d = np.clip(np.nan_to_num(depth), 0, ze[-1])
    frac = np.clip((d[None] - ze[:-1, None, None]) / np.diff(ze)[:, None, None], 0, 1)
    wet = frac > 1e-3
    k = np.maximum(wet.sum(0) - 1, 0)
    dz = np.take_along_axis(np.diff(ze)[:, None, None] * frac, k[None], 0)[0]
    return k, dz, ze


def main():
    ap = argparse.ArgumentParser(description="aggregate a fishNET run into FEISTY functional types")
    ap.add_argument("run")
    ap.add_argument("--setup", choices=sorted(SETUPS), default="basic")
    ap.add_argument("--window", type=float, default=365, help="days at the end of the run to average (0 = all)")
    ap.add_argument("--map", nargs="*", default=[], metavar="SPECIES=TYPE",
                    help="override a species' type (forage, mesopelagic, large_pelagic, midwater_predator, demersal, none)")
    ap.add_argument("--out", help="output folder (default runs/<name>/feisty)")
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()
    run = Path(a.run)
    out = Path(a.out) if a.out else run / "feisty"
    out.mkdir(parents=True, exist_ok=True)
    over = dict(x.split("=", 1) for x in a.map)

    F, S = nc.Dataset(run / "fields.nc"), nc.Dataset(run / "series.nc")
    for ds in (F, S):
        ds.set_auto_mask(False)
    cfg = ast.literal_eval(F.getncattr("config"))
    names = F.getncattr("species").split(",")
    lon, lat = F["lon"][:], F["lat"][:]
    mask = F["mask"][:] > 0.5
    ze_lat = np.linspace(*cfg["grid"]["lat"], len(lat) + 1)
    area = (6.371e6 ** 2 * np.radians(abs(lon[1] - lon[0])) * np.abs(np.diff(np.sin(np.radians(ze_lat)))))[:, None] \
        * np.ones(len(lon))[None] * mask
    A_ocean = area.sum()

    ftype = species_types(cfg, names, over)
    groups = list(dict.fromkeys(SETUPS[a.setup].values())) + ["not_in_FEISTY"]
    gi = np.array([groups.index(SETUPS[a.setup].get(t, "not_in_FEISTY")) for t in ftype])
    small = np.array([t in SMALL_TYPES for t in ftype])
    ng, ns = len(groups), len(SIZES)

    tf = F["time"][:]
    t1 = tf.max()
    use = np.nonzero(tf >= (t1 - a.window if a.window > 0 else -np.inf))[0]
    k_bed, dz_bed, ze = bed_cell(cfg, F["bottom_depth"][:])
    zc = 0.5 * (ze[1:] + ze[:-1])
    top = zc <= 100.0
    w_sink = float(cfg.get("npzd", {}).get("w_sink", 10.0))
    live_benthos = "benthos" in F.variables and str(cfg.get("npzd", {}).get("detritus_food", "")).lower() == "benthos"
    benthos_label = "benthic fauna pool" if live_benthos else "bed-cell detritus proxy"
    dze = np.diff(ze)[:, None, None]
    thick = dze * np.clip((np.nan_to_num(F["bottom_depth"][:])[None] - ze[:-1, None, None]) / dze, 0, 1)   # water per cell

    # ---- per-frame aggregation from fields.nc
    nt = len(tf)
    Bts = np.zeros((nt, ng, ns))                 # g, by group and size class
    SSBts = np.zeros((nt, ng))
    benth_ts = np.zeros(nt)
    Bmap = np.zeros((ng, ns, len(lat), len(lon)))
    SSBmap = np.zeros((ng, len(lat), len(lon)))
    env = {k: np.zeros((len(lat), len(lon))) for k in ("Tp", "Tb", "dfbot", "benthos", "bprod", "zoo_small", "zoo_large")}
    for n in range(nt):
        B = np.nan_to_num(F["fish_biomass"][n])           # (species, stage, lat, lon) g m-2
        N = np.nan_to_num(F["fish_numbers"][n])
        W = np.where(N > 0, B / np.where(N > 0, N, 1), 0)
        cls = np.searchsorted(SIZE_EDGES, W, side="right") - 1          # -1 = below 0.001 g (eggs)
        cls = np.where(small[:, None, None, None], np.minimum(cls, 1), cls)
        keep = (cls >= 0) & (np.arange(B.shape[1]) > 0)[None, :, None, None]   # FEISTY has no egg stage
        Bg = np.zeros((ng, ns, len(lat), len(lon)))
        idx = np.broadcast_to(gi[:, None, None, None], B.shape)
        np.add.at(Bg, (idx[keep], cls[keep], *np.nonzero(keep)[2:]), B[keep])
        ssb = np.zeros((ng, len(lat), len(lon)))
        np.add.at(ssb, gi, B[:, -1])
        D = np.nan_to_num(F["D"][n])
        Dbed = np.take_along_axis(D, k_bed[None], 0)[0]
        benth = (np.nan_to_num(F["benthos"][n]) if live_benthos else Dbed * dz_bed / RHO_N) * mask   # g wet m-2
        Bts[n] = (Bg * area).sum((-2, -1))
        SSBts[n] = (ssb * area).sum((-2, -1))
        benth_ts[n] = (benth * area).sum()
        if n in use:
            Bmap += Bg / len(use)
            SSBmap += ssb / len(use)
            T = np.nan_to_num(F["temp"][n])
            env["Tp"] += T[top].mean(0) / len(use)
            env["Tb"] += np.take_along_axis(T, k_bed[None], 0)[0] / len(use)
            df = Dbed * w_sink * 365 / RHO_N * mask                      # g wet m-2 yr-1 reaching the bed
            env["dfbot"] += df / len(use)
            env["bprod"] += BENTHIC_EFF * df / len(use)
            env["benthos"] += benth / len(use)
            for key, var in (("zoo_small", "Z"), ("zoo_large", "krill")):
                env[key] += (np.nan_to_num(F[var][n]) * thick).sum(0) / RHO_N * mask / len(use)

    # ---- yield (catch) by group, from series.nc
    ts = S["time"][:]
    catch = S["catch"][:]                                                # g d-1 (time, species)
    Y = np.zeros((len(ts), ng))
    for q in range(len(names)):
        Y[:, gi[q]] += catch[:, q]
    ys = ts >= (ts.max() - a.window if a.window > 0 else -np.inf)
    Ymean = Y[ys].mean(0) * 365                                          # g yr-1

    # ---- summary
    uw = np.isin(np.arange(nt), use)
    Bm, SSBm, benth_m = Bts[uw].mean(0), SSBts[uw].mean(0), benth_ts[uw].mean()
    rows = []
    for g, name in enumerate(groups):
        members = [n for n, x in zip(names, gi) if x == g]
        rows.append(dict(group=name, biomass_Mt=Bm[g].sum() / 1e12, biomass_g_m2=Bm[g].sum() / A_ocean,
                         **{f"{s}_g_m2": Bm[g, k] / A_ocean for k, s in enumerate(SIZES)},
                         SSB_g_m2=SSBm[g] / A_ocean, yield_kt_yr=Ymean[g] / 1e9, yield_g_m2_yr=Ymean[g] / A_ocean,
                         species=" ".join(members)))
    rows.append(dict(group="benthos", biomass_Mt=benth_m / 1e12, biomass_g_m2=benth_m / A_ocean,
                     species=f"{benthos_label}; production {float((env['bprod'] * area).sum() / A_ocean):.3g} "
                             f"g m-2 yr-1 = {BENTHIC_EFF} x detrital flux to the bed"))
    keys = ["group", "biomass_Mt", "biomass_g_m2", *[f"{s}_g_m2" for s in SIZES], "SSB_g_m2", "yield_kt_yr",
            "yield_g_m2_yr", "species"]
    with open(out / "feisty_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        w.writerows(rows)
    print(f"{run.name}: FEISTY {a.setup} types, model days {tf[use].min():.0f}-{tf[use].max():.0f} ({len(use)} frames), "
          f"ocean area {A_ocean:.3g} m2")
    print(f"{'group':20s} {'B Mt':>9s} {'B g/m2':>8s} {'S':>7s} {'M':>7s} {'L':>7s} {'SSB':>7s} {'yield kt/yr':>12s}")
    for r in rows:
        print(f"{r['group']:20s} {r['biomass_Mt']:9.3g} {r['biomass_g_m2']:8.3g} "
              + " ".join(f"{r.get(f'{s}_g_m2', np.nan):7.3g}" for s in SIZES)
              + f" {r.get('SSB_g_m2', np.nan):7.3g} {r.get('yield_kt_yr', np.nan):12.4g}")
    print("  " + "\n  ".join(f"{r['group']}: {r['species']}" for r in rows if r["species"]))

    # ---- netCDF
    with nc.Dataset(out / "feisty.nc", "w") as o:
        o.setncatts(dict(source=str(run), setup=a.setup, window_days=a.window, groups=",".join(groups),
                         sizes=",".join(SIZES), size_edges_g="0.001,0.5,250,125000",
                         mapping=";".join(f"{n}={t}" for n, t in zip(names, ftype)),
                         note="biomass in g wet mass; maps are time means over the window"))
        for d, v in (("time", None), ("group", ng), ("size", ns), ("lat", len(lat)), ("lon", len(lon))):
            o.createDimension(d, v)
        o.createVariable("time", "f8", ("time",))[:] = tf
        o["time"].units = F["time"].units
        o.createVariable("lat", "f8", ("lat",))[:] = lat
        o.createVariable("lon", "f8", ("lon",))[:] = lon

        def put(k, dims, x, units):
            v = o.createVariable(k, "f4", dims, zlib=True)
            v.units = units
            v[:] = x
        nanmask = lambda x: np.where(mask, x, np.nan)
        put("biomass", ("time", "group", "size"), Bts / A_ocean, "g m-2 (ocean mean)")
        put("ssb", ("time", "group"), SSBts / A_ocean, "g m-2 (ocean mean)")
        put("benthos", ("time",), benth_ts / A_ocean, "g m-2 (ocean mean)")
        put("biomass_map", ("group", "size", "lat", "lon"), nanmask(Bmap), "g m-2")
        put("ssb_map", ("group", "lat", "lon"), nanmask(SSBmap), "g m-2")
        put("yield", ("group",), Ymean / A_ocean, "g m-2 yr-1 (ocean mean)")
        for k, u in (("Tp", "degC, mean of levels above 100 m"), ("Tb", "degC, bed cell"), ("dfbot", "g m-2 yr-1"),
                     ("benthos", "g m-2"), ("bprod", "g m-2 yr-1"), ("zoo_small", "g m-2"), ("zoo_large", "g m-2")):
            put(k + "_map" if k == "benthos" else k, ("lat", "lon"), nanmask(env[k]), u)
        put("bottom_depth", ("lat", "lon"), nanmask(F["bottom_depth"][:]), "m")
    if not a.no_plot:
        plot(out, run.name, groups, tf, Bts / A_ocean, benth_ts / A_ocean, Bmap, lon, lat, mask, benthos_label)
    print(f"wrote {out / 'feisty.nc'}, {out / 'feisty_summary.csv'}" + ("" if a.no_plot else f", {out / 'feisty.png'}"))


def plot(out, name, groups, t, B, benth, Bmap, lon, lat, mask, benthos_label="benthos"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    show = [g for g in range(len(groups)) if B[:, g].sum() > 0]
    nmap = len(show)
    fig = plt.figure(figsize=(4.2 * max(nmap, 2), 8.2))
    gs = fig.add_gridspec(2, max(nmap, 2), height_ratios=[1, 1.1])
    ax = fig.add_subplot(gs[0, :])
    for c, g in enumerate(show):
        ax.semilogy(t, B[:, g].sum(1), lw=2, color=COLORS[c % len(COLORS)], label=groups[g])
    ax.semilogy(t, benth, lw=2, color="0.45", ls="--", label=f"benthos ({benthos_label})")
    ax.set(xlabel="model day", ylabel="g m$^{-2}$ (ocean mean)", title=f"{name}: biomass by FEISTY group")
    ax.legend(frameon=False, ncol=min(len(show) + 1, 4), fontsize=8)
    for c, g in enumerate(show):
        m = fig.add_subplot(gs[1, c])
        x = np.where(mask, Bmap[g].sum(0), np.nan)
        cm = plt.get_cmap("Blues").copy()
        cm.set_bad("0.85")
        im = m.pcolormesh(lon, lat, np.log10(np.maximum(x, 1e-6)), cmap=cm, shading="auto")
        m.set_title(f"{groups[g]}: log10 g m$^{{-2}}$", fontsize=9)
        m.set_aspect("equal")
        m.grid(False)
        fig.colorbar(im, ax=m, shrink=0.7)
    fig.tight_layout()
    fig.savefig(out / "feisty.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
