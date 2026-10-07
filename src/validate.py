#!/usr/bin/env python3
"""fishNET validation suite: small runs that check conservation, positivity, transport, genetics,
aggregation/disaggregation, determinism, restarts, the netCDF forcing path and IBM vs biomass consistency.

    python validate.py [namelist.toml] [--long]
"""
import copy, sys, tempfile, time
from pathlib import Path
import numpy as np
import netCDF4 as nc
import xarray as xr
import fishnet as fn

RESULTS = []


def check(name, ok, detail):
    RESULTS.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:34s} {detail}", flush=True)


def small(cfg, days=20, **over):
    c = copy.deepcopy(cfg)
    c["run"].update(days=days, dt_hours=12, verbose=0, spinup_days=0)
    c["grid"].update(nx=16, ny=12)
    c["ocean"].update(sw_spinup_days=90.0, sw_cache_dir=str(Path(tempfile.gettempdir()) / "fishnet_sw_cache"))
    c["hybrid"].update(mode="dynamic", dynamic_fraction=0.25, update_days=5, agents_per_class=4, max_agents=20000)
    for key, v in over.items():
        sect, k = key.split("__")
        c[sect][k] = v
    return c


def closed(cfg):
    return small(cfg, npzd__restore_days=0, npzd__upwelling_days=0, npzd__export_bottom=False)


def test_budget(cfg):
    m = fn.Model(closed(cfg)).run()
    err = (m.total_N() - m.total0) / m.total0
    check("N conservation (closed domain)", abs(err) < 1e-10 and m.ext == 0, f"rel. drift {err:+.2e} after {m.t:.0f} d")
    m = fn.Model(small(cfg, days=30)).run()
    err = (m.total_N() - m.total0 - m.ext) / m.total0
    check("N budget closure (open domain)", abs(err) < 1e-10, f"rel. error {err:+.2e}, external {m.ext / m.total0:+.2e}")
    a = m.ag
    neg = min(m.npzd.C.min(), m.EN.min(), m.EB.min(), m.EG.min(), a.n.min(initial=0), a.W.min(initial=1), a.E.min(initial=0))
    check("positivity of all state", neg >= -1e-12, f"min value {neg:.2e}")


def test_transport(cfg):
    for src, shelf in (("synthetic", 0), ("synthetic", 3), ("shallow_water", 3)):
        m = fn.Model(small(cfg, days=1, ocean__source=src, ocean__shelf_cells=shelf))
        f = m.ocean(0.0, m.dt)
        M = m.g.vol.copy()
        for _ in range(200):
            M = m.transport3d(M, f["Qx"], f["Qy"], 2000.0)
        C = M[m.g.wet] / m.g.vol[m.g.wet]
        check(f"uniform tracer preserved ({src}, shelf {shelf})", np.abs(C - 1).max() < 1e-12,
              f"max |C-1| = {np.abs(C - 1).max():.1e} (200 steps), max face speed "
              f"{np.abs(f['U'] / m.g.dy).max():.2f} m/s")
        X = m.rng.random(M.shape) * m.g.wet
        Y = m.transport3d(X, f["Qx"], f["Qy"], 2000.0)
        check(f"transport conserves mass ({src}, shelf {shelf})", abs(Y.sum() / X.sum() - 1) < 1e-13, f"rel. change {Y.sum() / X.sum() - 1:+.1e}")


def test_shallow_water(cfg):
    c = small(cfg, days=1, ocean__source="shallow_water")
    g = fn.Grid(c["grid"])
    fn.Ocean(dict(c["ocean"], source="synthetic"), g, fn.dtm.datetime(2000, 1, 1), 0)
    sw = fn.ShallowWater(c["ocean"], g, 0.0, 0.0, 0.0)
    sw.advance(365.0)
    ke = 0.5 * (sw.u[0] ** 2 + sw.v[0] ** 2)[1:-1, 1:-1]
    check("shallow water stable for a year", np.isfinite(sw.h).all() and np.sqrt(2 * ke.max()) < 3,
          f"max speed {np.sqrt(2 * ke.max()):.2f} m/s, h1 {sw.h[0, 1:-1, 1:-1].min():.0f}-{sw.h[0, 1:-1, 1:-1].max():.0f} m")
    m = fn.Model(small(cfg, days=30, ocean__source="shallow_water", npzd__restore_days=0, npzd__upwelling_days=0,
                       npzd__export_bottom=False))
    m.run()
    err = abs(m.diagnostics()["budget_error"])
    check("N conserved under shallow-water forcing", err < 1e-12, f"relative error {err:.1e} after 30 d")


def test_cmems(cfg):
    """End to end on a GLORYS-lookalike file: CMEMS names, descending latitude, packed ints, 1950 epoch."""
    import subprocess
    d = Path(tempfile.mkdtemp())
    r = subprocess.run([sys.executable, str(Path(__file__).parent / "make_test_forcing.py"), "--out", str(d),
                        "--res", "1.0", "--months", "4", "--levels", "10"], capture_output=True, text=True)
    if r.returncode:
        check("GLORYS-lookalike forcing", False, r.stderr.strip().splitlines()[-1][:90])
        return
    c = small(cfg, days=10)
    c["grid"].update(nx=18, ny=12, lon=[-66.0, -16.0], lat=[24.0, 56.0], depth_edges=[0, 30, 100, 300, 800])
    c["ocean"] = dict(source="netcdf", file=str(d / "forcing.nc"), bathymetry=str(d / "bathymetry.nc"), T_deep=4.0)
    c["run"]["start"] = "2010-01-15T00:00"
    m = fn.Model(c)
    g = m.g
    land, T = (~g.mask).sum(), m.f["T"][0][g.mask]
    ok = 0 < land < g.mask.size / 2 and 0 < T.min() and T.max() < 40 and (g.depth[g.mask] > 0).all()
    check("CMEMS-style file read correctly", ok, f"{g.mask.sum()} ocean columns, {land} land, SST "
          f"{T.min():.1f}-{T.max():.1f} degC, depth {g.depth[g.mask].min():.0f}-{g.depth[g.mask].max():.0f} m")
    m.run()
    e = m.diagnostics()["budget_error"]
    check("N conserved on reanalysis-style forcing", abs(e) < 1e-12, f"relative error {e:+.1e} after 10 d")


def test_skill(cfg):
    import skill
    rng = np.random.default_rng(1)
    o = np.sin(np.linspace(0, 6, 200)) + 5
    mt = skill.metrics(o, o)
    perfect = abs(mt["bias"]) < 1e-12 and mt["rmse"] < 1e-12 and abs(mt["corr"] - 1) < 1e-9 and abs(mt["nse"] - 1) < 1e-9
    shifted = skill.metrics(o + 0.5, o)
    scaled = skill.metrics(2 * o, o)
    ok = perfect and abs(shifted["bias"] - 0.5) < 1e-9 and abs(shifted["crms"]) < 1e-9 \
        and abs(scaled["sd_ratio"] - 2) < 1e-9
    check("skill metrics", ok, f"identical series scores bias 0, RMSE 0, corr 1, NSE 1; a +0.5 offset gives "
          f"bias {shifted['bias']:.2f} with zero centred RMS; doubling gives sd ratio {scaled['sd_ratio']:.1f}")
    ens = rng.normal(0, 1, 4000)
    c_at, c_off = skill.crps(ens, 0.0), skill.crps(ens, 3.0)
    ideal = 2 / np.sqrt(2 * np.pi) - 1 / np.sqrt(np.pi)                   # analytic CRPS of N(0,1) at its mean
    check("ensemble CRPS and rank", abs(c_at - ideal) < 0.02 and c_off > 2 and
          abs(skill.rank_histogram(ens, 0.0) / len(ens) - 0.5) < 0.05,
          f"CRPS {c_at:.3f} at the mean (ideal {ideal:.3f}), {c_off:.2f} three sigma away; median observation "
          "ranks mid-ensemble")
    t = np.arange(0, 1460, 5.0)
    sst = 15 + 3 * np.sin(2 * np.pi * t / 365)
    field = sst[:, None, None] + np.zeros((1, 4, 4))
    field[(t > 700) & (t < 760)] += 3.0                                   # a 60-day basin-wide warm event
    F = xr.Dataset({"temp": (("time", "depth", "lat", "lon"), field[:, None])}, coords={"time": t})
    hw = skill.find_heatwaves(F, threshold=1.5)
    got = [h for h in hw if abs(h["start"] - 700) < 30]
    check("heatwave detection", len(got) == 1 and got[0]["peak"] > 2.5,
          f"{len(hw)} event(s) found; the imposed one starts day {got[0]['start']:.0f} with peak "
          f"+{got[0]['peak']:.1f} degC" if got else "imposed event missed")
    b = np.r_[np.linspace(1, 10, 50), np.linspace(10, 2, 50)]
    ds = xr.Dataset({"biomass": (("time", "species", "stage"), b[:, None, None] * np.ones((1, 1, 1)))},
                    coords={"time": np.arange(100.0)})
    col = skill.find_collapses(ds, ["x"], drop=0.5)
    check("collapse detection", len(col) == 1 and col[0]["ratio"] < 0.5,
          f"one collapse found, to {col[0]['ratio']:.0%} of peak by day {col[0]['time']:.0f}" if col else "missed")


def test_calibration(cfg):
    import calibrate as cal
    rng = np.random.default_rng(3)
    tg = [dict(name="b", kind="biomass", value=100.0, sigma=0.3, log=True, weight=1.0),
          dict(name="c", kind="catch", value=10.0, sigma=0.5, log=True, weight=1.0)]
    exact, off = cal.loss(np.array([100.0, 10.0]), tg), cal.loss(np.array([100.0 * np.e ** 0.3, 10.0]), tg)
    ok = exact < 1e-12 and abs(off - 1) < 1e-9 and not np.isfinite(cal.loss(np.array([-1.0, 10.0]), tg))
    check("calibration loss", ok, f"0 at the observation, {off:.2f} at one sigma, infinite for impossible values")
    b, lg = np.array([[1e-3, 1e-1], [0.0, 10.0]]), np.array([True, False])
    u = rng.random((200, 2))
    ok = np.allclose(cal.to_unit(cal.from_unit(u, b, lg), b, lg), u) and \
        np.allclose(cal.from_unit(np.array([[0.5, 0.5]]), b, lg)[0], [1e-2, 5.0])
    check("prior transforms round-trip", ok, "log and linear priors map to the unit cube and back")
    d = cal.lhs(50, 3, rng)
    strat = all(((d[:, j] >= k / 50) & (d[:, j] < (k + 1) / 50)).sum() == 1 for j in range(3) for k in range(50))
    check("latin hypercube design", strat, "50 points, one per stratum in every dimension")
    X = rng.random((60, 2))
    truth = lambda x: 3 * (x[:, 0] - 0.3) ** 2 + 2 * (x[:, 1] - 0.7) ** 2
    gp = cal.GP(X, truth(X))
    Xs = rng.random((200, 2))
    pred, var = gp(Xs, var=True)
    err = np.abs(pred - truth(Xs)).max()
    check("emulator reproduces a known surface", err < 0.1 and gp.loo_r2() > 0.9 and (var >= 0).all(),
          f"max error {err:.3f} over 200 points, leave-one-out R2 {gp.loo_r2():.3f}")
    post, acc = cal.mcmc(lambda x: -0.5 * ((x - 0.35) ** 2).sum(-1) / 0.1 ** 2, 2, 40000, rng, step=0.25)
    mu, sd = post.mean(0), post.std(0)
    check("emulator MCMC samples the posterior", np.abs(mu - 0.35).max() < 0.02 and np.abs(sd - 0.1).max() < 0.02,
          f"recovered mean {mu.round(3)} and sd {sd.round(3)} for a Gaussian at 0.35 +/- 0.1, acceptance {acc:.2f}")


def test_multilocus(cfg):
    c = small(cfg, days=25)
    m = fn.Model(c)
    gen = m.gen
    if not gen.L:
        check("multilocus genetics enabled", False, "genetics.loci = 0")
        return
    s0 = np.zeros(4000, int)
    gl = gen.draw(s0, m.rng)
    A = gen.breeding(gl, s0)
    var_ok = np.abs(A.var(0) / m.sp.h2[0] - 1).max() < 0.15
    check("multilocus variance = h2", var_ok, f"Var(A) {A.var(0)[:3].round(3)} for h2 {m.sp.h2[0][:3]} "
          f"over {gen.L} loci")
    G = gen.g_matrix(A / np.sqrt(m.sp.h2[0]))
    corr = lambda i, j: G[i, j] / np.sqrt(G[i, i] * G[j, j])
    want = {tuple(sorted(k.split("-"))): v for k, v in c.get("genetics", {}).get("correlations", {}).items()}
    err = max(abs(corr(fn.TRAITS.index(a_), fn.TRAITS.index(b_)) - v) for (a_, b_), v in want.items()) if want else 1
    check("pleiotropy gives the requested G", err < 0.1,
          "; ".join(f"{a_}-{b_} {corr(fn.TRAITS.index(a_), fn.TRAITS.index(b_)):+.2f} (want {v:+.2f})"
                    for (a_, b_), v in want.items()))
    mum, dad = gen.draw(s0, m.rng), gen.draw(s0, m.rng)
    kid = gen.meiosis(mum, dad, s0, m.rng)
    Am, Ad, Ak = (gen.breeding(x, s0)[:, 0] for x in (mum, dad, kid))
    mid = 0.5 * (Am + Ad)
    slope = np.polyfit(mid, Ak, 1)[0]
    seg = np.var(Ak - mid) / m.sp.h2[0, 0]
    check("meiosis reproduces the infinitesimal model", abs(slope - 1) < 0.1 and abs(seg - 0.5) < 0.1,
          f"offspring on midparent slope {slope:.2f} (want 1), segregation variance {seg:.2f} h2 (want 0.5)")
    g2 = fn.Genetics(dict(loci=gen.L, mutation_rate=1.0, mutation_sd=0.3), m.sp.h2, m.rng)
    v0 = gen.breeding(gl, s0).var(0)[0]
    gm = gl.copy()
    for _ in range(10):
        gm = g2.meiosis(gm, gm[::-1].copy(), s0, m.rng)
    check("mutation adds genetic variance", g2.breeding(gm, s0).var(0)[0] > v0,
          f"variance {v0:.3f} -> {g2.breeding(gm, s0).var(0)[0]:.3f} after 10 generations of mutation")
    m.run()
    d = m.diagnostics()
    a = m.ag
    ok = len(a) and np.isfinite(d["fst"][np.isfinite(d["fst"])]).all() and (d["sel"] != 0).any()
    hot = d["sel"][:, fn.CAUSES.index("thermal"), fn.T_OPT]
    check("adaptation diagnostics", bool(ok), f"F_ST {np.round(d['fst'], 3)}, thermal selection on T_opt "
          f"{np.round(hot, 2)} SD, genetic variance {np.round(d['g_var'][:, 0], 2)}")


def spy_mating(config):
    """Run a model, recording every realised mating: the pair's cells and their distance in T_opt."""
    m = fn.Model(config)
    rec = {"same": 0, "cross": 0, "dz": []}
    orig = fn.Model.mate

    def mate(self, fem, males, col):
        dad = orig(self, fem, males, col)
        ok = dad >= 0
        if ok.any():
            mi = males[dad[ok]]
            same = col[fem[ok]] == col[mi]
            rec["same"] += int(same.sum())
            rec["cross"] += int((~same).sum())
            rec["dz"].append(np.abs(self.ag.A[fem[ok], fn.T_OPT] - self.ag.A[mi, fn.T_OPT]))
        return dad

    fn.Model.mate = mate
    try:
        m.run()
    finally:
        fn.Model.mate = orig
    rec["dz"] = np.concatenate(rec["dz"]) if rec["dz"] else np.array([np.nan])
    return m, rec


def test_mating(cfg):
    c = small(cfg, days=60)
    c.setdefault("genetics", {}).update(mate_search="column", mate_choice_sd=0.0, mate_choice_sample=5)
    m, local = spy_mating(c)
    check("mating is local (mate_search = column)", local["cross"] == 0 and local["same"] > 0,
          f"{local['same']} matings, none across cells")

    a = m.ag                                        # a female with no male in her own column
    if len(a) > 1:
        fem, males = np.array([0]), np.array([1])
        col = np.zeros(len(a), int)
        col[males] = 7
        a.sp[1] = a.sp[0]
        m.cfg["genetics"]["mate_search"] = "column"
        alone = int(m.mate(fem, males, col)[0])
        m.cfg["genetics"]["mate_search"] = "global"
        far = int(m.mate(fem, males, col)[0])
        check("isolation: no local male -> no spawning", alone < 0 and far >= 0,
              f"column mode returns {alone} (no mate), global mode {far} (a distant male)")

    own = fn.in_season((m.f["doy"] - a.Z[:, fn.SPAWN]) % 365, m.sp.spawn_doy[a.sp])
    base = fn.in_season(np.full(len(a), m.f["doy"]), m.sp.spawn_doy[a.sp])
    shifted = int((own != base).sum())
    check("spawn_shift moves an individual's season", shifted > 0 and a.Z[:, fn.SPAWN].std() > 0,
          f"{shifted}/{len(a)} agents in or out of season on their own timing, "
          f"shifts {a.Z[:, fn.SPAWN].min():+.0f}..{a.Z[:, fn.SPAWN].max():+.0f} d")

    ca = copy.deepcopy(c)
    ca["genetics"].update(mate_choice_sd=1.0)
    _, assort = spy_mating(ca)
    ratio = np.nanmean(assort["dz"]) / np.nanmean(local["dz"])
    check("assortative mating pairs like with like", ratio < 0.95,
          f"mean |dA(T_opt)| between mates {np.nanmean(local['dz']):.2f} -> {np.nanmean(assort['dz']):.2f} "
          f"({ratio:.2f}x) over {assort['same']} matings")


def test_diet_names(cfg):
    c = small(cfg, days=1)
    for s in c["species"]:
        s["fish_diet"] = {k.replace(" ", "_").upper(): v for k, v in s.get("fish_diet", {}).items()}
    loose = fn.Model(c).sp.PREF
    for s, orig in zip(c["species"], small(cfg, days=1)["species"]):
        s["fish_diet"] = orig.get("fish_diet", {})
    exact = fn.Model(c).sp.PREF
    check("diet species names match loosely", np.array_equal(loose, exact) and (exact[:, fn.ADULT, 3:] > 0).any(),
          "WARM_SARDINE / 'warm sardine' / warm_sardine resolve to the same species")


def test_oxygen(cfg):
    T = np.array([0.0, 10.0, 20.0, 30.0])
    sat, ref = fn.o2_saturation(T), np.array([357.0, 282.0, 231.0, 195.0])      # published solubility, 35 psu
    check("oxygen solubility", np.abs(sat / ref - 1).max() < 0.02,
          "mmol m-3 at 0/10/20/30 degC: " + " ".join(f"{x:.0f}" for x in sat) + " (published "
          + " ".join(f"{x:.0f}" for x in ref) + ")")
    m = fn.Model(small(cfg, days=5))
    p, g = m.npzd, m.g
    p.C[fn.O2_IDX] = 0.5 * fn.o2_saturation(m.f["T"])                                   # start the surface undersaturated
    before = p.C[fn.O2_IDX, 0][g.mask].mean()
    flux = p.air_sea(m.f["T"], 7.0, 1.0)
    after = p.C[fn.O2_IDX, 0][g.mask].mean()
    rate = (after - before) * g.dzx[0][g.mask].mean()
    check("air-sea oxygen flux", flux > 0 and before < after < fn.o2_saturation(m.f["T"][0])[g.mask].mean()
          and 20 < rate < 2000, f"undersaturated surface gains {rate:.0f} mmol m-2 d-1, toward saturation")
    m.run()
    d = m.diagnostics()
    C = p.C[fn.O2_IDX][g.wet]
    ok = (C >= 0).all() and C.max() < 1.5 * fn.o2_saturation(0.0) and abs(d["budget_error"]) < 1e-12
    check("oxygen stays physical", ok, f"range {C.min():.0f}-{C.max():.0f} mmol m-3, N budget {d['budget_error']:+.1e}")


def test_metabolic_index(cfg):
    W = np.array([10.0, 1000.0])
    phi = lambda O, T, w: fn.metabolic_index(O, T, w, 0.45, 0.4, -0.1)
    mono = phi(280.0, 15.0, W[0]) > phi(60.0, 15.0, W[0]) and phi(280.0, 15.0, W[0]) > phi(280.0, 15.0, W[1]) \
        and phi(280.0, 25.0, W[0]) < phi(280.0, 15.0, W[0])
    sat = phi(fn.o2_saturation(15.0), 15.0, 100.0)
    check("metabolic index responds correctly", mono and 2 < sat < 12,
          f"phi = {sat:.1f} for a 100 g fish in saturated 15 degC water; falls with low O2, big size and warmth")
    c = small(cfg, days=10, npzd__O2_utilised=0.97, npzd__restore_days=0, npzd__upwelling_days=0)
    m = fn.Model(c)
    m.run()
    d = m.diagnostics()
    hyp = d["losses"][:, fn.CAUSES.index("hypoxia")].sum()
    deep = m.npzd.C[fn.O2_IDX][m.g.wet].min()
    check("hypoxia kills fish in an oxygen minimum", hyp > 0 and deep < 60 and abs(d["budget_error"]) < 1e-12,
          f"min O2 {deep:.0f} mmol m-3, hypoxic loss {hyp / 1e6:.2g} t/d, mean index "
          + " ".join(f"{x:.1f}" for x in d["metabolic_index"]))


def test_export(cfg):
    m = fn.Model(small(cfg, days=20))
    m.run()
    d = m.diagnostics()
    A = (m.g.area * m.g.mask).sum()
    gc = lambda x: x * 12.011 / 1000 / A * 365                                  # mmol C d-1 -> g C m-2 yr-1
    tot = d["export_sinking"] + d["export_fish"] + d["export_respired"]
    ok = d["export_sinking"] > 0 and d["export_fish"] >= 0 and 0.1 < gc(tot) < 1000
    check("carbon export by pathway", ok, f"at {m.npzd.c.get('export_depth', 100):.0f} m: sinking "
          f"{gc(d['export_sinking']):.1f}, fish faeces/carcasses {gc(d['export_fish']):.2f}, fish respiration "
          f"{gc(d['export_respired']):.2f} g C m-2 yr-1")


def test_bathymetry(cfg):
    c = small(cfg, days=20, ocean__shelf_cells=3, ocean__shelf_depth=80.0)
    m = fn.Model(c)
    g = m.g
    vol = (g.dz[:, None, None] * g.area * np.clip((g.depth[None] - g.z_e[:-1, None, None]) / g.dz[:, None, None], 0, 1))
    shelf = (g.depth[g.mask] < g.z_e[-1] - 1).any()
    ok = shelf and np.allclose(g.vol, vol) and (g.vol[~g.wet] == 0).all() and g.wet[0][g.mask].all()
    check("partial bottom cells", ok, f"depths {np.unique(np.round(g.depth[g.mask]))[:4]} m, "
          f"wet cells {g.wet.sum()} of {g.wet.size}")
    dem = np.nonzero(m.sp.demersal)[0]
    pinned = all((m.LEV[d, fn.ADULT] == (g.bed & m.LEV[d, fn.ADULT])).all() for d in dem) if len(dem) else False
    inwater = bool((m.LEV <= g.wet[None, None]).all())
    check("habitat follows the sea bed", pinned and inwater,
          f"{len(dem)} demersal species pinned to the bed, no class in dry cells")
    m.run()
    d = m.diagnostics()
    C = m.npzd.C
    check("N conserved over bathymetry", abs(d["budget_error"]) < 1e-12 and not C[:, ~g.wet].any(),
          f"relative error {d['budget_error']:+.1e}, dry cells empty")


def test_ensemble(cfg):
    c = small(cfg, days=2)
    c["ensemble"] = dict(method="lhs", members=3, processes=2, sweep={"npzd.mu_max": {"min": 1.0, "max": 2.0},
                         "species.*.m0": {"scale": [0.5, 2.0]}, "behavior.engine": ["probabilistic", "softmax"]})
    mem = fn.ensemble_members(c)
    u = np.sort([m["npzd.mu_max"] for m in mem])
    strat = all(1 + k / 3 <= x < 1 + (k + 1) / 3 for k, x in enumerate(u))
    t = copy.deepcopy(c)
    fn.set_param(t, "species.*.m0", 2.0, scale=True)
    scaled = all(np.isclose(a["m0"], 2 * b["m0"]) for a, b in zip(t["species"], c["species"]))
    out = Path(tempfile.mkdtemp()) / "ens"
    fn.ensemble(c, out)
    ds = xr.open_dataset(out / "ensemble.nc", decode_times=False)
    ok = strat and scaled and ds.sizes["member"] == 3 and np.isfinite(ds.biomass.values).all()
    check("ensemble sweep (lhs, scale, parallel)", ok, f"LHS stratified {strat}, per-species scaling {scaled}, "
          f"ensemble.nc biomass {tuple(ds.biomass.shape)}")


def test_genetics():
    rng = np.random.default_rng(3)
    h2, n = 0.4, 400_000
    Am, Af = rng.normal(0, np.sqrt(h2), (2, n))
    zm, zf = Am + rng.normal(0, np.sqrt(1 - h2), n), Af + rng.normal(0, np.sqrt(1 - h2), n)
    Ao = fn.inherit(Am, Af, h2, rng)
    zo = Ao + rng.normal(0, np.sqrt(1 - h2), n)
    slope = np.polyfit(0.5 * (zm + zf), zo, 1)[0]
    check("parent-offspring regression = h2", abs(slope - h2) < 0.01, f"slope {slope:.3f} (h2 = {h2})")
    A = rng.normal(0, np.sqrt(h2), n)
    for _ in range(20):
        A = fn.inherit(A, rng.permutation(A), h2, rng)
    check("genetic variance stable (neutral)", abs(A.var() / h2 - 1) < 0.02, f"var ratio after 20 gens {A.var() / h2:.3f}")
    P = rng.gamma(2.0, 1.0, n) - 2.0                          # skewed parental breeding values
    S = np.stack([(P ** p).sum() for p in range(1, 5)])[None]
    raw = fn.offspring_moments(S, np.array([n]), h2)[0] * n
    mom = [x[0] for x in fn.central_moments(raw[None], np.array([n]))]
    Q = fn.inherit(P, rng.permutation(P), h2, rng)
    mc = [Q.mean(), Q.var(), ((Q - Q.mean()) ** 3).mean(), ((Q - Q.mean()) ** 4).mean()]
    err = max(abs(a - b) / max(abs(b), 0.1) for a, b in zip(mom, mc))
    check("offspring moments (biomass cells)", err < 0.03, f"max rel. error vs Monte Carlo {err:.3f}")
    x = fn.sample_moments(S, np.array([n]), 50_000, rng)[0, :, 0]
    cm = [x[0] for x in fn.central_moments(S, np.array([n]))]
    sk, sk0 = ((x - x.mean()) ** 3).mean() / x.var() ** 1.5, cm[2] / cm[1] ** 1.5
    check("Cornish-Fisher sampling", abs(x.mean() - cm[0]) < 1e-9 and abs(x.var() - cm[1]) < 1e-9 and abs(sk - sk0) < 0.15,
          f"mean/var exact, skew {sk:.2f} vs {sk0:.2f}")


def test_aggregation(cfg):
    m = fn.Model(small(cfg, days=1))
    mass0, n0 = m.fish_mass(), len(m.ag)
    m.merge(np.ones(len(m.ag), bool))
    ok1 = abs(m.fish_mass() / mass0 - 1) < 1e-12
    col = m.ibm & (m.EN.sum((0, 1)) > 0)
    mean0 = m.ES[:, :, :, 0] / np.maximum(m.EN, 1e-300)[:, :, None]
    var0 = fn.central_moments(np.moveaxis(m.ES, 3, -1), m.EN[:, :, None])[1]
    m.disaggregate(col)
    a = m.ag
    j, i, _ = m.cols()
    key = np.ravel_multi_index((a.sp, a.st, j, i), m.EN.shape)
    err = 0.0
    for kk in np.unique(key)[:200]:
        w = key == kk
        s, st, jj, ii = np.unravel_index(kk, m.EN.shape)
        mu = np.average(a.A[w], 0, a.n[w])
        va = np.average((a.A[w] - mu) ** 2, 0, a.n[w])
        err = max(err, np.abs(mu - mean0[s, st, :, jj, ii]).max(), np.abs(va - var0[s, st, :, jj, ii]).max())
    check("aggregate/disaggregate mass", ok1 and abs(m.fish_mass() / mass0 - 1) < 1e-12,
          f"{n0} agents -> fields -> {len(a)} agents, rel. mass change {m.fish_mass() / mass0 - 1:+.1e}")
    check("disaggregated trait mean/var exact", err < 1e-9, f"max |error| {err:.1e} over {min(200, len(np.unique(key)))} classes")


def test_generations(cfg):
    m = fn.Model(small(cfg, days=1))
    for X in (m.EN, m.EB, m.ER, m.EG, m.ES, m.EGEN):          # only the agents' fish in the biomass classes
        X[...] = 0
    m.ag.gen[:] = 3
    n0 = len(m.ag)
    m.merge(np.ones(len(m.ag), bool))
    live = m.EN > 0
    mean = m.EGEN[live] / m.EN[live]
    col = m.ibm & (m.EN.sum((0, 1)) > 0)
    m.disaggregate(col)
    ok = n0 and np.allclose(mean, 3) and len(m.ag) and (m.ag.gen == 3).all() and m.max_generation().max() == 3
    check("generations survive biomass round trip", ok,
          f"{n0} agents of generation 3 -> biomass (mean {mean.mean():.2f}) -> {len(m.ag)} agents of generation "
          f"{sorted(set(m.ag.gen.tolist()))}")


def test_spinup(cfg):
    with tempfile.TemporaryDirectory() as tmp:
        c = small(cfg, days=4, run__spinup_days=6, run__output_every_hours=24)
        m = fn.Model(c, Path(tmp) / "spun").run()
        with nc.Dataset(Path(tmp) / "spun" / "series.nc") as S:        # closing the output retires the agents,
            t, e = np.asarray(S["time"][:]), float(S["budget_error"][-1])  # so read the budget from the file
        plain = fn.Model(small(cfg, days=4)).run()
        ok = t[0] == 0 and abs(t[-1] - 4) < 1e-9 and abs(e) < 1e-10 and m.spin is None \
            and not np.isclose(m.fish_mass(), plain.fish_mass())
        check("spin-up before day 0", ok, f"output days {t[0]:g}..{t[-1]:g} after a 6 d spin-up, N budget {e:+.1e}, "
              f"fish {m.fish_mass() / 1e12:.4g} vs {plain.fish_mass() / 1e12:.4g} Mt without spin-up")


def test_warm_start(cfg):
    with tempfile.TemporaryDirectory() as tmp:
        c = small(cfg, days=6, run__restart_every_days=6, run__restart_keep=1)
        a = fn.Model(c, Path(tmp) / "a")
        a.run()
        st = fn.load_restart(fn.latest_restart(Path(tmp) / "a"))["model"]
        old_mass = (st["EB"].sum() + st["EG"].sum() + (st["ag"].n * (st["ag"].W + st["ag"].E + st["ag"].G)).sum())
        b = fn.Model(small(cfg, days=2, run__warm_start=str(Path(tmp) / "a")))
        mass = b.fish_mass()
        ok = abs(mass / old_mass - 1) < 1e-9 and b.t == 0 and (b.ag.tb <= 0).all() \
            and np.allclose(b.npzd.C, st["npzd"].C) and abs(b.run().diagnostics()["budget_error"]) < 1e-10
        check("warm start from a restart", ok, f"fish {old_mass / 1e12:.4g} Mt carried over (rel. diff {mass / old_mass - 1:+.1e}), "
              f"{len(st['ag'])} agents, clock reset to day 0, N budget closes afterwards")


def test_initial_state(cfg):
    lit = fn.find_file("observations/literature_biomass.toml")
    woa = fn.find_file("observations/woa/woa23_all_o00_01.nc")
    c = small(cfg, days=1)
    names = [s["name"] for s in c["species"]]
    if lit.exists():
        table = __import__("tomllib").loads(lit.read_text())
        gl = fn.find_file("global_species.toml")                 # the table is for the global species list
        c2 = fn.load_config(gl) if gl.exists() else copy.deepcopy(c)
        fn.literature_biomass(c2, verbose=0)
        hit = [s for s in c2["species"] if s["name"] in table]
        ok = len(hit) > 0 and all(np.isclose(s["biomass"], table[s["name"]]["biomass_t"] * 1e6 / 3.6e14) for s in hit)
        tot = sum(table[s["name"]]["biomass_t"] for s in hit)
        check("literature starting biomass", ok, f"{len(hit)} of {len(c2['species'])} species in {gl.name} set from the "
              f"table, {tot / 1e6:,.0f} Mt in all")
    if woa.exists():
        m = fn.Model(small(cfg, days=1, npzd__O2_init="woa"))
        O = m.npzd.C[fn.O2_IDX][m.g.wet]
        ok = np.isfinite(O).all() and O.min() >= 0 and O.max() < 450 and m.npzd.O2_ref is not None
        check("World Ocean Atlas oxygen", ok, f"{O.min():.0f}-{O.max():.0f} mmol m-3 on the test grid, "
              f"{(O < 100).sum()} of {O.size} cells below 100")
    m = fn.Model(small(cfg, days=1))
    S = m.sp
    fade = np.exp(-np.array([1.0, 30.0])[:, None] / S.crowd_L[None])
    check("crowding penalty fades with size", (fade[0] > fade[1]).all() and (S.crowd_L > 0).all(),
          "penalty at 1 cm vs 30 cm: " + ", ".join(f"{n} {a:.2f}/{b:.2f}" for n, a, b in zip(S.names, *fade)))


def test_realism_options(cfg):
    c = small(cfg, days=20, hybrid__juvenile_bins=6, npzd__detritus_food="bed_flux")
    for sp in c["species"]:
        sp.update(plankton_L_max=20.0, K_adult=5.0)
    m = fn.Model(c)
    m.run()
    d = m.diagnostics()
    gap = np.abs(m.NJ.sum(1) - m.EN[:, fn.JUV]).max() / max(m.EN[:, fn.JUV].max(), 1e-30)
    ok = abs(d["budget_error"]) < 1e-10 and gap < 1e-9 and m.NJ.shape[1] == 6 and m.fish_mass() > 0
    check("juvenile bins, bed detritus, adult crowding", ok, f"N budget {d['budget_error']:+.1e}, juvenile bins "
          f"match the class counts to {gap:.0e}, fish {m.fish_mass() / 1e12:.3g} Mt after 20 d")
    tp = m.sp.plankton_taper(np.array([5.0, 20.0, 60.0]), 0)
    check("plankton cut-off by length", tp[0] > 0.9 and abs(tp[1] - 0.5) < 1e-9 and tp[2] < 1e-6,
          f"plankton weight at 5/20/60 cm with plankton_L_max = 20: {tp[0]:.2f}/{tp[1]:.2f}/{tp[2]:.0e}")
    # a one-bin juvenile class leaks fish straight to adulthood; six bins hold them back
    leak = {}
    for K in (1, 6):
        mm = fn.Model(small(cfg, days=1, hybrid__juvenile_bins=K))
        mm.EN[:, fn.ADULT] = mm.EB[:, fn.ADULT] = mm.ER[:, fn.ADULT] = mm.EG[:] = 0
        mm.NJ[:] = 0
        mm.NJ[:, 0] = mm.EN[:, fn.JUV]                      # every juvenile has just arrived
        mm.Tclass = np.full(mm.EN.shape, 10.0)
        g = np.full(mm.EN.shape, 0.0)
        g[:, fn.JUV] = 0.05 * (mm.sp.weight(mm.sp.mu[:, fn.L_MAT]) - mm.sp.weight(mm.sp.L_juv))[:, None, None]
        n0 = mm.EN[:, fn.JUV].sum()
        for _ in range(3):
            mm.promote_euler(g)
        leak[K] = mm.EN[:, fn.ADULT].sum() / n0
    check("no shortcut to maturity with bins", leak[6] < 1e-6 < leak[1],
          f"share of new juveniles adult after 3 steps at 5% of the stage per step: 1 bin {leak[1]:.1%}, "
          f"6 bins {leak[6]:.1e}")


def test_fishing_map(cfg):
    with tempfile.TemporaryDirectory() as tmp:
        table = {"species": {"cod": {"mean": 0.4, "areas": {"27": {"years": {"1993": 0.6, "1995": 0.2}},
                                                              "21": {"years": {"1993": 0.3}}}}}}
        f = Path(tmp) / "F.json"
        f.write_text(__import__("json").dumps(table))
        c = small(cfg, days=1)
        c["fishing"] = dict(mode="ram", file=str(f), scale=0.5, unassessed=0.0)
        m = fn.Model(c)
        q = m.sp.names.index("cod")
        fao = __import__("obsval").fao_area(*np.meshgrid(m.g.lon, m.g.lat)) * m.g.mask
        vals = {}
        for yr in (1993, 1994, 1996):
            c["fishing"]["year"] = yr
            F = m.fishing_rate(np.arange(m.g.ny).repeat(m.g.nx), np.tile(np.arange(m.g.nx), m.g.ny),
                               np.full(m.g.ny * m.g.nx, q)).reshape(m.g.ny, m.g.nx)
            vals[yr] = {a: float(F[fao == a].max()) if (fao == a).any() else None for a in (21, 27, 31)}
        other = [n for n in m.sp.names if n != "cod"]
        none = all(m.fishing_rate(np.array([0]), np.array([0]), np.array([m.sp.names.index(n)]))[0] == 0 for n in other)
        ok = vals[1993][27] == 0.3 and vals[1993][21] in (0.15, None) and vals[1996][27] == 0.1 \
            and vals[1994][27] in (0.3, 0.1) and (vals[1993][31] in (0.0, None)) and none
        check("RAM fishing by area and year", ok, f"cod F x0.5 in area 27: 1993 {vals[1993][27]}, 1996 {vals[1996][27]}; "
              f"area 21 {vals[1993][21]}; unassessed area 31 {vals[1993][31]}; species without stocks unfished")


def test_engine(cfg):
    m = fn.Model(small(cfg, days=1))
    nc = len(fn.CUE)                                                   # therm, food, risk, crowd, light, home, ...
    w = np.array(([4.0, 4.0, -3.0, 1.0, -1.0, 2.0] + [0.0] * nc)[:nc])
    calm = np.array(([1.0, 1.0, 0.0, 1.0, 0.0, 0.0] + [0.0] * nc)[:nc])   # comfortable, fed, safe, schooled, dark
    P = lambda here, hunger=0.0, ripe=0.0: m.drive_probs(here, w, np.float64(hunger), np.float64(ripe))
    bump = lambda k, x: calm + np.eye(nc)[k] * x
    base, danger, hot = P(calm), P(bump(2, 0.9)), P(bump(0, -0.8))
    ok = (base[fn.RANDOM] > 0.99 and P(calm, hunger=1)[fn.FORAGE] > 0.9 and danger[fn.FLEE] > 0.9
          and hot[fn.HABITAT] > 0.9 and P(calm, ripe=1)[fn.HOME] > 0.85 and P(bump(4, 1.0))[fn.HIDE] > 0.8)
    check("drives respond to stimuli", ok, f"P(flee|risk 0.9) {danger[fn.FLEE]:.2f}, P(forage|hungry) "
          f"{P(calm, hunger=1)[fn.FORAGE]:.2f}, P(random|content) {base[fn.RANDOM]:.2f}")
    F = np.zeros((3, len(fn.CUE))); F[:, 1] = [0.1, 0.9, 0.3]; F[:, 2] = [0.5, 0.1, 0.9]; F[:, 0] = [0.2, 0.3, 0.8]; F[:, 4] = [0.9, 0.5, 0.1]
    th, tv = m.targets(F, w, np.ones(3, bool), False), m.targets(F, w, np.ones(3, bool), True)
    ok = th[fn.FORAGE] == 1 and th[fn.FLEE] == 1 and th[fn.HABITAT] == 2 and th[fn.HIDE] == 0 and tv[fn.HIDE] == 2
    check("targets: food, safety, habitat, dark", ok, f"forage->{th[fn.FORAGE]} flee->{th[fn.FLEE]} "
          f"habitat->{th[fn.HABITAT]} hide->{tv[fn.HIDE]} (depth)")
    rng = np.random.default_rng(5)
    nb = len(fn.BEHAVIORS)
    Pr, tg, valid = rng.dirichlet(np.ones(nb), 50), rng.integers(0, 5, (50, nb - 1)), rng.random((50, 5)) > 0.3
    valid[:, 0] = True
    tg = np.where(valid[np.arange(50)[:, None], tg], tg, 0)
    stay = rng.dirichlet(np.ones(5), 50) * valid; stay /= stay.sum(1, keepdims=True)
    p = m.option_probs(Pr, tg, valid, stay)
    n = 40000
    b = (Pr.cumsum(1)[:, None] < rng.random((50, n, 1))).sum(2)
    ro = (np.where(valid, 1.0, 0).cumsum(1) / valid.sum(1, keepdims=True))[:, None] < rng.random((50, n, 1))
    so = (stay.cumsum(1)[:, None] < rng.random((50, n, 1))).sum(2)
    o = np.select([b == fn.RANDOM, b == fn.HOME], [ro.sum(2), so],
                  np.take_along_axis(tg, np.minimum(b, nb - 2), 1))
    emp = np.stack([(o == q).mean(1) for q in range(5)], 1)
    check("agent sampling = biomass mean field", np.abs(emp - p).max() < 0.015 and np.allclose(p.sum(1), 1),
          f"max |empirical - expected| {np.abs(emp - p).max():.3f} over 50 fish x 5 options")


def test_determinism(cfg):
    runs = [fn.Model(small(cfg, days=10, run__seed=s)).run() for s in (7, 7, 8)]
    d = [r.diagnostics()["biomass"] for r in runs]
    check("same seed -> identical run", np.array_equal(d[0], d[1]) and len(runs[0].ag) == len(runs[1].ag),
          f"agents {len(runs[0].ag)} / {len(runs[1].ag)}")
    check("different seed -> different run", not np.array_equal(d[0], d[2]), f"agents {len(runs[2].ag)}")


def test_restart(cfg):
    """A run stopped and resumed must write exactly what an uninterrupted run writes: from a mid-run restart (as after
    a crash, with later output on file to discard) and from the end state of a shorter run (extending it)."""
    def same(a, b):
        diff = []
        for f in ("series.nc", "fields.nc", "agents.nc", "lifehist.nc"):
            with nc.Dataset(a / f) as A, nc.Dataset(b / f) as B:
                diff += [f"{f}:{k}" for k in A.variables
                         if not np.array_equal(np.ma.filled(A[k][:], np.nan), np.ma.filled(B[k][:], np.nan), equal_nan=True)]
        return diff
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        c = small(cfg, days=12, run__restart_every_days=4, run__restart_keep=0, run__output_every_hours=48,
                  run__agent_output_every_hours=48)
        fn.Model(c, tmp / "full").run()
        fn.Model(c, tmp / "crash").run()
        fn.Model.resume(c, tmp / "crash", tmp / "crash/restart/restart_day00004.0.pkl", verbose=0).run()
        d = same(tmp / "full", tmp / "crash")
        check("resume after a crash = no crash", not d, "bit-identical outputs" if not d else f"differ: {d[:4]}")
        fn.Model(small(cfg, days=6, **{f"run__{k}": v for k, v in c["run"].items() if k != "days"}), tmp / "ext").run()
        fn.Model.resume(c, tmp / "ext", fn.latest_restart(tmp / "ext"), verbose=0).run()
        d = same(tmp / "full", tmp / "ext")
        check("extending a finished run", not d, "bit-identical outputs" if not d else f"differ: {d[:4]}")


def test_netcdf(cfg):
    c = small(cfg, days=6, ocean__source="synthetic")     # reference ocean must be a pure function of time
    m = fn.Model(c)
    g, days = m.g, np.arange(0, 12.0, 1.0)
    tmp = Path(tempfile.mkdtemp()) / "forcing.nc"
    with nc.Dataset(tmp, "w") as ds:
        for k, v in (("time", len(days)), ("depth", g.nz), ("lat", g.ny), ("lon", g.nx)):
            ds.createDimension(k, v)
        for k, v in (("time", days), ("depth", g.z), ("lat", g.lat), ("lon", g.lon)):
            ds.createVariable(k, "f8", (k,))[:] = v
        T, u, v = (ds.createVariable(k, "f8", ("time", "depth", "lat", "lon")) for k in ("temp", "u", "v"))
        for n, t in enumerate(days):
            f = m.ocean(t, 1e-6)
            T[n] = np.where(g.mask, f["T"], np.nan)
            u[n] = 0.5 * (f["U"][..., :-1] + f["U"][..., 1:]) / g.dy
            v[n] = 0.5 * (f["V"][:, :-1] / g.lx[:-1, None] + f["V"][:, 1:] / g.lx[1:, None])
    c2 = copy.deepcopy(c)
    c2["ocean"].update(source="netcdf", file=str(tmp))
    m2 = fn.Model(c2).run()
    T1, T2 = m.ocean(3.0, 1e-6)["T"], m2.ocean(3.0, 1e-6)["T"]
    err = np.abs(T1 - T2)[:, g.mask].max()
    check("netCDF forcing reproduces synthetic", err < 1e-6 and (m2.g.mask == g.mask).all(), f"max |dT| {err:.1e} degC")
    e = (m2.total_N() - m2.total0 - m2.ext) / m2.total0
    check("N budget with netCDF forcing", abs(e) < 1e-10, f"rel. error {e:+.2e}")

    # real products (GLORYS, HYCOM) have many more levels than the model, and the bed is wherever the
    # data run out: check that vertical grids that do not match are read and turned into bathymetry
    zf = np.linspace(2.0, g.z_e[-1], 3 * g.nz + 2)
    bed = np.where(np.arange(g.nx)[None, :] < g.nx // 3, 90.0, g.z_e[-1])       # a shelf in the west
    tmp2 = Path(tempfile.mkdtemp()) / "forcing_fine.nc"
    with nc.Dataset(tmp2, "w") as ds:
        for k, v in (("time", len(days)), ("depth", len(zf)), ("lat", g.ny), ("lon", g.nx)):
            ds.createDimension(k, v)
        for k, v in (("time", days), ("depth", zf), ("lat", g.lat), ("lon", g.lon)):
            ds.createVariable(k, "f8", (k,))[:] = v
        T, u, v = (ds.createVariable(k, "f8", ("time", "depth", "lat", "lon")) for k in ("temp", "u", "v"))
        for n, t in enumerate(days):
            f = m.ocean(t, 1e-6)
            k = np.argmin(np.abs(zf[:, None] - g.z[None, :]), 1)                # nearest model level per file level
            deep = (zf[:, None, None] <= bed[None]) & g.mask[None]
            T[n] = np.where(deep, f["T"][k], np.nan)
            uu = 0.5 * (f["U"][..., :-1] + f["U"][..., 1:]) / g.dy
            vv = 0.5 * (f["V"][:, :-1] / g.lx[:-1, None] + f["V"][:, 1:] / g.lx[1:, None])
            u[n], v[n] = np.where(deep, uu[k], np.nan), np.where(deep, vv[k], np.nan)
    c3 = copy.deepcopy(c)
    c3["ocean"].update(source="netcdf", file=str(tmp2))
    m3 = fn.Model(c3).run()
    west = m3.g.depth[g.mask & (np.arange(g.nx)[None, :] < g.nx // 3)]
    east = m3.g.depth[g.mask & (np.arange(g.nx)[None, :] >= g.nx // 3)]
    dz = zf[1] - zf[0]
    ok = west.size and east.size and abs(west.mean() - 90.0) <= dz and east.mean() > west.mean() + 100
    check("forcing levels need not match the model", bool(ok),
          f"{len(zf)} file levels -> {g.nz} model levels, inferred bed {west.mean():.0f} m on the shelf "
          f"(90 m set) and {east.mean():.0f} m offshore")


def test_consistency(cfg, days):
    for engine in ("probabilistic", "softmax"):
        B = {}
        for mode in ("none", "all"):
            c = small(cfg, days=days, hybrid__mode=mode, hybrid__agents_per_class=12, behavior__mode="fixed",
                      behavior__engine=engine)
            B[mode] = fn.Model(c).run().diagnostics()["biomass"].sum(1)
        rel = np.abs(B["all"] - B["none"]) / B["none"]
        check(f"IBM vs biomass ({engine})", rel.max() < 0.35,
              "rel. diff by species " + " ".join(f"{x:.2f}" for x in rel) + f" after {days} d")


def main():
    cfg = fn.load_config(next((a for a in sys.argv[1:] if not a.startswith("--")), "namelist.toml"))
    t0 = time.time()
    print("fishNET validation")
    for name, fn_ in (("conservation", lambda: test_budget(cfg)), ("transport", lambda: test_transport(cfg)),
                      ("genetics", test_genetics), ("aggregation", lambda: test_aggregation(cfg)),
                      ("generations", lambda: test_generations(cfg)), ("spin-up", lambda: test_spinup(cfg)),
                      ("warm start", lambda: test_warm_start(cfg)), ("initial state", lambda: test_initial_state(cfg)),
                      ("realism options", lambda: test_realism_options(cfg)), ("fishing map", lambda: test_fishing_map(cfg)),
                      ("behaviour engine", lambda: test_engine(cfg)), ("shallow water", lambda: test_shallow_water(cfg)),
                      ("skill", lambda: test_skill(cfg)), ("calibration", lambda: test_calibration(cfg)), ("multilocus genetics", lambda: test_multilocus(cfg)),
                      ("mating and isolation", lambda: test_mating(cfg)), ("diet names", lambda: test_diet_names(cfg)),
                      ("oxygen", lambda: test_oxygen(cfg)), ("metabolic index", lambda: test_metabolic_index(cfg)),
                      ("carbon export", lambda: test_export(cfg)), ("bathymetry", lambda: test_bathymetry(cfg)), ("cmems forcing", lambda: test_cmems(cfg)), ("ensemble", lambda: test_ensemble(cfg)), ("determinism", lambda: test_determinism(cfg)), ("restart", lambda: test_restart(cfg)), ("netcdf", lambda: test_netcdf(cfg)),
                      ("consistency", lambda: test_consistency(cfg, 90 if "--long" in sys.argv else 30))):
        print(f"- {name}")
        fn_()
    bad = [r for r in RESULTS if not r[1]]
    print(f"{len(RESULTS) - len(bad)}/{len(RESULTS)} checks passed in {time.time() - t0:.0f} s")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
