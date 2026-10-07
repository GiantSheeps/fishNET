#!/usr/bin/env python3
"""fishNET: hybrid Lagrangian-Eulerian marine ecosystem model.

Fish are tagged super-individuals inside IBM cells and stage-structured biomass (carrying
four moments of every heritable trait) everywhere else. Both representations feed, respire,
excrete, die and reproduce through the same equations and couple two-way to an NPZD ocean
driven by prescribed physics.

    python fishnet.py namelist.toml              # one run
    python fishnet.py namelist.toml --resume     # continue it from its newest restart file (or --resume FILE)
    python fishnet.py namelist.toml --ensemble   # parameter sweep from [ensemble]
"""
import os, sys, time, pickle, tomllib, hashlib, itertools, datetime as dtm
from pathlib import Path
import numpy as np
import netCDF4 as nc

R_EARTH, RHO_N, Q10, YOLK = 6.371e6, 1.8, 2.0, 0.5   # m; mmol N per g wet mass; metabolic Q10; egg yolk fraction
O2N, C2N = 10.6, 106 / 16                            # mmol O2 per mmol N respired; mol C per mol N (Redfield)
STAGES = ["egg", "larva", "juvenile", "adult"]
EGG, LARVA, JUV, ADULT = range(4)
TRAITS = ["T_opt", "T_width", "speed", "metab", "L_mat", "egg_mass",
          "w_therm", "w_food", "w_risk", "w_crowd", "w_light", "w_home", "w_shelf", "spawn_shift"]
T_OPT, T_WID, SPEED, METAB, L_MAT, EGG_M = range(6)
BEH = slice(6, 13)
SPAWN = 13                                   # days, heritable shift of the species spawning window
TMIN = np.array([-2.0, 0.5, 0.05, 0.2, 0.5, 1e-5] + [-np.inf] * 8)
TRAIT_DEFAULT = {"spawn_shift": [0.0, 0.0, 0.0]}   # traits a species file may omit: [mean, sd, h2]
PLANKTON = ["phyto", "zoo", "krill", "detritus"]
PLANKTON_L = np.array([0.005, 0.2, 2.5, 1.0])      # effective prey length (cm): phyto, meso-zoo, krill, detritus as benthos proxy
N_PLANK = len(PLANKTON)
N_NUTR = 5                                          # N, P, Z, K, D
O2_IDX = 5                                          # oxygen tracer index in NPZD.C
CAUSES = ["alive", "predation", "fishing", "natural", "thermal", "starvation", "hypoxia", "senescence", "merged"]
BEHAVIORS = ["forage", "seek_habitat", "flee", "school", "hide", "home", "seek_shelf", "random"]
FORAGE, HABITAT, FLEE, SCHOOL, HIDE, HOME, SHELF, RANDOM = range(8)
CUE = [1, 0, 2, 3, 4, 5, 6]       # cue (and drive-weight trait) behind each non-random behaviour
DIRS = np.array([[0, 0], [1, 0], [-1, 0], [0, 1], [0, -1]])   # stay, E, W, N, S
ACTIVITY = np.array([0.3, 1.0, 1.0, 1.0])                    # routine metabolism by stage
GAPE = 0.8                                                   # max fish prey/predator length ratio (class means)

# ---------------------------------------------------------------- biogeographic sectors
# Broad named regions so a species file can state where it starts without writing raw boxes.
# They are deliberately coarse: temperature, oxygen and bathymetry do the fine-grained work,
# and these only stop a species from being seeded in an ocean it has never lived in.
LAT_SECTORS = {"s_polar": [-90.0, -50.0], "s_temperate": [-50.0, -23.5], "tropical": [-23.5, 23.5],
               "n_temperate": [23.5, 50.0], "n_polar": [50.0, 90.0]}
LON_SECTORS = {"pacific_c": [-180.0, -135.0],   # central N/S Pacific, Hawaii, eastern Bering
               "pacific_e": [-135.0, -90.0],    # California Current, Gulf of Alaska (the Humboldt coast is east of 90W)
               "atlantic_w": [-90.0, -45.0],    # NW Atlantic shelf, Caribbean, Brazil, Patagonia
               "atlantic_e": [-45.0, 0.0],      # mid and NE Atlantic, Iberia, Rockall
               "africa": [0.0, 45.0],           # Mediterranean, Benguela, Agulhas, Red Sea
               "indian_w": [45.0, 90.0],        # Arabian Sea, western Indian Ocean
               "indian_e": [90.0, 135.0],       # Bay of Bengal, SE Asia, western Australia
               "pacific_w": [135.0, 180.0]}     # Kuroshio, Coral Sea, western Bering
DEPTH_SECTORS = {"surface": [0.0, 50.0], "epipelagic": [0.0, 200.0],
                 "mesopelagic": [200.0, 1000.0], "bathypelagic": [1000.0, 4000.0]}


def expand_sectors(s):
    """Turn `lat_sectors`/`lon_sectors`/`depth_sectors` in a species file into the raw boxes and
    `depth_range` the model already understands. Listing sectors on one axis and not the other means
    "all of the other axis". An explicit range always wins.

    Sectors set `seed_range`: where the species is placed at t=0, and nothing more. Fish are then free
    to swim, drift, spawn and shift their distribution anywhere the ocean allows, so ranges can move
    with climate and with the evolution of the movement traits. Set `range_is_hard = true` to write
    `native_range` instead, which is a wall nothing crosses for the whole run."""
    def look(table, keys, what):
        bad = [k for k in keys if k not in table]
        if bad:
            sys.exit(f"species '{s.get('name', '?')}': unknown {what} {bad}; choose from {sorted(table)}")
        return [table[k] for k in keys]
    lat = look(LAT_SECTORS, s.get("lat_sectors", []), "lat_sectors")
    lon = look(LON_SECTORS, s.get("lon_sectors", []), "lon_sectors")
    key = "native_range" if s.get("range_is_hard", False) else "seed_range"
    if (lat or lon) and key not in s:
        s[key] = [[lo0, lo1, la0, la1]                                 # one box per sector pair
                  for lo0, lo1 in (lon or [[-180.0, 180.0]])
                  for la0, la1 in (lat or [[-90.0, 90.0]])]
    dep = look(DEPTH_SECTORS, s.get("depth_sectors", []), "depth_sectors")
    if dep and "depth_range" not in s:
        s["depth_range"] = [min(d[0] for d in dep), max(d[1] for d in dep)]
    return s


def load_config(path):
    path = Path(path)
    cfg = tomllib.loads(path.read_text())
    find = lambda f: path.parent / f if (path.parent / f).exists() else path.parent / Path(f).name  # flat folders too
    cfg["species"] = [expand_sectors({**tomllib.loads(find(s["file"]).read_text()), **s}) for s in cfg["species"]]
    return cfg


def size_pref(ratio, beta, sigma):
    return np.exp(-0.5 * ((np.log(ratio) - np.log(beta)) / sigma) ** 2)


def key_name(x):
    """Species names are matched loosely: 'warm sardine', 'warm_sardine' and 'Warm-Sardine' are the same
    species in diet tables and parameter paths."""
    return "".join(c for c in str(x).lower() if c.isalnum())


def in_season(doy, window):
    d0, d1 = window[..., 0], window[..., 1]
    return np.where(d0 <= d1, (doy >= d0) & (doy <= d1), (doy >= d0) | (doy <= d1))


def softmax(x, axis):
    """Softmax that tolerates rows with no legal option (every entry -inf), returning zeros there rather
    than NaN. Such rows occur wherever the geometry leaves a cell with no reachable neighbour or depth."""
    mx = np.max(np.where(np.isfinite(x), x, -np.inf), axis=axis, keepdims=True)
    e = np.exp(np.where(np.isfinite(mx), x - np.where(np.isfinite(mx), mx, 0.0), -np.inf))
    tot = e.sum(axis=axis, keepdims=True)
    return np.where(tot > 0, e / np.where(tot > 0, tot, 1.0), 0.0)


def central_moments(S, N):
    """Mean and central moments 2-4 from power sums S[..., p-1] = sum x^p."""
    r = S / np.maximum(N, 1e-300)[..., None]
    m = r[..., 0]
    m2 = np.maximum(r[..., 1] - m ** 2, 0)
    m3 = r[..., 2] - 3 * m * r[..., 1] + 2 * m ** 3
    m4 = np.maximum(r[..., 3] - 4 * m * r[..., 2] + 6 * m ** 2 * r[..., 1] - 3 * m ** 4, 0)
    return m, m2, m3, m4


def raw_moments(m, m2, m3, m4):
    return np.stack([m, m2 + m ** 2, m3 + 3 * m * m2 + m ** 3, m4 + 4 * m * m3 + 6 * m ** 2 * m2 + m ** 4], -1)


class Genetics:
    """Multilocus architecture shared by all species: `loci` unlinked diploid loci whose allelic values are
    standardised, with a fixed pleiotropy matrix B mapping loci to traits. Breeding values are A = B (g1 + g2),
    so genetic correlations come from loci that affect several traits at once, and the additive genetic
    covariance G follows from the allele distribution rather than being imposed. Agents carry their alleles;
    biomass classes keep the moment (infinitesimal) description and are resampled on disaggregation."""

    def __init__(gen, c, h2, rng):
        gen.L = int(c.get("loci", 0))
        gen.rate, gen.sd = c.get("mutation_rate", 0.0), c.get("mutation_sd", 0.1)
        gen.h2 = h2                                        # (species, trait) heritabilities
        if not gen.L:
            return
        nt = h2.shape[1]
        R = np.eye(nt)
        for key, val in c.get("correlations", {}).items():
            a, b = [TRAITS.index(x.strip()) for x in key.split("-")]
            R[a, b] = R[b, a] = float(val)
        w, V = np.linalg.eigh(R)                           # nearest positive-definite correlation matrix
        R = V @ np.diag(np.maximum(w, 1e-6)) @ V.T
        d = np.sqrt(np.diag(R))
        gen.R = R / np.outer(d, d)
        Q = np.linalg.qr(rng.standard_normal((gen.L, nt)))[0].T       # nt x L with orthonormal rows
        gen.B = np.linalg.cholesky(gen.R) @ Q / np.sqrt(2)            # A = B (g1 + g2), unit-variance alleles
        gen.Binv = np.linalg.pinv(gen.B)

    def draw(gen, s, rng, n=None):
        """Allele values (rows, 2, L) for new fish of species s, giving Var(A) = h2 per trait."""
        n = len(s) if n is None else n
        return rng.standard_normal((n, 2, gen.L))

    def breeding(gen, gl, s):
        """Breeding values from alleles, scaled to each species' heritability."""
        return (gl.sum(1) @ gen.B.T) * np.sqrt(gen.h2[s])

    def match(gen, A, s, rng, k=1):
        """Alleles for target breeding values (used when biomass classes become individuals). Biomass classes
        store one moment per trait and no covariance, so the genotypes are drawn from the loci (which carry the
        genetic correlations) and then rescaled, per group of k siblings, to the stored mean and variance."""
        gl = gen.draw(s, rng, len(A))
        want = A / np.maximum(np.sqrt(gen.h2[s]), 1e-12)
        if k > 1:
            w, d = want.reshape(-1, k, want.shape[-1]), (gl.sum(1) @ gen.B.T).reshape(-1, k, want.shape[-1])
            z = (d - d.mean(1, keepdims=True)) / np.maximum(d.std(1, keepdims=True), 1e-12)
            want = (w.mean(1, keepdims=True) + z * w.std(1, keepdims=True)).reshape(want.shape)
        gl += 0.5 * ((want - gl.sum(1) @ gen.B.T) @ gen.Binv.T)[:, None, :]
        return gl

    def meiosis(gen, gm, gf, s, rng):
        """One gamete from each parent (free recombination), then mutation."""
        pick = lambda gp: np.take_along_axis(gp, rng.integers(0, 2, (len(gp), 1, gen.L)), 1)[:, 0]
        gl = np.stack([pick(gm), pick(gf)], 1)
        if gen.rate > 0:
            hit = rng.random(gl.shape) < gen.rate
            gl = gl + hit * rng.standard_normal(gl.shape) * gen.sd
        return gl

    fst_min = 20                                                     # agents a latitude band needs to count

    def fst(gen, gl, group, ngroup):
        """Wright's F_ST across groups (regions) from allele values: between-group variance over the total."""
        x = gl.mean(1)                                              # individual mean allelic value per locus
        if len(x) < 2 or ngroup < 2:
            return np.nan
        n = np.bincount(group, None, ngroup)[:, None]
        mu = np.bincount(np.repeat(group, gen.L), x.ravel(), ngroup * gen.L).reshape(ngroup, gen.L) / np.maximum(n, 1)
        use = (n >= gen.fst_min).ravel()                            # too few agents in a band say nothing about it
        if use.sum() < 2:
            return np.nan
        between = np.average(mu[use], 0, weights=n[use].ravel()).reshape(1, -1)
        vb = np.average((mu[use] - between) ** 2, 0, weights=n[use].ravel()).mean()
        vt = x.var(0).mean()
        return float(np.clip(vb / max(vt, 1e-12), 0, 1))

    def g_matrix(gen, A):
        """Realised additive genetic (co)variance matrix from a sample of breeding values."""
        return np.cov(A.T) if len(A) > len(TRAITS) else np.full((len(TRAITS), len(TRAITS)), np.nan)


def inherit(Am, Af, h2, rng):
    """Infinitesimal model for one offspring: midparent breeding value + segregation N(0, h2/2) (standardised units)."""
    return 0.5 * (Am + Af) + rng.standard_normal(np.shape(Am)) * np.sqrt(h2 / 2)


def offspring_moments(S, N, h2):
    """Infinitesimal model with random mating: offspring = midparent + N(0, h2/2) segregation."""
    m, v, m3, m4 = central_moments(S, N)
    s2 = h2 / 2
    return raw_moments(m, v / 2 + s2, m3 / 4, (m4 + 3 * v ** 2) / 8 + 3 * v * s2 + 3 * s2 ** 2)


def sample_moments(S, N, k, rng):
    """k samples per row matching mean and variance exactly, skew/kurtosis via Cornish-Fisher."""
    m, v, m3, m4 = [x[:, None, :] for x in central_moments(S, N[:, None])]
    sd = np.sqrt(v)
    g1 = np.clip(m3 / np.maximum(sd ** 3, 1e-12), -2, 2)
    g2 = np.clip(m4 / np.maximum(v ** 2, 1e-12) - 3, -1, 6)
    z = rng.standard_normal((S.shape[0], k, S.shape[1]))
    z = z + (z ** 2 - 1) * g1 / 6 + (z ** 3 - 3 * z) * g2 / 24 - (2 * z ** 3 - 5 * z) * g1 ** 2 / 36
    if k > 1:
        z = (z - z.mean(1, keepdims=True)) / np.maximum(z.std(1, keepdims=True), 1e-12)
    else:
        z[:] = 0
    return m + sd * z


class Species:
    """All species parameters packed into arrays indexed by species."""
    KEYS = ("L_inf L_juv lifespan_days egg_days T_ref cmax alpha K_food f_crit c_swim reserve_max kappa_R "
            "gonad_frac m0 m_early K_nursery m_thermal prey_size_ratio prey_size_sd larva_depth_max").split()

    def __init__(sp, specs, heritable=True, behavior="evolve"):
        sp.names, sp.n = [s["name"] for s in specs], len(specs)
        sp.color = [s.get("color", f"C{i}") for i, s in enumerate(specs)]
        for k in sp.KEYS:
            setattr(sp, k, np.array([float(s[k]) for s in specs]))
        sp.biomass = np.array([float(s.get("biomass", 0)) for s in specs])
        sp.crowd_L = np.array([float(s.get("crowd_length", np.nan)) for s in specs])
        sp.F = np.array([float(s.get("fishing_F", 0)) for s in specs])
        sp.lw_a, sp.lw_b = np.array([s["length_weight"] for s in specs], float).T
        sp.L50, sp.L95 = np.array([s.get("fishing_L50_L95", [1e9, 2e9]) for s in specs], float).T
        sp.depth = np.array([s["depth_range"] for s in specs], float)
        sp.demersal = np.array([bool(s.get("demersal", False)) for s in specs])
        sp.boxes = [np.array(s.get("native_range", []), float).reshape(-1, 4) for s in specs]   # lon0 lon1 lat0 lat1
        sp.seed_boxes = [np.array(s.get("seed_range", []), float).reshape(-1, 4) for s in specs]  # t=0 placement only
        sp.bottom_max = np.array([float(s.get("bottom_max_m", 1e9)) for s in specs])            # keep off the abyss
        sp.A_o = np.array([float(s.get("A_o", 0.4)) for s in specs])         # metabolic index: supply per kPa at W=1 g
        sp.E_o = np.array([float(s.get("E_o", 0.4)) for s in specs])         # eV, temperature sensitivity of demand
        sp.eps_o = np.array([float(s.get("eps_o", -0.1)) for s in specs])    # mass exponent of the index
        sp.phi_crit = np.array([float(s.get("phi_crit", 2.0)) for s in specs])   # index needed for active life
        sp.m_hypoxia = np.array([float(s.get("m_hypoxia", 0.2)) for s in specs])
        # length (cm) above which a fish stops eating phyto- and zooplankton (smooth cut-off); default: never
        sp.plankton_L_max = np.array([float(s.get("plankton_L_max", np.inf)) for s in specs])
        # adult crowding (g m-2 of same-species adults in a column at which adult natural mortality doubles)
        sp.K_adult = np.array([float(s.get("K_adult", np.inf)) for s in specs])
        sp.spawn_doy = np.array([s["spawn_doy"] for s in specs], float)
        sp.grounds = [np.array(s.get("spawning_grounds", []), float).reshape(-1, 2) for s in specs]
        tr = np.array([[s["traits"].get(t, TRAIT_DEFAULT[t]) if t in TRAIT_DEFAULT else s["traits"][t]
                        for t in TRAITS] for s in specs], float)
        sp.mu, sp.sd, sp.h2 = tr[..., 0], tr[..., 1], tr[..., 2]
        # length (cm) over which the juvenile crowding penalty falls by e; default a quarter of the maturity length
        sp.crowd_L = np.where(np.isnan(sp.crowd_L), 0.25 * sp.mu[:, L_MAT], sp.crowd_L)
        if not heritable:
            sp.h2[:] = 0
        if behavior == "fixed":
            sp.sd[:, BEH] = 0
        sp.W_inf = sp.weight(sp.L_inf)
        sp.rmet = sp.alpha * sp.cmax * sp.f_crit * sp.W_inf ** (2 / 3 - 0.8) / (1 + sp.c_swim)
        sp.ntype = N_PLANK + 4 * sp.n
        sp.taxo = np.zeros((sp.n, 4, sp.ntype))          # plankton by stage; fish prey filtered by size
        known, sp.diet_unmatched = {key_name(n) for n in sp.names}, []
        for i, s in enumerate(specs):
            sp.taxo[i, :2, :N_PLANK] = [s.get("larva_diet", s["diet"]).get(p, 0) for p in PLANKTON]
            sp.taxo[i, 2:, :N_PLANK] = [s["diet"].get(p, 0) for p in PLANKTON]
            diet = {key_name(k): v for k, v in s.get("fish_diet", {}).items()}
            for j, name in enumerate(sp.names):
                sp.taxo[i, :, N_PLANK + 4 * j: N_PLANK + 4 * (j + 1)] = diet.get(key_name(name), 0)
            # a prey key naming no species in this run would otherwise vanish in silence; the banner reports it
            sp.diet_unmatched += [(s["name"], k, v) for k, v in s.get("fish_diet", {}).items()
                                  if key_name(k) not in known and v]
            # Egg cannibalism: filter feeders / planktivores ingest planktonic eggs; conspecific eggs targeted
            is_planktivore = (sp.taxo[i, ADULT, 1] > 0) or (sp.taxo[i, ADULT, 2] > 0)
            if is_planktivore:
                for j in range(sp.n):
                    sp.taxo[i, JUV:, N_PLANK + 4 * j + EGG] = np.maximum(sp.taxo[i, JUV:, N_PLANK + 4 * j + EGG], 0.2)
            can_w = float(s.get("egg_cannibalism", 0.4 if is_planktivore else 0.0))
            if can_w > 0:
                sp.taxo[i, JUV:, N_PLANK + 4 * i + EGG] = np.maximum(sp.taxo[i, JUV:, N_PLANK + 4 * i + EGG], can_w)
        sp.taxo[:, EGG] = 0
        Legg, Lmat = sp.length(sp.mu[:, EGG_M]), sp.mu[:, L_MAT]
        sp.Lref = np.stack([Legg, (Legg + sp.L_juv) / 2, (sp.L_juv + Lmat) / 2, (Lmat + sp.L_inf) / 2], 1)
        sp.PREF = sp.taxo * sp.size_window(np.r_[PLANKTON_L, sp.Lref.ravel()][None, None] / sp.Lref[..., None],
                                           np.arange(sp.n)[:, None, None])
        sp.PREF[..., :N_PLANK - 1] *= sp.plankton_taper(sp.Lref, np.arange(sp.n)[:, None])[..., None]
        sp.RISK = sp.PREF.reshape(4 * sp.n, sp.ntype)[:, N_PLANK:].T   # [prey class, predator class]

    def plankton_taper(sp, L, s):
        """1 for small fish, falling to 0 around plankton_L_max: big fish do not live on copepods."""
        Lx = sp.plankton_L_max[s]
        with np.errstate(over="ignore", invalid="ignore"):
            return np.where(np.isfinite(Lx), 1 / (1 + np.exp((L - Lx) / np.maximum(0.1 * Lx, 1e-9))), 1.0)

    def size_window(sp, ratio, s):
        """Log-normal preference for prey/predator length ratio; fish prey (types N_PLANK+) are also gape-limited."""
        w = size_pref(ratio, sp.prey_size_ratio[s], sp.prey_size_sd[s])
        w[..., N_PLANK:] *= ratio[..., N_PLANK:] <= GAPE
        return w

    def _ab(sp, x, s):
        if s is None:
            return [v.reshape(-1, *[1] * (np.ndim(x) - 1)) for v in (sp.lw_a, sp.lw_b)]
        return sp.lw_a[s], sp.lw_b[s]

    def weight(sp, L, s=None):
        a, b = sp._ab(L, s)
        return a * L ** b

    def length(sp, W, s=None):
        a, b = sp._ab(W, s)
        return (np.maximum(W, 1e-12) / a) ** (1 / b)

    def box_mask(sp, g, boxes):
        """Columns inside any of each species' lon/lat boxes. No boxes at all means everywhere."""
        lon = ((g.lon[None, :] + 180) % 360) - 180
        out = np.zeros((sp.n, g.ny, g.nx), bool)
        for q, box in enumerate(boxes):
            if not len(box):
                out[q] = True
                continue
            for lo0, lo1, la0, la1 in box:
                inlon = (lon >= lo0) & (lon <= lo1) if lo0 <= lo1 else (lon >= lo0) | (lon <= lo1)
                out[q] |= inlon & (g.lat[:, None] >= la0) & (g.lat[:, None] <= la1)
        return out

    def home_range(sp, g):
        """Columns each species may occupy for its whole life: inside any `native_range` box and over water
        no deeper than `bottom_max_m`. A species with no boxes roams freely, which is the default -- the
        sectors in a species file set `seed_range` instead, which only places it at t=0."""
        return sp.box_mask(g, sp.boxes) & g.mask[None] & (g.depth[None] <= sp.bottom_max[:, None, None])

    def seed_range(sp, g):
        """Columns a species is seeded into at t=0. Nothing constrains it there afterwards."""
        return sp.box_mask(g, sp.seed_boxes) & g.mask[None]

    def levels_3d(sp, g):
        """Allowed levels per column (species, stage, nz, ny, nx): the species' depth band intersected with
        the water column and its home range, demersal stages pinned to the sea bed, and the shallowest wet
        cell as a fallback."""
        lev = sp.levels(g.z)[:, :, :, None, None] & g.wet[None, None] & sp.home_range(g)[:, None, None]
        dem = sp.demersal[:, None, None, None, None] & (np.arange(4) >= JUV)[None, :, None, None, None]
        lev = np.where(dem, lev & g.bed[None, None], lev)          # demersal stages live on the bed
        top = g.wet & ~np.roll(g.wet, 1, 0)                       # surface cell of each wet column
        top[0] |= g.wet[0]
        home = sp.home_range(g)[:, None, None]
        # A column with nothing in the species' depth band is simply unsuitable for a bed-dwelling species:
        # it cannot follow the sea floor into the abyss. Pelagic stages fall back to the surface cell.
        fallback = top[None, None] & g.wet[None, None] & home & ~sp.demersal[:, None, None, None, None]
        return np.where(lev.any(2, keepdims=True), lev, fallback)

    def levels(sp, z):
        """Allowed depth levels [species, stage, level]."""
        lev = np.zeros((sp.n, 4, len(z)), bool)
        lev[:, EGG, 0] = True
        lev[:, LARVA] = z[None] <= np.maximum(sp.larva_depth_max[:, None], z[0])
        ok = (z[None] >= sp.depth[:, :1]) & (z[None] <= sp.depth[:, 1:])
        nearest = z[None] == z[np.argmin(np.abs(z[None] - sp.depth[:, :1]), 1)][:, None]
        ok = np.where(ok.any(1, keepdims=True), ok, nearest)
        lev[:, JUV] = lev[:, ADULT] = ok
        return lev


class Grid:
    def __init__(g, c):
        g.nx, g.ny, g.periodic = c["nx"], c["ny"], c.get("periodic_x", False)
        g.lon_e, g.lat_e = np.linspace(*c["lon"], g.nx + 1), np.linspace(*c["lat"], g.ny + 1)
        g.lon, g.lat = 0.5 * (g.lon_e[1:] + g.lon_e[:-1]), 0.5 * (g.lat_e[1:] + g.lat_e[:-1])
        g.z_e = np.asarray(c["depth_edges"], float)
        g.z, g.dz, g.nz = 0.5 * (g.z_e[1:] + g.z_e[:-1]), np.diff(g.z_e), len(g.z_e) - 1
        g.dlon, g.dlat = g.lon_e[1] - g.lon_e[0], g.lat_e[1] - g.lat_e[0]
        rl = np.radians(g.dlon)
        g.area = np.repeat((R_EARTH ** 2 * rl * np.diff(np.sin(np.radians(g.lat_e))))[:, None], g.nx, 1)
        g.dx = np.repeat((R_EARTH * np.cos(np.radians(g.lat)) * rl)[:, None], g.nx, 1)
        g.dy = R_EARTH * np.radians(g.dlat)
        g.lx = R_EARTH * np.cos(np.radians(g.lat_e)) * rl        # zonal face lengths (ny+1)
        g.vol = g.dz[:, None, None] * g.area
        g.ncol, g.nc3 = g.ny * g.nx, g.nz * g.ny * g.nx
        g.flen = np.stack([np.full((g.ny, g.nx), g.dy), np.full((g.ny, g.nx), g.dy),
                           np.repeat(g.lx[1:, None], g.nx, 1), np.repeat(g.lx[:-1, None], g.nx, 1)])

    def set_bathymetry(g, depth):
        """Bottom depth per column (m). Cells are wet in proportion to the water they hold (partial bottom
        cells), so volumes, areas and the bed follow real topography instead of a flat-bottom box."""
        g.depth = np.where(np.isfinite(depth), np.clip(depth, 0.0, g.z_e[-1]), 0.0)
        frac = np.clip((g.depth[None] - g.z_e[:-1, None, None]) / g.dz[:, None, None], 0, 1)
        g.wet = frac > 1e-3                                  # (nz, ny, nx) water present
        g.vol = g.dz[:, None, None] * g.area * frac
        g.dzc = g.dz[:, None, None] * frac                   # water thickness in each cell
        g.volx = np.maximum(g.vol, 1e-9)                     # safe denominator (dry cells hold nothing)
        tx = np.zeros((g.nz, g.ny, g.nx + 1)); ty = np.zeros((g.nz, g.ny + 1, g.nx))
        thk = g.dz[:, None, None] * frac
        tx[..., 1:-1] = np.minimum(thk[..., 1:], thk[..., :-1])   # open thickness of each face
        ty[:, 1:-1] = np.minimum(thk[:, 1:], thk[:, :-1])
        if g.periodic:
            tx[..., 0] = tx[..., -1] = np.minimum(thk[..., 0], thk[..., -1])
        g.tx, g.ty = tx, ty
        pad = np.zeros((g.nz, g.ny + 2, g.nx + 2))
        pad[:, 1:-1, 1:-1] = thk
        if g.periodic:
            pad[:, 1:-1, 0], pad[:, 1:-1, -1] = thk[..., -1], thk[..., 0]
        g.tcorner = np.minimum(np.minimum(pad[:, 1:, 1:], pad[:, :-1, 1:]),
                               np.minimum(pad[:, 1:, :-1], pad[:, :-1, :-1]))   # (nz, ny+1, nx+1)
        g.dzx = np.maximum(g.dzc, 1e-9)
        g.kbot = np.maximum(g.wet.sum(0) - 1, 0)             # index of the deepest wet level
        g.bed = g.wet & ~np.roll(g.wet, -1, 0)               # deepest wet cell of each column
        g.bed[-1] |= g.wet[-1]
        g.set_mask(g.wet[0])

    def set_mask(g, mask):
        g.mask = mask
        if not hasattr(g, "wet"):
            g.set_bathymetry(np.where(mask, g.z_e[-1], 0.0))
            return
        jj, ii = np.mgrid[:g.ny, :g.nx]
        nj, ni = jj[None] + DIRS[:, 1, None, None], ii[None] + DIRS[:, 0, None, None]
        if g.periodic:
            ni %= g.nx
        ok = (nj >= 0) & (nj < g.ny) & (ni >= 0) & (ni < g.nx)
        g.nj, g.ni = np.clip(nj, 0, g.ny - 1), np.clip(ni, 0, g.nx - 1)
        g.nvalid = ok & mask[g.nj, g.ni] & mask[None]
        g.wetf = [g.wet & np.roll(g.wet, -d, ax) for ax, d in ((2, 1), (2, -1), (1, 1), (1, -1))]   # open faces E W N S
        if not g.periodic:
            g.wetf[0][..., -1] = g.wetf[1][..., 0] = False
        g.wetf[2][:, -1] = g.wetf[3][:, 0] = False

    def cell(g, lon, lat):
        i = np.floor((lon - g.lon_e[0]) / g.dlon).astype(int)
        j = np.floor((lat - g.lat_e[0]) / g.dlat).astype(int)
        return j, i

    def inside(g, j, i):
        ok = (j >= 0) & (j < g.ny) & (i >= 0) & (i < g.nx)
        return ok & g.mask[np.clip(j, 0, g.ny - 1), np.clip(i, 0, g.nx - 1)]


# Shallow-water scheme from:
# SHALLOW WATER MODEL
# Copyright (c) 2017 by Paul Connolly
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
def lax_wendroff(dx, dy, dt, g, u, v, h, u_tendency, v_tendency):
    """One Lax-Wendroff step of the shallow water equations (arrays indexed [x, y]; returns the interior)."""
    uh = u*h
    vh = v*h
    h_mid_xt = 0.5*(h[1:,:]+h[0:-1,:]) - (0.5*dt/dx)*(uh[1:,:]-uh[0:-1,:])
    h_mid_yt = 0.5*(h[:,1:]+h[:,0:-1]) - (0.5*dt/dy)*(vh[:,1:]-vh[:,0:-1])
    Ux = uh*u+0.5*g*h**2.
    Uy = uh*v
    uh_mid_xt = 0.5*(uh[1:,:]+uh[0:-1,:]) - (0.5*dt/dx)*(Ux[1:,:]-Ux[0:-1,:])
    uh_mid_yt = 0.5*(uh[:,1:]+uh[:,0:-1]) - (0.5*dt/dy)*(Uy[:,1:]-Uy[:,0:-1])
    Vx = Uy
    Vy = vh*v+0.5*g*h**2.
    vh_mid_xt = 0.5*(vh[1:,:]+vh[0:-1,:]) - (0.5*dt/dx)*(Vx[1:,:]-Vx[0:-1,:])
    vh_mid_yt = 0.5*(vh[:,1:]+vh[:,0:-1]) - (0.5*dt/dy)*(Vy[:,1:]-Vy[:,0:-1])
    h_new = h[1:-1,1:-1] - (dt/dx)*(uh_mid_xt[1:,1:-1]-uh_mid_xt[0:-1,1:-1]) \
        - (dt/dy)*(vh_mid_yt[1:-1,1:]-vh_mid_yt[1:-1,0:-1])
    Ux_mid_xt = uh_mid_xt*uh_mid_xt/h_mid_xt + 0.5*g*h_mid_xt**2.
    Uy_mid_yt = uh_mid_yt*vh_mid_yt/h_mid_yt
    uh_new = uh[1:-1,1:-1] - (dt/dx)*(Ux_mid_xt[1:,1:-1]-Ux_mid_xt[0:-1,1:-1]) \
        - (dt/dy)*(Uy_mid_yt[1:-1,1:]-Uy_mid_yt[1:-1,0:-1]) + dt*u_tendency*0.5*(h[1:-1,1:-1]+h_new)
    Vx_mid_xt = uh_mid_xt*vh_mid_xt/h_mid_xt
    Vy_mid_yt = vh_mid_yt*vh_mid_yt/h_mid_yt + 0.5*g*h_mid_yt**2.
    vh_new = vh[1:-1,1:-1] - (dt/dx)*(Vx_mid_xt[1:,1:-1]-Vx_mid_xt[0:-1,1:-1]) \
        - (dt/dy)*(Vy_mid_yt[1:-1,1:]-Vy_mid_yt[1:-1,0:-1]) + dt*v_tendency*0.5*(h[1:-1,1:-1]+h_new)
    return uh_new/h_new, vh_new/h_new, h_new


class ShallowWater:
    """Two active layers over a motionless abyss (2.5-layer reduced gravity) in a closed basin on a grid
    `sw_refine` times finer than the ecosystem grid. Driven by a double-gyre wind and by diabatic heating
    that moves mass between the layers (meridional gradient, seasonal and diurnal cycles)."""

    def __init__(sw, c, g, t0, doy0, hour0):
        r = c.get("sw_refine", 3)
        sw.c, sw.r, sw.ny, sw.nx, sw.doy0, sw.hour0, sw.t = c, r, g.ny * r, g.nx * r, doy0, hour0, t0
        lat = g.lat_e[0] + (np.arange(sw.ny) + 0.5) * g.dlat / r
        sw.lon = g.lon_e[0] + (np.arange(sw.nx) + 0.5) * g.dlon / r
        sw.dx = R_EARTH * np.cos(np.radians(g.lat.mean())) * np.radians(g.dlon) / r
        sw.dy = R_EARTH * np.radians(g.dlat) / r
        sw.f = 2 * 7.292e-5 * np.sin(np.radians(lat))[:, None]
        sw.yp = ((np.arange(sw.ny) + 0.5) / sw.ny)[:, None]
        sw.H = np.array([c.get("sw_H1", 400.0), c.get("sw_H2", 800.0)])
        sw.gp = (c.get("sw_gprime1", 0.05), c.get("sw_gprime2", 0.02))
        sw.tau = -c.get("sw_wind_stress", 0.12) * np.cos(2 * np.pi * sw.yp)          # N m-2, westerlies mid-basin
        sw.nl = c.get("sw_layers", 1)
        cmax = np.sqrt((sw.gp[0] + sw.gp[1] * (sw.nl == 2)) * sw.H[0] + sw.gp[1] * sw.H[1] * (sw.nl == 2)) + 1.5
        sw.dt = c.get("sw_cfl", 0.4) * min(sw.dx, sw.dy) / cmax
        shape = (2, sw.ny + 2, sw.nx + 2)
        sw.u, sw.v, sw.h = np.zeros(shape), np.zeros(shape), np.zeros(shape)
        noise = 1 + 0.02 * np.random.default_rng(c.get("sw_seed", 0)).standard_normal((sw.ny, sw.nx))
        sw.h[0, 1:-1, 1:-1] = sw.h1_eq(t0) * noise
        sw.h[1, 1:-1, 1:-1] = sw.H.sum() - sw.h[0, 1:-1, 1:-1]
        sw.vol = sw.h[:, 1:-1, 1:-1].sum()
        sw.walls()

    def h1_eq(sw, t):
        c, doy = sw.c, (sw.doy0 + t) % 365
        return sw.H[0] * (1 + c.get("sw_heating_gradient", 0.5) * (0.5 - sw.yp)) \
            * (1 + c.get("sw_heating_seasonal", 0.15) * np.sin(2 * np.pi * (doy - 135) / 365))

    def sun(sw, t, lon):
        """Diurnal insolation shape max(0, cos(hour angle)) in [0, 1]."""
        return np.maximum(0, np.cos(np.radians(15 * ((sw.hour0 + t * 24) % 24 + lon / 15 - 12))))

    def walls(sw):
        """Closed basin: zero-gradient thickness and no-slip walls (normal and tangential velocity reflected)."""
        for q, sx, sy in ((sw.h, 1, 1), (sw.u, -1, -1), (sw.v, -1, -1)):
            q[:, :, 0], q[:, :, -1] = sx * q[:, :, 1], sx * q[:, :, -2]
            q[:, 0], q[:, -1] = sy * q[:, 1], sy * q[:, -2]

    def step(sw, dt):
        c, h, u, v, (g1, g2) = sw.c, sw.h, sw.u, sw.v, sw.gp
        nu, rb = c.get("sw_viscosity", 900.0), c.get("sw_drag", 1e-7)
        ddx = lambda q: (q[1:-1, 2:] - q[1:-1, :-2]) / (2 * sw.dx)
        ddy = lambda q: (q[2:, 1:-1] - q[:-2, 1:-1]) / (2 * sw.dy)
        lap = lambda q: (q[1:-1, 2:] + q[1:-1, :-2] - 2 * q[1:-1, 1:-1]) / sw.dx ** 2 \
            + (q[2:, 1:-1] + q[:-2, 1:-1] - 2 * q[1:-1, 1:-1]) / sw.dy ** 2

        nl = sw.nl

        def layers(hp):                        # own-layer pressure in the flux, other layer's (from hp) as a tendency
            out = []
            for k in range(nl):
                gx = g2 * (nl == 2)
                ut, vt = -gx * ddx(hp[1 - k]) + nu * lap(u[k]), -gx * ddy(hp[1 - k]) + nu * lap(v[k])
                if k == 0 and nl == 1:                                          # interfacial drag
                    ut, vt = ut - rb * u[0, 1:-1, 1:-1], vt - rb * v[0, 1:-1, 1:-1]
                if k == 0:
                    ut = ut + sw.tau / (1025.0 * np.maximum(h[0, 1:-1, 1:-1], c.get("sw_wind_depth", 100.0)))
                else:
                    ut, vt = ut - rb * u[1, 1:-1, 1:-1], vt - rb * v[1, 1:-1, 1:-1]
                out.append([x.T for x in lax_wendroff(sw.dx, sw.dy, dt, (g1 + g2 * (nl == 2)) if k == 0 else g2,
                                                       u[k].T, v[k].T, h[k].T, ut.T, vt.T)])
            return out
        new = layers(h)                        # predictor, then corrector with cross-layer pressure at t + dt/2
        if nl == 2:
            hmid = np.pad(0.5 * (h[:, 1:-1, 1:-1] + np.stack([n[2] for n in new])), ((0, 0), (1, 1), (1, 1)), mode="edge")
            new = layers(hmid)
        th = sw.f * dt                          # exact inertial rotation (stable for any f dt)
        for k, (un, vn, hn) in enumerate(new):
            u[k, 1:-1, 1:-1] = un * np.cos(th) + vn * np.sin(th)
            v[k, 1:-1, 1:-1] = -un * np.sin(th) + vn * np.cos(th)
            h[k, 1:-1, 1:-1] = np.maximum(hn, c.get("sw_h_min", 20.0))
        he, te = c.get("sw_entrain_depth", 150.0), c.get("sw_entrain_days", 2.0) * 86400   # entrain before outcropping
        q = (sw.h1_eq(sw.t) - h[0, 1:-1, 1:-1]) / (c.get("sw_relax_days", 120.0) * 86400) \
            + c.get("sw_diurnal_w", 0.1) / 86400 * (np.pi * sw.sun(sw.t, sw.lon) - 1) \
            + (np.maximum(he - h[0, 1:-1, 1:-1], 0) - np.maximum(he - h[1, 1:-1, 1:-1], 0)) / te
        q = np.clip(q * dt, 10 - h[0, 1:-1, 1:-1], h[1, 1:-1, 1:-1] - 10)
        h[0, 1:-1, 1:-1] += q
        if nl == 2:                             # diabatic exchange between layers; restore any volume leak
            h[1, 1:-1, 1:-1] -= q
            h[1, 1:-1, 1:-1] += (sw.vol - h[:, 1:-1, 1:-1].sum()) / (sw.ny * sw.nx)
        else:                                   # 1.5 layers: the abyss is at rest, its upper interface flat
            h[1] = sw.H.sum() - h[0]
        sw.walls()
        sw.t += dt / 86400

    def advance(sw, t, verbose=False):
        report = sw.t
        while sw.t < t - 1e-9:
            sw.step(min(sw.dt, (t - sw.t) * 86400))
            if sw.t >= report and not np.isfinite(sw.h).all():
                sys.exit(f"shallow-water ocean became unstable near day {sw.t:.0f}: lower sw_cfl or raise sw_viscosity")
            if verbose and sw.t >= report:
                ke = 0.5 * (sw.u[0] ** 2 + sw.v[0] ** 2)[1:-1, 1:-1]
                print(f"  shallow-water spin-up day {sw.t:8.1f}: upper-layer mean KE {ke.mean():.2e} m2 s-2, "
                      f"max speed {np.sqrt(2 * ke.max()):.2f} m/s, h1 {sw.h[0].min():.0f}-{sw.h[0].max():.0f} m", flush=True)
            if sw.t >= report:
                report += 60

    def coarse(sw, q):
        return q[1:-1, 1:-1].reshape(sw.ny // sw.r, sw.r, sw.nx // sw.r, sw.r).mean((1, 3))


class Ocean:
    """Prescribed physics: temperature, face transports, vertical diffusivity and surface PAR."""

    def __init__(o, c, g, start, verbose=1):
        o.c, o.g, o.start, o.verbose = c, g, start, verbose
        o.doy0 = start.timetuple().tm_yday - 1 + (start.hour + start.minute / 60) / 24
        o.hour0 = start.hour + start.minute / 60
        xe = (g.lon_e - g.lon_e[0]) / (g.lon_e[-1] - g.lon_e[0])
        o.ye = (g.lat_e - g.lat_e[0]) / (g.lat_e[-1] - g.lat_e[0])
        o.yc = 0.5 * (o.ye[1:] + o.ye[:-1])
        if c["source"] == "netcdf":
            o._load()
        elif c.get("geometry") == "global":
            o._build_global()
        else:
            w = c.get("coast_cells", 0)
            land = np.zeros((g.ny, g.nx), bool)
            if w:
                edge = g.nx - w - np.round(1.5 * np.sin(2 * np.pi * 1.5 * o.yc)).astype(int)
                land = np.arange(g.nx)[None] >= edge[:, None]
            shelf = c.get("shelf_cells", 0)
            dep = np.full((g.ny, g.nx), g.z_e[-1], float)
            if w and shelf:                                  # continental shelf and slope off the east coast
                dist = (edge[:, None] - np.arange(g.nx)[None]).astype(float)
                dep = np.clip(c.get("shelf_depth", 100.0) * (1 + 9 * np.clip((dist - shelf) / max(shelf, 1), 0, 1) ** 2),
                              0, g.z_e[-1])
            g.set_bathymetry(np.where(land, 0.0, dep))
            shape = (1 - xe) * (1 - np.exp(-xe / c.get("wbc_width", 0.05)))
            psi = c.get("gyre_psi", 3e4) * np.sin(2 * np.pi * o.ye)[:, None] * (shape / shape.max())[None]
            corner = np.zeros((g.ny + 1, g.nx + 1), bool)
            for dj in (0, 1):
                for di in (0, 1):
                    corner[dj:dj + g.ny, di:di + g.nx] |= land
            psi[corner] = 0
            decay = np.exp(-g.z / c.get("current_decay_m", 500))[:, None, None]
            o.psi3d = psi[None] * decay                  # corner streamfunction per level (m2/s)
            coast = np.where(land.any(1), np.argmax(land, 1), 10 ** 6)
            o.dcoast = (coast[:, None] - np.arange(g.nx)[None] - 0.5) * g.dx / 1e3
        o.upwell_zone = np.exp(-np.maximum(o.dcoast, 0) / c.get("upwelling_width_km", 150)) if hasattr(o, "dcoast") \
            else np.zeros((g.ny, g.nx))
        o.projector()
        if c["source"] == "shallow_water":
            o.sw = ShallowWater(c, g, -c.get("sw_spinup_days", 720.0), o.doy0, o.hour0)
            key = hashlib.md5(repr((sorted((k, v) for k, v in c.items() if k.startswith("sw_")), g.nx, g.ny,
                                    list(g.lon_e[[0, -1]]), list(g.lat_e[[0, -1]]), o.doy0)).encode()).hexdigest()[:12]
            cache = Path(c.get("sw_cache_dir", "sw_cache")) / f"spinup_{key}.npz"
            if cache.exists():
                z = np.load(cache)
                o.sw.h, o.sw.u, o.sw.v, o.sw.t = z["h"], z["u"], z["v"], float(z["t"])
            else:
                if verbose:
                    print(f"shallow-water ocean: {o.sw.nx}x{o.sw.ny}x2 layers, dt {o.sw.dt / 60:.0f} min, "
                          f"spinning up {-o.sw.t:.0f} days (cached afterwards in {cache})", flush=True)
                o.sw.advance(0.0, verbose)
                cache.parent.mkdir(parents=True, exist_ok=True)
                np.savez(cache, h=o.sw.h, u=o.sw.u, v=o.sw.v, t=o.sw.t)
            o.projector()

    def _build_global(o):
        """Global ocean: real coastlines, depth from distance offshore, and a wind-driven surface circulation
        (trades, westerlies and a circumpolar current) projected onto the geometry."""
        from scipy.ndimage import distance_transform_edt, label
        c, g = o.c, o.g
        try:
            from global_land_mask import globe
        except ImportError:
            sys.exit("geometry = 'global' needs a land mask: pip install global-land-mask")
        LON, LAT = np.meshgrid(g.lon, g.lat)
        land = globe.is_land(np.clip(LAT, -89.9, 89.9), ((LON + 180) % 360) - 180)
        if c.get("fill_lakes", True):                        # drop basins too small to matter at this resolution
            lab, n = label(~land, structure=np.ones((3, 3)))
            if n:
                size = np.bincount(lab.ravel())
                size[0] = 0
                land |= (lab != size.argmax()) & ~land
        cells = distance_transform_edt(~land, sampling=(g.dy / 1e3, g.dx.mean() / 1e3))   # km offshore
        shelf, slope = c.get("shelf_km", 150.0), c.get("slope_km", 250.0)
        deep, shallow = min(c.get("abyss_m", 4000.0), g.z_e[-1]), c.get("shelf_depth", 150.0)
        dep = shallow + (deep - shallow) * np.clip((cells - shelf) / slope, 0, 1) ** 2
        g.set_bathymetry(np.where(land, 0.0, dep))
        o.dcoast = np.where(land, 0.0, cells)
        o.projector()
        psi = o.sverdrup()
        o.psi3d = psi[None] * np.exp(-g.z / c.get("current_decay_m", 500))[:, None, None]

    def sverdrup(o):
        """Wind-driven circulation as a corner streamfunction. Sverdrup transport is integrated westward from
        each basin's eastern boundary, so every basin gets its own subtropical and subpolar gyres, and closed
        by a western boundary layer; a circumpolar term is added on top."""
        c, g = o.c, o.g
        phi = np.radians(g.lat)
        tau = -c.get("tau0", 0.1) * np.cos(3 * np.clip(phi, -np.pi / 3, np.pi / 3))   # trades, westerlies, polar
        curl = -np.gradient(tau, g.dy)                                               # zonal stress only
        beta = 2 * 7.292e-5 * np.cos(phi) / R_EARTH
        V = curl / (1025.0 * np.maximum(beta, 1e-13))                                # m2/s meridional transport
        wet, psi, lwbc = g.mask, np.zeros((g.ny, g.nx)), max(c.get("wbc_cells", 2.0), 0.5)
        for j in range(g.ny):
            row = wet[j]
            if row.all() or not row.any():                                           # circumpolar or dry
                continue
            order = (np.arange(g.nx) + np.argmax(~row)) % g.nx
            seg, run = [], []
            for i in order:
                if row[i]:
                    run.append(i)
                elif run:
                    seg.append(run)
                    run = []
            if run:
                seg.append(run)
            for sgm in seg:                                                      # sgm runs west to east
                acc = 0.0
                for i in reversed(sgm):                                          # psi = 0 at the eastern wall
                    psi[j, i] = acc
                    acc -= V[j] * g.dx[j, 0]                                     # psi(x) = -int_x^xe V dx'
        psi /= max(c.get("current_decay_m", 500.0), 1.0)     # depth-integrated transport -> per unit depth
        if c.get("acc_u", 0.0):
            band = np.exp(-((g.lat - c.get("acc_lat", -55.0)) / c.get("acc_width", 10.0)) ** 2)
            psi = psi - np.cumsum(c["acc_u"] * band * g.dy)[:, None]
        u = -np.gradient(psi, g.dy, axis=0) * wet                                # interior flow implied by psi
        v = np.gradient(psi, axis=1) / g.dx * wet
        return o.streamfunction(u, v)                       # projected: zero on coasts, gyres closed smoothly

    def projector(o):
        """Least-squares streamfunction fit (psi = 0 on coasts and walls) giving non-divergent face transports."""
        from scipy.sparse import coo_matrix
        from scipy.sparse.linalg import splu
        g = o.g
        pad = np.pad(~g.mask, 1, constant_values=True)
        o.free = ~(pad[1:, 1:] | pad[:-1, 1:] | pad[1:, :-1] | pad[:-1, :-1])
        idx = np.full(o.free.shape, -1)
        idx[o.free] = np.arange(o.free.sum())
        o.mf = np.zeros((g.ny, g.nx + 1), bool); o.mf[:, 1:-1] = g.mask[:, 1:] & g.mask[:, :-1]
        o.zf = np.zeros((g.ny + 1, g.nx), bool); o.zf[1:-1] = g.mask[1:] & g.mask[:-1]
        jm, im = np.nonzero(o.mf); jz, iz = np.nonzero(o.zf)
        nm = len(jm)
        rows = np.r_[np.arange(nm), np.arange(nm), nm + np.arange(len(jz)), nm + np.arange(len(jz))]
        cidx = np.r_[idx[jm + 1, im], idx[jm, im], idx[jz, iz + 1], idx[jz, iz]]
        vals = np.r_[-np.ones(nm), np.ones(nm), np.ones(len(jz)), -np.ones(len(jz))]
        use = cidx >= 0
        o.A = coo_matrix((vals[use], (rows[use], cidx[use])), shape=(nm + len(jz), o.free.sum())).tocsr()
        o.lu = splu((o.A.T @ o.A).tocsc())

    def streamfunction(o, uc, vc):
        """psi (ny+1, nx+1) whose face transports best match the cell-centred velocities uc, vc."""
        g = o.g
        U = np.zeros((g.ny, g.nx + 1)); V = np.zeros((g.ny + 1, g.nx))
        U[:, 1:-1] = 0.5 * (uc[:, 1:] + uc[:, :-1]) * g.dy
        V[1:-1] = 0.5 * (vc[1:] + vc[:-1]) * g.lx[1:-1, None]
        psi = np.zeros(o.free.shape)
        psi[o.free] = o.lu.solve(o.A.T @ np.r_[U[o.mf], V[o.zf]])
        return psi

    def _load(o):
        """Read forcing from netCDF. Variable and coordinate names are auto-detected (CMEMS/GLORYS, ROMS and
        plain names) unless `ocean.names` overrides them; longitude convention, latitude order, packed
        integers and the surrounding box are all handled here."""
        import xarray as xr
        c, g = o.c, o.g
        ds = xr.open_mfdataset(c["file"], combine="by_coords") if any(ch in str(c["file"]) for ch in "*?[") \
            else xr.open_dataset(c["file"])
        alt = dict(temp=["temp", "thetao", "votemper", "TEMP", "temperature", "water_temp"],
                   u=["u", "uo", "vozocrtx", "u_eastward", "water_u", "UVEL"],
                   v=["v", "vo", "vomecrty", "v_northward", "water_v", "VVEL"],
                   lon=["lon", "longitude", "nav_lon", "x"], lat=["lat", "latitude", "nav_lat", "y"],
                   depth=["depth", "deptht", "lev", "z", "depthu"], time=["time", "time_counter"])
        names = {}
        for k, cand in alt.items():
            given = c.get("names", {}).get(k)
            found = given or next((n for n in cand if n in ds.variables or n in ds.dims), None)
            if found is None:
                sys.exit(f"netcdf forcing: no variable found for '{k}' (tried {cand}); set ocean.names")
            names[k] = found
        ds = ds.rename({names[k]: k for k in ("lon", "lat", "depth", "time") if names[k] in ds.dims or names[k] in ds.coords})
        if "depth" not in ds.dims:                              # surface-only forcing: one level
            ds = ds.expand_dims(depth=[0.0])
        lon = ds["lon"].values
        if float(np.nanmax(lon)) > 180 and g.lon_e[0] < 0:      # 0..360 file, -180..180 grid
            ds = ds.assign_coords(lon=(((lon + 180) % 360) - 180)).sortby("lon")
        elif float(np.nanmin(lon)) < 0 and g.lon_e[-1] > 180:
            ds = ds.assign_coords(lon=(lon % 360)).sortby("lon")
        for k in ("lat", "depth"):                              # CMEMS ships descending latitude
            if ds[k].size > 1 and float(ds[k][1] - ds[k][0]) < 0:
                ds = ds.isel({k: slice(None, None, -1)})
        pad = 2 * max(g.dlon, g.dlat)
        ds = ds.sel(lon=slice(g.lon_e[0] - pad, g.lon_e[-1] + pad), lat=slice(g.lat_e[0] - pad, g.lat_e[-1] + pad))
        if min(ds.sizes["lon"], ds.sizes["lat"]) < 2:
            sys.exit("netcdf forcing: the file does not cover the model grid")
        ds = ds.sel(depth=slice(0, g.z_e[-1] + 200)).load()
        t = ds["time"].values
        o.tf = (t - np.datetime64(o.start)) / np.timedelta64(1, "D") if np.issubdtype(t.dtype, np.datetime64) \
            else t.astype(float)
        if o.verbose and o.tf[-1] - o.tf[0] > 400 and (o.tf[0] > 31 or o.tf[-1] < 0):   # a dated multi-year file
            print(f"  warning: the forcing covers days {o.tf[0]:.0f}..{o.tf[-1]:.0f} after run.start, so the run "
                  "starts outside it and the file is looped: check run.start", flush=True)
        src = lambda v: ds[names[v]].transpose("time", "depth", "lat", "lon")
        wet = xr.DataArray(np.isfinite(src("temp").values[0, 0]).astype(float), coords={"lat": ds.lat, "lon": ds.lon})
        if not np.isfinite(src("temp").values[0, 0]).any():
            sys.exit("netcdf forcing: the first time/level of the temperature field is entirely missing")
        wetc = wet.interp(lat=g.lat, lon=g.lon, method="nearest", kwargs={"fill_value": None}).values > 0.5
        if o.c.get("bathymetry"):                            # depth from a file variable (e.g. deptho)
            bz = xr.open_dataset(o.c["bathymetry"]) if isinstance(o.c["bathymetry"], str) else ds
            bz = bz.rename({n: k for k in ("lon", "lat") for n in alt[k] if n in bz.dims or n in bz.coords})
            for k in ("lat", "lon"):                          # match the forcing conventions
                if bz[k].size > 1 and float(bz[k][1] - bz[k][0]) < 0:
                    bz = bz.isel({k: slice(None, None, -1)})
            if float(bz.lon.max()) > 180 and g.lon_e[0] < 0:
                bz = bz.assign_coords(lon=(((bz.lon.values + 180) % 360) - 180)).sortby("lon")
            bv = bz[o.c.get("bathymetry_var", "deptho")]
            dep = bv.interp(lat=g.lat, lon=g.lon, kwargs={"fill_value": None}).values
        else:                                                # infer from the deepest level with valid data
            col = np.isfinite(src("temp").values[0]).astype(float)
            kk = xr.DataArray(col, dims=("depth", "lat", "lon"),
                              coords={"depth": ds[names["depth"]], "lat": ds[names["lat"]], "lon": ds[names["lon"]]})
            kk = kk.interp(lat=g.lat, lon=g.lon, method="nearest", kwargs={"fill_value": None}).values > 0.5
            zf = np.asarray(ds["depth"].values, float)       # the file's own levels, not the model's
            dep = np.where(kk.any(0), zf[np.clip(kk.sum(0) - 1, 0, len(zf) - 1)], 0.0)
        g.set_bathymetry(np.where(wetc, np.nan_to_num(dep), 0.0))

        def get(v):                                # fill land/below-bottom first so coastal cells interpolate cleanly
            a, out = src(v), []
            for i in range(0, a.sizes["time"], 12):   # a year of months at a time: fill_nan makes ~7 temporary
                b = a.isel(time=slice(i, i + 12))     # copies, which OOMs a 20-year 1-deg file loaded whole
                x = b.values.copy()
                for k in range(1, x.shape[1]):
                    x[:, k] = np.where(np.isnan(x[:, k]), x[:, k - 1], x[:, k])
                out.append(b.copy(data=fill_nan(x)).interp(lon=g.lon, lat=g.lat, depth=g.z,
                                                           kwargs={"fill_value": None}).values)
            return np.concatenate(out)
        o.Tf, uf, vf = get("temp"), get("u"), get("v")
        o.Tf, uf, vf = [np.nan_to_num(a) for a in (o.Tf, uf, vf)]
        uE = 0.5 * (uf + np.roll(uf, -1, -1)) * g.dy
        o.Uf = np.concatenate([np.roll(uE, 1, -1)[..., :1], uE], -1)
        vN = 0.5 * (vf[..., :-1, :] + vf[..., 1:, :]) * g.lx[None, None, 1:-1, None]
        z = np.zeros_like(vN[..., :1, :])
        o.Vf = np.concatenate([z, vN, z], -2)
        o.U, o.V = o.Uf[0].copy(), o.Vf[0].copy()

    def __call__(o, t, dt):
        c, g = o.c, o.g
        doy = (o.doy0 + t) % 365
        z = g.z[:, None, None]
        if c["source"] == "netcdf":
            span = o.tf[-1] - o.tf[0] + (o.tf[1] - o.tf[0] if len(o.tf) > 1 else 1)
            tt = o.tf[0] + (t - o.tf[0]) % span
            i1 = np.searchsorted(o.tf, tt, "right") % len(o.tf)
            i0 = (i1 - 1) % len(o.tf)
            dt_seg = o.tf[i1] - o.tf[i0] if i1 != 0 else o.tf[0] + span - o.tf[-1]   # the wrap-around gap
            w = np.clip((tt - o.tf[i0]) / max(dt_seg, 1e-9), 0, 1)
            T = (1 - w) * o.Tf[i0] + w * o.Tf[i1]
            U = (1 - w) * o.Uf[i0] + w * o.Uf[i1]
            V = (1 - w) * o.Vf[i0] + w * o.Vf[i1]
            uc = 0.5 * (U[..., 1:] + U[..., :-1]) / g.dy                  # face transports -> cell velocities
            vc = 0.5 * (V[:, 1:] / g.lx[1:, None] + V[:, :-1] / g.lx[:-1, None])
            psi3d = np.stack([o.streamfunction(uc[k], vc[k]) for k in range(g.nz)])
            upwell = np.zeros((g.ny, g.nx))
        else:
            seas = np.sin(2 * np.pi * (doy - 135) / 365)
            if c.get("geometry") == "global":                 # symmetric about the equator, seasons out of phase
                la = np.radians(g.lat)
                sst = c.get("T_equator", 28.0) - (c.get("T_equator", 28.0) - c.get("T_polar", -1.5)) * np.sin(la) ** 2 \
                    + c["T_season"] * np.sin(la) * np.abs(np.sin(la)) * seas
            else:
                sst = c["T_south"] - (c["T_south"] - c["T_north"]) * o.yc + c["T_season"] * (0.4 + o.yc) * seas
            upwell = o.upwell_zone * max(0.0, np.sin(2 * np.pi * (doy - 100) / 365))
            sst = sst[:, None] - c.get("upwelling_dT", 0) * upwell
            mld = c["mld_min"] + (c["mld_max"] - c["mld_min"]) * 0.5 * (1 + np.cos(2 * np.pi * (doy - 60) / 365))
            T = c["T_deep"] + (sst - c["T_deep"]) * np.exp(-np.maximum(z - mld, 0) / c.get("thermocline_m", 150))
            psi3d = o.psi3d
            if c["source"] == "shallow_water":
                T, psi3d, upwell = o.shallow(t, doy, sst, upwell, mld, z)
        T = np.maximum(T + c.get("warming_per_year", 0.0) * t / 365 * np.exp(-z / 300), -1.8)   # seawater freezes
        for hw in c.get("heatwave", []):
            if hw["start_day"] <= t < hw["start_day"] + hw["days"]:
                d = haversine(g.lon[None], g.lat[:, None], hw["lon"], hw["lat"])
                ramp = np.sin(np.pi * (t - hw["start_day"]) / hw["days"]) ** 2
                T = T + hw["amplitude"] * ramp * np.exp(-(d / hw["radius_km"]) ** 2) * np.exp(-z / 100)
        grad = np.maximum(T[:-1] - T[1:], 0) / np.diff(g.z)[:, None, None]
        Kz = c.get("kz_bg", 1e-5) + c.get("kz_ml", 0.02) * np.exp(-grad / c.get("kz_grad", 0.005))
        Psi = psi3d * g.tcorner                                   # m3/s: zero wherever a corner cell is dry
        Qx, Qy = -(Psi[:, 1:] - Psi[:, :-1]), Psi[:, :, 1:] - Psi[:, :, :-1]   # non-divergent by construction
        U, V = Qx / np.maximum(g.tx, 1e-9), Qy / np.maximum(g.ty, 1e-9)        # m2/s per unit depth
        divH = (Qx[..., 1:] - Qx[..., :-1]) + (Qy[:, 1:, :] - Qy[:, :-1, :])
        Qw = np.zeros((g.nz + 1, g.ny, g.nx))
        Qw[1:] = np.cumsum(divH, axis=0)
        Wface = Qw / np.maximum(g.area, 1e-9)
        W = 0.5 * (Wface[:-1] + Wface[1:]) * g.wet
        return dict(T=T, U=U, V=V, W=W, Qx=Qx, Qy=Qy, Kz=Kz, par=o.par(t, dt), doy=doy, upwell=upwell)

    def shallow(o, t, doy, sst, upwell, mld, z):
        """Temperature, currents and nutrient upwelling from the stacked shallow-water ocean."""
        c, g, sw = o.c, o.g, o.sw
        sw.advance(t)
        h1, h2 = sw.coarse(sw.h[0]), sw.coarse(sw.h[1])
        h1eq = sw.h1_eq(t)[:, 0].reshape(g.ny, sw.r).mean(1)[:, None]
        sst = sst + c.get("sw_dT_thickness", 6.0) * (h1 - h1eq) / sw.H[0]       # warm rings deep, cold rings shallow
        upwell = np.maximum(upwell, np.clip((h1eq - h1) / (c.get("sw_upwelling_scale", 0.3) * sw.H[0]), 0, 1))
        T2 = 0.5 * (c["T_south"] - (c["T_south"] - c["T_north"]) * o.yc[:, None] + c["T_deep"])
        sig = lambda x: 0.5 * (1 + np.tanh(x))
        T = c["T_deep"] + (T2 - c["T_deep"]) * sig((h1 + h2 - z) / 100) \
            + (sst - T2) * np.exp(-np.maximum(z - mld, 0) / c.get("thermocline_m", 150)) * sig((h1 - z) / 30) \
            + c.get("diurnal_dT", 0.5) * sw.sun(t, g.lon)[None, None] * np.exp(-z / 10)
        pc = lambda x: np.pad(x, 1, mode="edge")
        h1c, h2c = [0.25 * (p[1:, 1:] + p[:-1, 1:] + p[1:, :-1] + p[:-1, :-1]) for p in (pc(h1), pc(h1 + h2))]
        ov = lambda a, b: np.clip(np.minimum(g.z_e[1:, None, None], b) - np.maximum(g.z_e[:-1, None, None], a), 0, None) \
            / g.dz[:, None, None]
        psi = ov(0, h1c) * o.streamfunction(sw.coarse(sw.u[0]), sw.coarse(sw.v[0])) \
            + ov(h1c, h2c) * o.streamfunction(sw.coarse(sw.u[1]), sw.coarse(sw.v[1]))
        return T, psi, upwell

    def mask_faces(o, U, V):
        nv = o.g.nvalid
        U[..., 1:] *= nv[1]
        U[..., 0] *= nv[2][:, 0]
        V[..., 1:, :] *= nv[3]
        V[..., 0, :] *= nv[4][0]
        return U, V

    def par(o, t, dt):
        """Surface PAR (W m-2) averaged over [t, t+dt] from solar geometry."""
        g = o.g
        ts = t + (np.arange(8) + 0.5) / 8 * dt
        doy = o.doy0 + ts
        decl = np.radians(23.44) * np.sin(2 * np.pi * (284 + doy) / 365)
        hour = (o.hour0 + ts * 24) % 24
        h = np.radians(((hour[:, None] + g.lon[None] / 15) - 12) * 15)[:, None, :]
        phi = np.radians(g.lat)[None, :, None]
        cosz = np.sin(phi) * np.sin(decl)[:, None, None] + np.cos(phi) * np.cos(decl)[:, None, None] * np.cos(h)
        return 0.43 * 1361 * 0.7 * np.maximum(cosz, 0).mean(0)


def fill_nan(x):
    """Extend values into NaN cells by repeated 4-neighbour averaging over the last two axes."""
    for _ in range(sum(x.shape[-2:])):
        bad = np.isnan(x)
        if not bad.any():
            break
        p = np.pad(x, [(0, 0)] * (x.ndim - 2) + [(1, 1), (1, 1)], constant_values=np.nan)
        nb = np.stack([p[..., 1:-1, :-2], p[..., 1:-1, 2:], p[..., :-2, 1:-1], p[..., 2:, 1:-1]])
        cnt = (~np.isnan(nb)).sum(0)
        x = np.where(bad & (cnt > 0), np.nansum(nb, 0) / np.maximum(cnt, 1), x)
    return x


def metabolic_index(O2, T, W, A_o, E_o, eps_o):
    """Deutsch et al. (2015): ratio of oxygen supply to resting demand. phi < 1 cannot be sustained,
    and active behaviour (feeding, spawning, swimming) needs phi above a species' phi_crit."""
    pO2 = np.maximum(O2, 0) / np.maximum(o2_saturation(T), 1e-9) * 21.0      # kPa partial pressure
    kB, Tk = 8.617e-5, T + 273.15
    return A_o * np.maximum(W, 1e-9) ** eps_o * pO2 * np.exp(E_o / kB * (1 / Tk - 1 / 288.15))   # demand rises with T


def o2_saturation(T, S=35.0):
    """Oxygen solubility in equilibrium with air (mmol m-3), Garcia and Gordon (1992)."""
    Ts = np.log((298.15 - T) / (273.15 + T))
    A = [5.80871, 3.20291, 4.17887, 5.10006, -0.0986643, 3.80369]
    B = [-0.00701577, -0.00770028, -0.0113864, -0.00951519]
    lnC = sum(a * Ts ** i for i, a in enumerate(A)) + S * sum(b * Ts ** i for i, b in enumerate(B)) - 2.75915e-7 * S ** 2
    return np.exp(lnC) * 1.025                      # umol/kg -> mmol m-3


def schmidt_o2(T):
    """Schmidt number for oxygen in seawater (Wanninkhof 2014)."""
    return 1920.4 - 135.6 * T + 5.2122 * T ** 2 - 0.10939 * T ** 3 + 0.00093777 * T ** 4


def haversine(lon1, lat1, lon2, lat2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 2 * R_EARTH / 1e3 * np.arcsin(np.sqrt(np.minimum(a, 1)))


def find_file(path):
    """A data file named in a namelist: as given (relative to the working directory), else next to fishnet.py."""
    p = Path(path)
    return p if p.exists() or p.is_absolute() else (Path(__file__).resolve().parent / p)


def woa_field(path, g, var="o_an", scale=1.025):
    """A World Ocean Atlas annual field (e.g. WOA23 dissolved oxygen, umol/kg) averaged onto the model grid:
    layer-weighted in depth, box-averaged in the horizontal, then gaps (coasts, cells deeper than the
    data) filled from wet neighbours and from the level above. `scale` converts umol/kg to mmol m-3."""
    with nc.Dataset(find_file(path)) as D:
        D.set_auto_mask(False)
        o = np.asarray(D[var][0], float)
        fill = getattr(D[var], "_FillValue", None)
        o = np.where((o == fill) | (np.abs(o) > 1e20), np.nan, o) if fill is not None else np.where(np.abs(o) > 1e20, np.nan, o)
        lat, lon = np.asarray(D["lat"][:], float), np.asarray(D["lon"][:], float)
        zb = np.asarray(D["depth_bnds"][:], float) if "depth_bnds" in D.variables else None
        z = np.asarray(D["depth"][:], float)
    if zb is None:
        mid = 0.5 * (z[1:] + z[:-1])
        zb = np.stack([np.r_[z[0], mid], np.r_[mid, z[-1] + (z[-1] - mid[-1])]], 1)
    ov = np.clip(np.minimum(g.z_e[1:, None], zb[None, :, 1]) - np.maximum(g.z_e[:-1, None], zb[None, :, 0]), 0, None)
    ok = np.isfinite(o)
    Oz = np.einsum("kl,lji->kji", ov, np.nan_to_num(o)) / np.maximum(np.einsum("kl,lji->kji", ov, ok * 1.0), 1e-12)
    Oz = np.where(np.einsum("kl,lji->kji", ov, ok * 1.0) > 0, Oz, np.nan)          # (nz, lat, lon) of the file
    lon = (lon - g.lon_e[0]) % 360 + g.lon_e[0] if g.periodic or g.lon_e[-1] - g.lon_e[0] > 300 else lon
    jj = np.floor((lat - g.lat_e[0]) / g.dlat).astype(int)
    ii = np.floor((lon - g.lon_e[0]) / g.dlon).astype(int)
    J, I = np.meshgrid(jj, ii, indexing="ij")
    inside = (J >= 0) & (J < g.ny) & (I >= 0) & (I < g.nx)
    tot, cnt = np.zeros((g.nz, g.ny, g.nx)), np.zeros((g.nz, g.ny, g.nx))
    for k in range(g.nz):
        v = Oz[k]
        w = inside & np.isfinite(v)
        np.add.at(tot[k], (J[w], I[w]), v[w])
        np.add.at(cnt[k], (J[w], I[w]), 1)
    out = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)
    for _ in range(max(g.nx, g.ny)):                  # grow into wet cells the data miss, level by level
        gap = g.wet & ~np.isfinite(out)
        if not gap.any():
            break
        pad = np.pad(out, ((0, 0), (1, 1), (1, 1)), constant_values=np.nan)
        if g.periodic:
            pad[:, :, 0], pad[:, :, -1] = pad[:, :, -2], pad[:, :, 1]
        nb = np.stack([pad[:, 1:-1, :-2], pad[:, 1:-1, 2:], pad[:, :-2, 1:-1], pad[:, 2:, 1:-1]])
        n = np.isfinite(nb).sum(0)
        grow = gap & (n > 0)
        if not grow.any():
            break
        out[grow] = (np.nansum(nb, 0) / np.maximum(n, 1))[grow]
    for k in range(1, g.nz):                          # anything still missing takes the level above
        out[k] = np.where(np.isfinite(out[k]), out[k], out[k - 1])
    return out * scale


class NPZD:
    """Nutrient-based NPZD (mmol N m-3) as C[0:5] = N, P, Z, K, D, C[5] = O2 on the (nz, ny, nx) grid."""

    def __init__(p, c, g):
        p.c, p.g = c, g
        z = g.z[:, None, None] * np.ones((1, g.ny, g.nx))
        N = c["N_deep"] - (c["N_deep"] - c["N_surf"]) * np.exp(-z / 150)
        K0 = c.get("K0", c["Z0"] * 0.5)
        p.C = np.stack([N, c["P0"] * np.exp(-z / 50), c["Z0"] * np.exp(-z / 50), K0 * np.exp(-z / 80),
                        c["D0"] * np.exp(-z / 100),
                        o2_saturation(np.full_like(z, c.get("T_init_O2", 12.0))) * (1 - c.get("O2_utilised", 0.35)
                        * (1 - np.exp(-z / c.get("O2_scale_m", 250.0))))])
        p.O2_ref = None                              # deep oxygen the sea bed relaxes to (None = O2_deep)
        if str(c.get("O2_init", "profile")).lower() == "woa":
            O = woa_field(c.get("O2_file", "observations/woa/woa23_all_o00_01.nc"), g)
            p.C[O2_IDX] = np.where(np.isfinite(O), O, p.C[O2_IDX])
            p.O2_ref = p.C[O2_IDX].copy()
        p.C *= g.wet
        p.h = np.diff(g.z)

    def light(p, par):
        k = p.c["kw"] + p.c["kc"] * p.C[1]
        tau = k * p.g.dzx
        top = np.cumsum(tau, 0) - tau
        return par * np.exp(-top) * (1 - np.exp(-tau)) / tau

    def biology(p, T, I, dt):
        c = p.c
        N, P, Z, K, D, O = p.C
        fP, fZ = 1.066 ** (T - 20), Q10 ** ((T - 20) / 10)
        hyp = np.maximum(O, 0) / (c.get("kO2_remin", 8.0) + np.maximum(O, 0))   # respiration slows as oxygen runs out
        up = c["mu_max"] * fP * N / (c["kN"] + N) * I / np.sqrt(I ** 2 + c["Ik"] ** 2) * P
        gr = c["g_max"] * fZ * P ** 2 / (c["kP"] ** 2 + P ** 2) * Z
        pm = c["mP"] * P
        ze, zq = c["mZ"] * fZ * Z * hyp, c["mZ2"] * Z ** 2
        g_max_k = c.get("g_max_k", 0.6 * c["g_max"])
        kP_k = c.get("kP_k", c["kP"])
        gr_kp = g_max_k * fZ * (P ** 2 / (kP_k ** 2 + P ** 2 + 1e-12)) * K
        gr_kz = c.get("g_max_kz", 0.4 * c["g_max"]) * fZ * (Z ** 2 / (kP_k ** 2 + Z ** 2 + 1e-12)) * K
        ke = c.get("mK", 0.04) * fZ * K * hyp
        kq = c.get("mK2", 0.02) * K ** 2
        beta_k = c.get("beta_k", c["beta"])
        rm = c["rD"] * fZ * D * hyp
        lim = lambda X, *f: [x * np.minimum(1, X / np.maximum(dt * sum(f), 1e-30)) for x in f]
        up, = lim(N, up)
        gr, gr_kp, pm = lim(P, gr, gr_kp, pm)
        gr_kz, ze, zq = lim(Z, gr_kz, ze, zq)
        ke, kq = lim(K, ke, kq)
        rm, = lim(D, rm)
        dN = ze + ke + rm - up
        dP = up - gr - gr_kp - pm
        dZ = c["beta"] * gr - ze - zq - gr_kz
        dK = beta_k * (gr_kp + gr_kz) - ke - kq
        dD = (1 - c["beta"]) * gr + (1 - beta_k) * (gr_kp + gr_kz) + pm + zq + kq - rm
        p.C[:5] += dt * np.stack([dN, dP, dZ, dK, dD])
        p.C[O2_IDX] = np.maximum(p.C[O2_IDX] + dt * O2N * (up - ze - ke - rm), 0)   # photosynthesis makes it, respiration spends it

    def air_sea(p, T, wind, dt):
        """Oxygen exchange with the atmosphere in the surface layer; returns the flux into the ocean (mmol)."""
        g = p.g
        kw = 0.251 * wind ** 2 * (660 / schmidt_o2(T[0])) ** 0.5 / 100 * 24        # m/d (Wanninkhof 2014)
        dO = (o2_saturation(T[0]) - p.C[O2_IDX, 0]) * np.minimum(1, kw * dt / g.dzx[0]) * g.wet[0]
        p.C[O2_IDX, 0] += dO
        return (dO * g.vol[0]).sum()

    def vdiff(p, Kz, dt):
        """Implicit, conservative vertical diffusion (Thomas algorithm over all columns)."""
        g, C = p.g, p.C
        a = np.zeros_like(C); b = np.zeros_like(C)
        face = g.wet[1:] & g.wet[:-1]                        # interior interfaces that hold water
        coef = dt * Kz * face / (0.5 * (g.dzx[1:] + g.dzx[:-1]))
        a[:, 1:] = -coef / g.dzx[1:]
        b[:, :-1] = -coef / g.dzx[:-1]
        d = 1 - a - b
        cp, dp = np.empty_like(C), np.empty_like(C)
        cp[:, 0], dp[:, 0] = b[:, 0] / d[:, 0], C[:, 0] / d[:, 0]
        for k in range(1, g.nz):
            den = d[:, k] - a[:, k] * cp[:, k - 1]
            cp[:, k], dp[:, k] = b[:, k] / den, (C[:, k] - a[:, k] * dp[:, k - 1]) / den
        C[:, -1] = dp[:, -1]
        for k in range(g.nz - 2, -1, -1):
            C[:, k] = dp[:, k] - cp[:, k] * C[:, k + 1]

    def sink(p, dt):
        """Upwind detritus sinking; returns export out of the bottom (mmol)."""
        g, D = p.g, p.C[4]
        n = int(np.ceil(p.c["w_sink"] * dt / np.maximum(g.dzc[g.wet].min(), 1e-3) / 0.5))
        export, keep = 0.0, not p.c.get("export_bottom", True)
        kref = int(np.searchsorted(g.z_e[1:], p.c.get("export_depth", 100.0)))
        p.export_ref = 0.0
        for _ in range(n):
            flux = D * p.c["w_sink"] * dt / n                       # mmol m-2
            D -= flux / g.dzx
            D[1:] += (flux[:-1] * g.wet[1:]) / g.dzx[1:]            # into the cell below, if there is one
            if kref < g.nz:
                p.export_ref += (flux[kref] * g.area * g.wet[kref]).sum()   # sinking past the reference depth
            hit = flux * g.bed                                      # reaching the sea bed
            export += (hit * g.area).sum()
            if keep:
                D += hit / g.dzx
        return 0.0 if keep else export

    def restore(p, upwell, dt):
        """Deep-water and coastal-upwelling nutrient supply; returns net input (mmol)."""
        c, g = p.c, p.g
        dN = np.zeros_like(p.C[0])
        if c.get("restore_days", 0) > 0:
            dN += (c["N_deep"] - p.C[0]) * (dt / c["restore_days"]) * g.bed
            ref = c.get("O2_deep", 220.0) if getattr(p, "O2_ref", None) is None else p.O2_ref
            dO = (ref - p.C[O2_IDX]) * (dt / c["restore_days"]) * g.bed * g.wet
            p.C[O2_IDX] += dO
            p.o2_in = (dO * g.vol).sum()
        if c.get("upwelling_days", 0) > 0:
            up = (g.z < c.get("upwelling_depth", 100))[:, None, None] * upwell[None]
            dN += (c["N_deep"] - p.C[0]) * np.minimum(1, up * dt / c["upwelling_days"])
        dN *= g.wet
        p.C[0] += dN
        return (dN * g.vol).sum()


def exchange(M, f, periodic):
    """Conservatively move fractions f = (E, W, N, S) of amounts M (..., ny, nx) to neighbours."""
    out = [M * fi for fi in f]
    new = M - sum(out)
    new[..., :, 1:] += out[0][..., :, :-1]
    new[..., :, :-1] += out[1][..., :, 1:]
    new[..., 1:, :] += out[2][..., :-1, :]
    new[..., :-1, :] += out[3][..., 1:, :]
    if periodic:
        new[..., :, 0] += out[0][..., :, -1]
        new[..., :, -1] += out[1][..., :, 0]
    return new


def pick(rng, keys, query, k=None):
    """For each query, a random index i with keys[i] == query (or -1 if none).
    With k, returns k independent draws per query (shape (len(query), k)): the candidates a chooser compares."""
    shape = (len(query),) if k is None else (len(query), k)
    if not len(keys):
        return np.full(shape, -1)
    order = np.argsort(keys, kind="stable")
    ks = keys[order]
    lo, hi = np.searchsorted(ks, query, "left"), np.searchsorted(ks, query, "right")
    if k is not None:
        lo, hi = lo[:, None], hi[:, None]
    idx = np.minimum(lo + np.floor(rng.random(shape) * (hi - lo)).astype(int), len(ks) - 1)
    return np.where(hi > lo, order[idx], -1)


class Agents:
    """Structure-of-arrays container for tagged super-individuals (n = fish represented)."""
    INT = "id sp st sex k gen mother father origin beh hold".split()
    FLT = "lon lat n W E G age dev amax tb lon0 lat0 base r0".split()

    def __init__(a, ntr, nbeh, nloci=0):
        d = {k: np.zeros(0, np.int64) for k in a.INT} | {k: np.zeros(0) for k in a.FLT}
        d |= dict(A=np.zeros((0, ntr)), Z=np.zeros((0, ntr)), dl=np.zeros((0, nbeh)), grad=np.zeros((0, nbeh)))
        if nloci:
            d["gl"] = np.zeros((0, 2, nloci))
        object.__setattr__(a, "d", d)
        object.__setattr__(a, "next_id", 0)

    def __getattr__(a, k):
        if k == "d" or k.startswith("__"):                   # unpickling probes these before `d` exists
            raise AttributeError(k)
        return a.d[k]

    def __setattr__(a, k, v):
        if k in a.d:
            a.d[k] = np.asarray(v, a.d[k].dtype)
        else:
            object.__setattr__(a, k, v)

    def __len__(a):
        return len(a.d["n"])

    def add(a, **cols):
        m = len(cols["n"])
        cols["id"] = np.arange(a.next_id, a.next_id + m)
        a.next_id += m
        for k, v in a.d.items():
            new = np.broadcast_to(np.asarray(cols.get(k, 0), v.dtype), (m, *v.shape[1:]))
            a.d[k] = np.concatenate([v, new])

    def keep(a, mask):
        for k in a.d:
            a.d[k] = a.d[k][mask]


class Model:
    def __init__(m, cfg, out_dir=None, verbose=None):
        m.cfg, r = cfg, cfg["run"]
        m.hy, m.be = cfg["hybrid"], cfg["behavior"]
        m.rng = np.random.default_rng(r.get("seed", 0))
        m.verbose = r.get("verbose", 1) if verbose is None else verbose
        m.start = dtm.datetime.fromisoformat(r.get("start", "2000-01-01T00:00"))
        m.dt = r["dt_hours"] / 24
        m.dt_s = m.dt * 86400
        m.nsteps = int(round(r["days"] / m.dt))
        m.g = g = Grid(cfg["grid"])
        m.ocean = Ocean(cfg["ocean"], g, m.start, m.verbose)
        if str(r.get("initial_biomass", "namelist")).lower() == "literature":
            literature_biomass(cfg, m.verbose)
        m.sp = S = Species(cfg["species"], cfg.get("genetics", {}).get("heritable", True), m.be["mode"])
        m.npzd = NPZD(cfg["npzd"], g)
        m.HOME = S.home_range(g)
        m.SEED = S.seed_range(g)                 # where each species starts; not a constraint afterwards
        m.shelf = 1.0 / (1.0 + g.depth / m.be.get("shelf_scale_m", 400.0))      # 1 inshore, ~0 over the abyss
        m.LEV = S.levels_3d(g)
        m.gdir = m.ground_directions()
        m.has_ground = np.array([len(x) > 0 for x in S.grounds])
        m.gen = Genetics(cfg.get("genetics", {}), S.h2, m.rng)
        m.ag = Agents(len(TRAITS), BEH.stop - BEH.start, m.gen.L)
        m.EN, m.EB = np.zeros((S.n, 4, g.ny, g.nx)), np.zeros((S.n, 4, g.ny, g.nx))
        m.ER = np.zeros_like(m.EB)             # reserve part of EB (EB = structure + reserve)
        m.EG = np.zeros((S.n, g.ny, g.nx))
        m.ES = np.zeros((S.n, 4, len(TRAITS), 4, g.ny, g.nx))
        m.EGEN = np.zeros_like(m.EN)           # fish x generation per class, so lineages survive biomass <-> agent swaps
        # Juveniles in biomass cells carry `juvenile_bins` progress bins toward maturity (numbers per bin; they sum to
        # EN[:, JUV]). With one bin a fixed fraction matures each step whatever its age, so some fish mature almost at
        # once; K bins make the time to maturity Erlang(K) with the same mean.
        m.NJ = np.zeros((S.n, max(1, int(m.hy.get("juvenile_bins", 1))), g.ny, g.nx))
        m.Fv = m.LEV / np.maximum(m.LEV.sum(2, keepdims=True), 1e-300)
        m.ibm, m.bref = np.zeros((g.ny, g.nx), bool), np.full(S.n, np.nan)
        m.t, m.step_i, m.ext, m.o2_ext, m.obits = 0.0, 0, 0.0, 0.0, []
        m.spin = None                            # (length, cycle) in days while spinning up, else None
        m.exp_sink = m.exp_fish = m.exp_resp = m.exp_bed = 0.0
        m.reset_step_tallies()
        m.init_fish()
        if r.get("warm_start"):
            m.warm_start(r["warm_start"], r.get("warm_start_plankton", True))
        m.force()
        m.fields()
        m.vertical()
        m.regulate()
        m.total0 = m.total_N()
        m.out = Output(m, out_dir) if out_dir else None

    # ------------------------------------------------------------------ setup
    def ground_directions(m):
        """Where ripe fish head. With `spawning_grounds` listed in a species file this is the nearest one;
        otherwise it is the local gradient of ripe conspecific biomass, so spawning aggregations form and
        persist by themselves. Agents additionally home to their own birth place (natal fidelity)."""
        g, S = m.g, m.sp
        gd = np.zeros((S.n, 2, g.ny, g.nx))
        lon, lat = np.meshgrid(g.lon, g.lat)
        for s, gr in enumerate(S.grounds):
            if len(gr):
                dx = (gr[:, 0, None, None] - lon) * np.cos(np.radians(lat))
                dy = gr[:, 1, None, None] - lat
                near = np.argmin(dx ** 2 + dy ** 2, 0)[None]
                v = np.stack([np.take_along_axis(dx, near, 0)[0], np.take_along_axis(dy, near, 0)[0]])
                gd[s] = v / np.maximum(np.hypot(*v), 0.5)
        return gd

    def spawn_attraction(m):
        """Unit vector up the gradient of spawning-ready conspecific biomass (smoothed), per species."""
        g, S = m.g, m.sp
        B = m.EB[:, ADULT] + m.EG
        if len(m.ag):
            a = m.ag
            ad = a.st == ADULT
            j, i = g.cell(a.lon[ad], a.lat[ad])
            B = B + np.bincount((a.sp[ad] * g.ny + j) * g.nx + i, a.n[ad] * (a.W + a.E + a.G)[ad],
                                S.n * g.ncol).reshape(B.shape)
        f = np.log1p(B / np.maximum(g.area, 1.0))
        for _ in range(m.be.get("spawn_smooth", 2)):                       # sense the neighbourhood, not one cell
            f = 0.25 * (np.roll(f, 1, -1) + np.roll(f, -1, -1) + np.roll(f, 1, -2) + np.roll(f, -1, -2))
            f = f * g.mask
        dy, dx = np.gradient(f, axis=(-2, -1))
        v = np.stack([dx, dy], 1)
        return v / np.maximum(np.hypot(v[:, 0], v[:, 1]), 1e-9)[:, None]

    def init_fish(m):
        S, g = m.sp, m.g
        T = m.ocean(0.0, m.dt)["T"]
        for s in range(S.n):
            phi = np.exp(-0.5 * ((T - S.mu[s, T_OPT]) / S.mu[s, T_WID]) ** 2)
            suit = (phi * m.LEV[s, ADULT]).max(0) ** 2 * g.mask * m.SEED[s]
            if S.biomass[s] > 0 and not (suit > 0).any():
                print(f"  warning: {S.names[s]} has no suitable cell in its seed_range (check its sectors, "
                      f"depth_range and T_opt {S.mu[s, T_OPT]:.1f}+-{S.mu[s, T_WID]:.1f} degC); it starts empty")
            Btot = S.biomass[s] * suit / max((suit * g.area).sum() / (g.area * g.mask).sum(), 1e-30) * g.area
            for st, frac in ((JUV, 0.3), (ADULT, 0.7)):
                m.EB[s, st] = frac * Btot
                m.EN[s, st] = frac * Btot / (S.weight(S.Lref[s, st], s) * (1 + 0.5 * S.reserve_max[s]))
                m.ER[s, st] = m.EB[s, st] * 0.5 * S.reserve_max[s] / (1 + 0.5 * S.reserve_max[s])   # half-full
        m.NJ[:] = m.EN[:, JUV, None] / m.NJ.shape[1]          # start spread evenly over the progress bins
        h2 = S.h2[:, None, :, None, None]
        m.ES[:, :, :, 1] = m.EN[:, :, None] * h2
        m.ES[:, :, :, 3] = 3 * m.EN[:, :, None] * h2 ** 2

    def fishing_rate(m, j, i, s):
        """Fishing mortality (1/yr) where each fish is. [fishing] mode = "constant" (default): each species'
        fishing_F everywhere. mode = "ram": RAM Legacy F for the calendar year (or the fixed `year`), per FAO major
        area, so only the areas with assessed stocks are fished (`unassessed` x the species' mean F elsewhere),
        multiplied by `scale`. Years outside the table take the nearest year it has."""
        fc = m.cfg.get("fishing", {})
        if str(fc.get("mode", "constant")).lower() != "ram":
            return m.sp.F[s]
        year = int(fc.get("year", 0)) or (m.start + dtm.timedelta(days=float(m.t))).year
        if getattr(m, "Fyear", None) != year:
            m.Fmap, m.Fyear = m.fishing_map(year, fc), year
        return m.Fmap[s, j, i]

    def fishing_map(m, year, fc):
        import json
        import obsval                                        # FAO areas live with the observation code
        g, S = m.g, m.sp
        if not hasattr(m, "_ram_F"):
            path = find_file(fc.get("file", "observations/ram_fishing_series.json"))
            m._ram_F = json.loads(Path(path).read_text())["species"]
            m._fao = obsval.fao_area(*np.meshgrid(g.lon, g.lat)) * g.mask
        out = np.zeros((S.n, g.ny, g.nx))
        scale, rest = float(fc.get("scale", 1.0)), float(fc.get("unassessed", 0.0))
        for q, name in enumerate(S.names):
            e = m._ram_F.get(name)
            if not e:
                continue
            out[q] = rest * (e.get("mean") or 0.0) * g.mask
            for a, rec in e["areas"].items():
                ys = np.array([int(y) for y in rec["years"]])
                y = ys[np.argmin(np.abs(ys - year))]                # nearest year with data for this area
                out[q][m._fao == int(a)] = rec["years"][str(y)]
        if m.verbose:
            fished = [f"{n} {out[q][g.mask].max():.2f}" for q, n in enumerate(S.names) if out[q].any()]
            print(f"  fishing {year}: RAM Legacy F x {scale:g} by FAO area (max F per species: {', '.join(fished)})",
                  flush=True)
        return out * scale

    def warm_start(m, src, plankton=True):
        """Start from the fish (and, if `plankton`, the NPZD and oxygen) of an earlier run's restart file, with
        the clock back at 0. `src` is a restart .pkl or a run folder (its newest restart). Species are matched
        by name, so the namelist may add, drop or retune species: new ones keep their fresh start. The grid
        must match. Agents keep their ages, genes and generations; birth days shift to the new clock. With
        npzd.O2_init = "woa", oxygen comes from the World Ocean Atlas rather than the old run."""
        path = Path(src)
        path = latest_restart(path) if path.is_dir() else path
        st = load_restart(path)
        old, day = st["model"], st["day"]
        g0 = old["g"]
        if (g0.nx, g0.ny, g0.nz) != (m.g.nx, m.g.ny, m.g.nz) or not np.allclose(g0.z_e, m.g.z_e):
            sys.exit(f"warm_start: {path} is on a {g0.nx}x{g0.ny}x{g0.nz} grid, this run is {m.g.nx}x{m.g.ny}x{m.g.nz}")
        names = list(old["sp"].names)
        pairs = [(q, names.index(n)) for q, n in enumerate(m.sp.names) if n in names]
        new, got = [q for q, _ in pairs], [o for _, o in pairs]
        for k in ("EN", "EB", "ER", "EG", "ES", "EGEN"):
            if k in old:
                getattr(m, k)[new] = old[k][got]
            elif k == "EGEN":
                m.EGEN[new] = 0
        K = m.NJ.shape[1]                                   # progress bins: carried over if the count matches
        if "NJ" in old and old["NJ"].shape[1] == K:
            m.NJ[new] = old["NJ"][got]
        else:
            m.NJ[new] = m.EN[new, JUV, None] / K
        m.bref[new] = old["bref"][got]
        a0 = old["ag"]
        remap = np.full(len(names), -1)
        remap[got] = new
        same = set(a0.d) == set(m.ag.d) and all(a0.d[k].shape[1:] == m.ag.d[k].shape[1:] for k in a0.d)
        if same and len(a0):
            keep = remap[a0.sp] >= 0
            m.ag.d = {k: v[keep].copy() for k, v in a0.d.items()}
            m.ag.sp = remap[m.ag.sp]
            m.ag.tb = m.ag.tb - day
            m.ag.next_id = a0.next_id
        elif len(a0) and m.verbose:
            print("  warm start: the agents' columns differ (genetics settings changed?); they are left out and "
                  "the fish they carried are not carried over", flush=True)
        woa = m.npzd.O2_ref is not None
        if plankton:                                        # with O2_init = "woa" the atlas oxygen is kept
            O2 = m.npzd.C[O2_IDX].copy()
            m.npzd.C = old["npzd"].C.copy()
            if woa:
                m.npzd.C[O2_IDX] = O2
        if m.verbose:
            miss = [n for n in m.sp.names if n not in names]
            print(f"warm start from {path} (day {day:g}): {len(pairs)} species, {len(m.ag):,d} agents"
                  + ((", nutrients and plankton too (oxygen from the World Ocean Atlas)" if woa else
                      ", plankton, nutrients and oxygen too") if plankton else "")
                  + (f"; new species start fresh: {', '.join(miss)}" if miss else ""), flush=True)

    def tally_selection(m, sp, killed, causes):
        """Accumulate the breeding values of the fish killed this step, by species and cause."""
        A, nC = m.ag.A[:len(sp)], len(CAUSES)
        for d, c in zip(killed, causes):
            key = sp * nC + CAUSES.index(c)
            m.sel_sum += np.stack([np.bincount(key, d * A[:, t], m.sp.n * nC) for t in range(len(TRAITS))],
                                  -1).reshape(m.sel_sum.shape)
            m.sel_n += np.bincount(key, d, m.sp.n * nC).reshape(m.sel_n.shape)

    def adaptation(m):
        """Population genetic diagnostics: realised G, F_ST between regions, and the selection differential
        (mean breeding value of the fish that died this step, minus the population mean, per cause)."""
        S, a, g = m.sp, m.ag, m.g
        out = dict(fst=np.full(S.n, np.nan), g_var=np.full((S.n, len(TRAITS)), np.nan),
                   sel=np.zeros((S.n, len(CAUSES), len(TRAITS))))
        if not len(a):
            return out
        nb = m.hy.get("fst_bands", 4)
        band = np.clip(((a.lat - g.lat_e[0]) / (g.lat_e[-1] - g.lat_e[0]) * nb).astype(int), 0, nb - 1)
        for q in range(S.n):
            w = a.sp == q
            if w.sum() < 2 * len(TRAITS):
                continue
            out["g_var"][q] = a.A[w].var(0)
            if m.gen.L:
                out["fst"][q] = m.gen.fst(a.gl[w], band[w], nb)
        dead = m.sel_sum / np.maximum(m.sel_n, 1e-300)[..., None]
        pop = np.stack([np.average(a.A[a.sp == q], 0, weights=a.n[a.sp == q]) if (a.sp == q).sum() > 1
                        else np.zeros(len(TRAITS)) for q in range(S.n)])
        out["sel"] = np.where(m.sel_n[..., None] > 0, dead - pop[:, None], 0.0)
        return out

    def reset_step_tallies(m):
        n = m.sp.n
        m.loss, m.catch, m.eggs = np.zeros((n, len(CAUSES))), np.zeros(n), np.zeros(n)
        m.mi_sum, m.mi_n = np.zeros(n), np.zeros(n)
        m.sel_sum = np.zeros((n, len(CAUSES), len(TRAITS)))
        m.sel_n = np.zeros((n, len(CAUSES)))
        m.born, m.died = 0, 0
        m.toN, m.toD = np.zeros(m.g.nc3), np.zeros(m.g.nc3)

    # ------------------------------------------------------------- utilities
    def force(m):
        g = m.g
        tf = m.t if m.spin is None else (m.t + m.spin[0]) % m.spin[1]     # spin-up replays the first days of forcing
        m.f = f = m.ocean(tf, m.dt)
        m.light = m.npzd.light(f["par"])
        m.uc = 0.5 * (f["U"][..., :-1] + f["U"][..., 1:]) / g.dy
        m.vc = 0.5 * (f["V"][:, :-1] / g.lx[:-1, None] + f["V"][:, 1:] / g.lx[1:, None])
        m.wc = f["W"]

    def cols(m):
        j, i = m.g.cell(m.ag.lon, m.ag.lat)
        return j, i, (m.ag.k * m.g.ny + j) * m.g.nx + i

    def fish_mass(m):
        a = m.ag
        return m.EB.sum() + m.EG.sum() + (a.n * (a.W + a.E + a.G)).sum()

    def total_N(m):
        return (m.npzd.C[:N_NUTR] * m.g.vol).sum() + RHO_N * m.fish_mass()

    def zbar(m):
        S = m.sp
        mean = m.ES[:, :, :, 0] / np.maximum(m.EN, 1e-300)[:, :, None]
        Z = S.mu[:, None, :, None, None] + S.sd[:, None, :, None, None] * mean
        return np.maximum(Z, TMIN[None, None, :, None, None])

    def newborn_genes(m, A, s, k=1):
        """Columns describing a new fish's genotype and phenotype for a target breeding value."""
        out = dict(A=A, Z=m.phenotype(s, A))
        if m.gen.L:
            out["gl"] = m.gen.match(A, s, m.rng, k)
            out["A"] = A = m.gen.breeding(out["gl"], s)
            out["Z"] = m.phenotype(s, A)
        return out

    def phenotype(m, s, A):
        S = m.sp
        e = m.rng.standard_normal(A.shape) * np.sqrt(1 - S.h2[s])
        return np.maximum(S.mu[s] + S.sd[s] * (A + e), TMIN)

    def migrating(m):
        """Homing is a pre-spawning migration: active from `migration_lead_days` before the season to its end."""
        win = m.sp.spawn_doy - np.array([m.be.get("migration_lead_days", 30), 0])
        return in_season(m.f["doy"], win % 365).astype(float)

    def choose(m, p):
        idx = (np.cumsum(p, 1) < m.rng.random(len(p))[:, None]).sum(1)
        return np.minimum(idx, p.shape[1] - 1)

    def drive_probs(m, here, w, hunger, ripe):
        """Probabilistic engine: P(behaviour) proportional to |drive weight| x urgency from the local stimuli.
        here: cues where the fish is (..., 6); w: heritable drive weights (..., 6); hunger, ripe in [0, 1]."""
        up = w >= 0
        urg = np.stack([hunger, 1 - here[..., 0], here[..., 2],
                        np.where(up[..., 3], 1 - here[..., 3], here[..., 3]),
                        np.where(up[..., 4], 1 - here[..., 4], here[..., 4]), ripe,
                        np.where(up[..., 6], 1 - here[..., 6], here[..., 6])], -1)
        wt = np.abs(w[..., CUE]) * np.clip(urg, 0, 1)
        wt = np.concatenate([wt, np.full(wt.shape[:-1] + (1,), m.be.get("random_drive", 0.2))], -1)
        return wt / np.maximum(wt.sum(-1, keepdims=True), 1e-12)

    def targets(m, F, w, valid, vertical):
        """Option each behaviour heads for: most food, best temperature, least risk, most (or least) company,
        darkest (or brightest) depth, spawning ground. F: cues at options (..., nopt, 6). Returns (..., 6).
        Hiding is vertical-only and homing horizontal-only; callers treat the other axis as staying put."""
        sg = np.where(w >= 0, 1.0, -1.0)[..., CUE]
        S = F[..., CUE] * sg[..., None, :]
        if not vertical:
            S[..., HIDE] = np.arange(F.shape[-2]) == 0
        else:
            S[..., SHELF] = -np.inf                                     # shelf is a horizontal choice only
        return np.argmax(np.where(valid[..., None], S, -np.inf), -2)

    def option_probs(m, P, tgt, valid, stay=None):
        """Expected choice over options implied by P(behaviour) and targets (the mean field of agent sampling).
        stay: current distribution over options for the behaviour that keeps position on this axis."""
        p = (P[..., :RANDOM, None] * (tgt[..., None] == np.arange(valid.shape[-1]))).sum(-2)
        if stay is not None:
            p += P[..., HOME, None] * (stay - (tgt[..., HOME, None] == np.arange(valid.shape[-1])))
        return p + P[..., RANDOM:] * valid / np.maximum(valid.sum(-1, keepdims=True), 1e-300)

    def ripe_euler(m):
        return np.clip(m.EG / np.maximum(m.sp.gonad_frac[:, None, None] * (m.EB - m.ER)[:, ADULT], 1e-300), 0, 1) \
            * m.migrating()[:, None, None]

    def fractions(m, UE, UW, VN, VS, Kh, extra=None, faces=None):
        """Outflow fractions (E, W, N, S) from upwind advection, diffusion and optional swimming.
        faces: open-face mask per direction (3-D wet faces for tracers, 2-D columns for depth-integrated fish)."""
        g, dt = m.g, m.dt_s
        faces = g.nvalid[1:] if faces is None else faces
        kd = [Kh * dt / g.dx ** 2] * 2 + [Kh * dt * g.lx[1:, None] / (g.dy * g.area),
                                          Kh * dt * g.lx[:-1, None] / (g.dy * g.area)]
        adv = [np.maximum(UE, 0), np.maximum(-UW, 0), np.maximum(VN, 0), np.maximum(-VS, 0)]
        f = [(a * dt / g.area + d) * v for v, (a, d) in zip(faces, zip(adv, kd))]
        if extra is not None:
            f = [fi + e for fi, e in zip(f, extra)]
        n = max(1, int(np.ceil(np.max(sum(f)) / 0.5)))
        return [fi / n for fi in f], n

    def transport3d(m, M, Qx, Qy, Kh):
        """Advection and diffusion of a 3-D tracer mass over topography. Face fluxes scale with the open
        thickness and the leftover column divergence becomes a vertical velocity (zero through the sea bed),
        so a uniform concentration stays uniform even with partial bottom cells."""
        g, dt = m.g, m.dt_s
        adv = [np.maximum(Qx[..., 1:], 0), np.maximum(-Qx[..., :-1], 0),
               np.maximum(Qy[:, 1:], 0), np.maximum(-Qy[:, :-1], 0)]
        dif = [Kh * g.tx[..., 1:] * g.dy / g.dx, Kh * g.tx[..., :-1] * g.dy / g.dx,          # m3/s
               Kh * g.ty[:, 1:] * g.lx[1:, None] / g.dy, Kh * g.ty[:, :-1] * g.lx[:-1, None] / g.dy]
        f = [(a + d) * dt / g.volx * v for v, a, d in zip(g.wetf, adv, dif)]
        n = max(1, int(np.ceil(np.max(sum(f)) / 0.5)))
        f = [x / n for x in f]
        for _ in range(n):
            M = exchange(M, f, g.periodic)
        return M

    def transport(m, M, fr, expand=0):
        f, n = fr
        f = [fi.reshape(fi.shape[:-2] + (1,) * expand + fi.shape[-2:]) for fi in f]
        for _ in range(n):
            M = exchange(M, f, m.g.periodic)
        return M

    # ------------------------------------------------------------- behaviour
    def features(m, s, st, k, j, i, Topt, Tw, vertical, home=0.0):
        """[habitat suitability (thermal x oxygen), food, predation risk, crowding, light, homing, shelf]
        at candidate locations. Incorporates metabolic index Phi and hypoxia avoidance."""
        T = m.f["T"][k, j, i]
        O2 = m.npzd.C[O2_IDX, k, j, i]
        Wref = m.sp.weight(m.sp.Lref[s, st], s)
        mi = metabolic_index(O2, T, Wref, m.sp.A_o[s], m.sp.E_o[s], m.sp.eps_o[s])
        o2_suit = np.clip((mi - 1.0) / np.maximum(m.sp.phi_crit[s] - 1.0, 1e-6), 0.0, 1.0)
        therm_suit = np.exp(-0.5 * ((T - Topt) / Tw) ** 2)
        hab_suit = therm_suit * o2_suit
        return np.stack(np.broadcast_arrays(hab_suit, m.food[s, st, k, j, i],
                                            m.risk[s, st, k, j, i], m.crowd[s, k, j, i],
                                            m.lightf[k, j, i] * vertical, home, m.shelf[j, i]), -1)

    def fields(m):
        """Food, risk and crowding fields seen by each species/stage (from the current fish+plankton state)."""
        S, g, a = m.sp, m.g, m.ag
        nf, shp = 4 * S.n, (S.n, 4, g.nz, g.ny, g.nx)
        adult = (np.arange(4) == ADULT)[None, :, None, None]
        fish = ((m.EB + m.EG[:, None] * adult)[:, :, None] * m.Fv).reshape(nf, g.nz, g.ny, g.nx)
        if len(a):
            j, i, cell = m.cols()
            fish = fish + np.bincount((a.sp * 4 + a.st) * g.nc3 + cell, a.n * (a.W + a.E + a.G),
                                      nf * g.nc3).reshape(fish.shape)
        dens = np.concatenate([m.npzd.C[1:1+N_PLANK] / RHO_N, fish / g.volx])
        K = S.K_food[:, None, None, None, None]
        A = np.einsum("xt,tkji->xkji", S.PREF.reshape(nf, -1), dens).reshape(shp)
        r = np.einsum("xq,qkji->xkji", S.RISK, dens[N_PLANK:]).reshape(shp)
        c = dens[N_PLANK:].reshape(shp)[:, JUV:].sum(1)
        m.food, m.risk = A / (K + A), r / (K + r)
        m.crowd = c / (K[:, 0] + c)
        m.lightf = m.light / (m.light + 10)

    def vertical(m):
        """Vertical distribution of every biomass class (softmax of the mean innate policy)."""
        S, g = m.sp, m.g
        m.Zb = Zb = m.zbar()
        s, st = np.arange(S.n)[:, None, None, None, None], np.arange(4)[None, :, None, None, None]
        k = np.arange(g.nz)[None, None, :, None, None]
        j, i = np.mgrid[:g.ny, :g.nx]
        F = m.features(s, st, k, j, i, Zb[:, :, None, T_OPT], Zb[:, :, None, T_WID], True)
        w = np.moveaxis(Zb[:, :, BEH], 2, -1)
        if m.be.get("engine", "softmax") == "probabilistic":
            here = (F * m.Fv[..., None]).sum(2)
            ripe = m.ripe_euler()[:, None] * (np.arange(4) == ADULT)[None, :, None, None]
            P = m.drive_probs(here, w, 1 - here[..., 1], ripe)
            valid = np.moveaxis(m.LEV, 2, -1)
            tgt = m.targets(np.moveaxis(F, 2, -2), w, valid, True)
            Fv = np.moveaxis(m.option_probs(P, tgt, valid, np.moveaxis(m.Fv, 2, -1)), -1, 2) * m.LEV
            tot = Fv.sum(2, keepdims=True)             # every fish must sit somewhere it is allowed to be
            m.Fv = np.where(tot > 0, Fv / np.where(tot > 0, tot, 1.0), m.LEV / np.maximum(m.LEV.sum(2, keepdims=True), 1))
            return
        logit = np.where(m.LEV, m.be["beta"] * (F * w[:, :, None]).sum(-1), -np.inf)
        m.Fv = np.where(m.LEV.any(2, keepdims=True), softmax(logit, 2), 0.0)

    def move_euler(m):
        S, g, f, Zb, Fv = m.sp, m.g, m.f, m.Zb, m.Fv
        w = lambda X: np.einsum("sckji,kji->scji", Fv, X)
        UE, UW, VN, VS = w(f["U"][..., 1:]), w(f["U"][..., :-1]), w(f["V"][:, 1:]), w(f["V"][:, :-1])
        s, st = np.arange(S.n)[:, None, None, None, None, None], np.array([JUV, ADULT])[None, :, None, None, None, None]
        k = np.arange(g.nz)[None, None, None, :, None, None]
        Za = Zb[:, JUV:]
        ripe = m.ripe_euler()
        gd = np.where(m.has_ground[:, None, None, None], m.gdir, m.spawn_attraction())
        home = np.einsum("od,sdji->soji", DIRS, gd) * ripe[:, None]
        H = np.stack([np.zeros_like(home), home], 1)[:, :, :, None]
        F = m.features(s, st, k, g.nj[:, None], g.ni[:, None], Za[:, :, None, None, T_OPT],
                       Za[:, :, None, None, T_WID], False, H)
        Fh = (F * Fv[:, JUV:, None, ..., None]).sum(3)
        valid = g.nvalid[None] & m.HOME[:, None, None]    # stay inside the species' range
        valid = np.broadcast_to(valid, (S.n, 1, 5, g.ny, g.nx)).copy()
        valid[:, :, 0] = m.HOME[:, None]                 # staying put is always allowed in range
        w = np.moveaxis(Za[:, :, BEH], 2, -1)
        if m.be.get("engine", "softmax") == "probabilistic":
            here = Fh[:, :, 0]
            P = m.drive_probs(here, w, 1 - here[..., 1], ripe[:, None] * np.array([0, 1])[None, :, None, None])
            vo = np.moveaxis(valid, 2, -1)
            p = np.moveaxis(m.option_probs(P, m.targets(np.moveaxis(Fh, 2, -2), w, vo, False), vo), -1, 2)
        else:
            p = softmax(np.where(valid, m.be["beta"] * (Fh * w[:, :, None]).sum(-1), -np.inf), 2)
        L = S.length((m.EB - m.ER)[:, JUV:] / np.maximum(m.EN[:, JUV:], 1e-300))
        swim = Za[:, :, SPEED] * L / 100 * m.dt_s
        extra = np.zeros((4, S.n, 4, g.ny, g.nx))
        extra[:, :, JUV:] = np.moveaxis(p[:, :, 1:] * np.minimum(1, swim[:, :, None] * g.flen / g.area), 2, 0)
        fr = m.fractions(UE, UW, VN, VS, m.be.get("fish_kh", 500.0), list(extra))
        dest = [np.roll(m.HOME, -1, -1), np.roll(m.HOME, 1, -1),                 # E, W, N, S destination in range
                np.roll(m.HOME, -1, -2), np.roll(m.HOME, 1, -2)]
        if not g.periodic:
            dest[0][..., -1] = dest[1][..., 0] = False
        dest[2][:, -1] = dest[3][:, 0] = False
        fr = ([f * d[:, None] for f, d in zip(fr[0], dest)], fr[1])              # no drift out of the home range
        m.EN, m.EB, m.ER = m.transport(m.EN, fr), m.transport(m.EB, fr), m.transport(m.ER, fr)
        m.ES = m.transport(m.ES, fr, expand=2)
        m.EGEN = m.transport(m.EGEN, fr)
        m.NJ = m.transport(m.NJ, ([f[:, JUV, None] for f in fr[0]], fr[1]))
        m.EG = m.transport(m.EG, ([fi[:, ADULT] for fi in fr[0]], fr[1]))

    def move_agents(m):
        a, S, g, be = m.ag, m.sp, m.g, m.be
        if not len(a):
            return
        j, i, _ = m.cols()
        s, st, ar = a.sp, a.st, np.arange(len(a))
        th = a.Z[:, BEH] + a.dl
        Topt, Tw = a.Z[:, T_OPT, None], a.Z[:, T_WID, None]
        a.r0 = np.log(np.maximum(a.n * (a.W + a.E + a.G), 1e-300))
        act, swimmer, prob = st >= LARVA, st >= JUV, be.get("engine", "softmax") == "probabilistic"
        ripe = np.clip(a.G / (S.gonad_frac[s] * a.W), 0, 1) * (st == ADULT) * m.migrating()[s]
        natal = np.stack([(a.lon0 - a.lon) * np.cos(np.radians(a.lat)), a.lat0 - a.lat])   # back to where it was born
        natal = natal / np.maximum(np.hypot(*natal), m.be.get("natal_tol_deg", 2.0))
        gd = np.where(m.has_ground[s][None], m.gdir[s, :, j, i].T, natal)     # (2, nagent)
        home = ripe[:, None] * (gd.T @ DIRS.T)
        valid = g.nvalid[:, j, i].T & m.HOME[s[:, None], g.nj[:, j, i].T, g.ni[:, j, i].T]
        valid[:, 0] = True
        F = m.features(s[:, None], st[:, None], np.arange(g.nz)[None], j[:, None], i[:, None], Topt, Tw, True)
        if prob:                                # one behaviour per decision, held for 1..inertia_steps steps
            P = m.drive_probs(F[ar, a.k], th, np.clip(1 - a.E / (S.reserve_max[s] * a.W), 0, 1), ripe)
            new = a.hold <= 0
            a.beh = np.where(act, np.where(new, m.choose(P), a.beh), -1)
            a.hold = np.where(new, m.rng.integers(0, be.get("inertia_steps", 1), len(a)), a.hold - 1)
            b = np.maximum(a.beh, 0)
            lev = m.LEV[s, st, :, j, i]
            tgt = m.targets(F, th, lev, True)
            k = np.select([b == HOME, b == RANDOM], [a.k, m.choose(lev / np.maximum(lev.sum(1, keepdims=True), 1e-300))],
                          tgt[ar, np.minimum(b, HOME)])
            a.k = np.where(act, k, a.k)
            F = m.features(s[:, None], st[:, None], a.k[:, None], g.nj[:, j, i].T, g.ni[:, j, i].T, Topt, Tw, False, home)
            tgt = m.targets(F, th, valid, False)
            o = np.where(b == RANDOM, m.choose(valid / valid.sum(1, keepdims=True)), tgt[ar, np.minimum(b, HOME)])
            o = np.where(swimmer, o, 0)
            a.grad = np.zeros_like(a.grad)
        else:
            p = softmax(np.where(m.LEV[s, st, :, j, i], be["beta"] * (F * th[:, None]).sum(-1), -np.inf), 1)
            k = m.choose(p)
            a.k = np.where(act, k, a.k)
            grad = act[:, None] * (F[ar, k] - (p[..., None] * F).sum(1))
            F = m.features(s[:, None], st[:, None], a.k[:, None], g.nj[:, j, i].T, g.ni[:, j, i].T, Topt, Tw, False, home)
            p = softmax(np.where(valid, be["beta"] * (F * th[:, None]).sum(-1), -np.inf), 1)
            o = np.where(swimmer, m.choose(p), 0)
            a.grad = grad + swimmer[:, None] * (F[ar, o] - (p[..., None] * F).sum(1))
            a.beh = np.full(len(a), -1)
        dist = a.Z[:, SPEED] * S.length(a.W, s) / 100 * m.dt_s * swimmer
        noise = m.rng.standard_normal((2, len(a))) * np.sqrt(2 * be.get("fish_kh", 500.0) * m.dt_s)
        de = DIRS[o, 0] * dist + m.uc[a.k, j, i] * m.dt_s + noise[0]
        dn = DIRS[o, 1] * dist + m.vc[a.k, j, i] * m.dt_s + noise[1]
        lon = a.lon + np.degrees(de / (R_EARTH * np.cos(np.radians(a.lat))))
        lat = a.lat + np.degrees(dn / R_EARTH)
        if g.periodic:
            lon = g.lon_e[0] + (lon - g.lon_e[0]) % (g.lon_e[-1] - g.lon_e[0])
        ok = g.inside(*g.cell(lon, lat))
        jt, it = g.cell(np.where(ok, lon, a.lon), np.where(ok, lat, a.lat))
        ok &= m.HOME[s, np.clip(jt, 0, g.ny - 1), np.clip(it, 0, g.nx - 1)]       # drift stops at the range edge
        jn, iN = g.cell(np.where(ok, lon, a.lon), np.where(ok, lat, a.lat))
        a.k = np.minimum(a.k, g.kbot[np.clip(jn, 0, g.ny - 1), np.clip(iN, 0, g.nx - 1)])
        a.lon, a.lat = np.where(ok, lon, a.lon), np.where(ok, lat, a.lat)

    # ------------------------------------------------ aggregation / disaggregation
    def to_euler(m, s, st, j, i, n, M, G, A, E, gen=0):
        """Add fish (n per row; soma+reserve M, reserve E and gonad G per fish; standardised breeding values A;
        generation gen)."""
        g, sh = m.g, m.EN.shape
        idx = np.ravel_multi_index((s, st, j, i), sh)
        m.EN += np.bincount(idx, n, m.EN.size).reshape(sh)
        m.EGEN += np.bincount(idx, n * np.broadcast_to(gen, np.shape(n)), m.EN.size).reshape(sh)
        juv = st == JUV
        if juv.any():                                          # juveniles join the bin matching their progress
            S, K = m.sp, m.NJ.shape[1]
            Wj, Wm = S.weight(S.L_juv[s[juv]], s[juv]), S.weight(S.mu[s[juv], L_MAT], s[juv])
            W = np.broadcast_to(M - E, np.shape(n))[juv]
            b = np.clip(((W - Wj) / np.maximum(Wm - Wj, 1e-12) * K).astype(int), 0, K - 1)
            np.add.at(m.NJ, (s[juv], b, j[juv], i[juv]), np.broadcast_to(n, np.shape(st))[juv])
        m.EB += np.bincount(idx, n * M, m.EN.size).reshape(sh)
        m.ER += np.bincount(idx, n * E, m.EN.size).reshape(sh)
        m.EG += np.bincount(np.ravel_multi_index((s, j, i), m.EG.shape), n * G, m.EG.size).reshape(m.EG.shape)
        ntr = A.shape[1]
        ix = np.broadcast_arrays(s[:, None, None], st[:, None, None], np.arange(ntr)[None, :, None],
                                 np.arange(4)[None, None], j[:, None, None], i[:, None, None])
        P = n[:, None, None] * A[:, :, None] ** np.arange(1, 5)
        m.ES += np.bincount(np.ravel_multi_index(ix, m.ES.shape).ravel(), P.ravel(), m.ES.size).reshape(m.ES.shape)

    def remove(m, mask, cause):
        """Drop agents, recording their life histories."""
        a = m.ag
        if not mask.any():
            return
        S = m.sp
        idx = np.nonzero(mask)[0]
        rec = {k: a.d[k][idx] for k in ("id", "sp", "origin", "gen", "mother", "father", "tb", "lon0", "lat0",
                                        "lon", "lat", "st", "n")}
        rec.update(t_end=np.full(len(idx), m.t), cause=np.broadcast_to(cause, mask.shape)[idx].astype(np.int8),
                   L=S.length(a.W[idx], a.sp[idx]), A=S.mu[a.sp[idx]] + S.sd[a.sp[idx]] * a.A[idx])
        m.obits.append(rec)
        m.died += len(idx)
        a.keep(~mask)

    def merge(m, mask):
        a = m.ag
        if mask.any():
            j, i, _ = m.cols()
            m.to_euler(a.sp[mask], a.st[mask], j[mask], i[mask], a.n[mask], (a.W + a.E)[mask], a.G[mask], a.A[mask], a.E[mask],
                        a.gen[mask])
            m.remove(mask, CAUSES.index("merged"))

    def sync(m):
        """Agents outside IBM cells or depleted below the minimum super-individual size become biomass."""
        a = m.ag
        if len(a):
            j, i, _ = m.cols()
            small = a.n * (a.W + a.E + a.G) < np.nan_to_num(m.hy.get("min_agent_frac", 1e-4) * m.bref[a.sp], nan=-1)
            m.merge(~m.ibm[j, i] | small)

    def disaggregate(m, colmask):
        S, g, a, rng = m.sp, m.g, m.ag, m.rng
        k = m.hy.get("agents_per_class", 8)
        s, st, j, i = np.nonzero(m.EN * colmask >= k)
        N = m.EN[s, st, j, i]
        M, E = m.EB[s, st, j, i] / N, m.ER[s, st, j, i] / N
        G = np.where(st == ADULT, m.EG[s, j, i] / np.maximum(m.EN[s, ADULT, j, i], 1e-300), 0)
        ne = np.floor(N / k)
        ok = ~(ne * (M + G) < m.hy.get("min_agent_frac", 1e-4) * m.bref[s] * 3)
        room = max(0, m.hy.get("max_agents", 10 ** 5) - len(a)) // k
        
        if m.hy.get("budget_by_biomass", True):
            # Share the agent budget in proportion to each species' biomass, so the individuals sample where
            # the fish actually are, and fill each species' quota with its largest classes first.
            B = m.EB.sum((1, 2, 3)) + (np.bincount(a.sp, a.n * (a.W + a.E + a.G), m.sp.n) if len(a) else 0)
            share = np.maximum(B / max(B.sum(), 1e-300), m.hy.get("min_species_share", 0.05))
            cap = share / share.sum() * m.hy.get("max_agents", 10 ** 5) / k
            have = np.bincount(a.sp, minlength=m.sp.n) / k if len(a) else np.zeros(m.sp.n)
            big = np.argsort(-ne * (M + G))                      # biggest classes first
            rank = np.zeros(len(s))
            for q in range(m.sp.n):
                w = ok[big] & (s[big] == q)
                rank[big[w]] = np.cumsum(w)[w] + have[q]
            ok &= rank <= cap[s]
        else:
            ok &= np.cumsum(ok) <= room
            
        s, st, j, i, N, M, E, G, ne = [x[ok] for x in (s, st, j, i, N, M, E, G, ne)]
        if not len(s):
            return
        A = sample_moments(m.ES[s, st, :, :, j, i], N, k, rng)
        frac = ne * k / N
        gm = m.EGEN[s, st, j, i] / N                          # the class's mean generation
        m.EGEN[s, st, j, i] -= ne * k * gm
        juv = st == JUV
        m.NJ[s[juv], :, j[juv], i[juv]] *= (1 - (ne * k / N)[juv])[:, None]
        m.EN[s, st, j, i] -= ne * k
        m.EB[s, st, j, i] -= ne * k * M
        m.ER[s, st, j, i] -= ne * k * E
        m.ES[s, st, :, :, j, i] *= (1 - frac)[:, None, None]
        ad = st == ADULT
        m.EG[s[ad], j[ad], i[ad]] -= (ne * k * G)[ad]
        r = np.repeat(np.arange(len(s)), k)
        sr, str_ = s[r], st[r]
        W = M[r] - E[r]
        u = rng.random(len(r))
        gr = gm[r]
        gen = np.floor(gr) + (rng.random(len(r)) < gr - np.floor(gr))   # integer, with the class mean preserved
        age = np.select([str_ == EGG, str_ == LARVA, str_ == JUV],
                        [0 * u, S.egg_days[sr] * (1 + 5 * u), 30 + 300 * u], S.lifespan_days[sr] * (0.15 + 0.35 * u))
        lon = g.lon_e[i[r]] + rng.random(len(r)) * g.dlon
        lat = g.lat_e[j[r]] + rng.random(len(r)) * g.dlat
        for q in np.unique(sr):
            if np.isnan(m.bref[q]):
                m.bref[q] = np.mean((ne * (M + G))[s == q])
        a.add(sp=sr, st=str_, sex=np.arange(len(r)) % 2, k=m.choose(m.Fv[s, st, :, j, i][r]), n=ne[r], W=W,
              E=M[r] - W, G=G[r], age=age, dev=np.where(str_ == EGG, 0.5 * u, 0),
              amax=np.maximum(S.lifespan_days[sr] * np.exp(0.15 * rng.standard_normal(len(r))), age + 30),
              tb=m.t - age, lon=lon, lat=lat, lon0=lon, lat0=lat, origin=1, mother=-1, father=-1, gen=gen,
              **m.newborn_genes(A.reshape(-1, A.shape[-1]), sr, k))
        m.born += len(r)

    def regulate(m):
        """Choose IBM cells (static regions, dynamic fronts/blooms, all or none) and re-balance representations."""
        g, hy = m.g, m.hy
        new = np.full((g.ny, g.nx), hy["mode"] == "all")
        if hy["mode"] in ("static", "dynamic"):
            for r in hy.get("regions", []):
                new |= haversine(g.lon[None], g.lat[:, None], r["lon"], r["lat"]) <= r["radius_km"]
        if hy["mode"] == "dynamic":
            if not hasattr(m, "_near"):                      # nearest ocean cell, for filling land once
                from scipy.ndimage import distance_transform_edt
                m._near = tuple(distance_transform_edt(~g.mask, return_distances=False, return_indices=True))
            score = 0
            for x in (m.f["T"][0], np.log(np.maximum(m.npzd.C[1, 0], 1e-6))):
                gx = np.hypot(*np.gradient(x[m._near]))      # land filled from the coast: no artificial step
                score = score + gx / max(gx[g.mask].std(), 1e-12)
            if hy.get("front_species_weight", 0.0):          # follow where each species aggregates, not just
                fsc = 0                                       # whichever species has the most total biomass
                for q in range(m.sp.n):
                    fish = np.log(np.maximum(m.EB[q].sum(0) + 1e-6, 1e-6))
                    gx = np.hypot(*np.gradient(fish[m._near]))
                    fsc = fsc + gx / max(gx[g.mask].std(), 1e-12)      # normalised per species, then combined
                score = score + hy["front_species_weight"] * fsc / m.sp.n
            frac = hy.get("dynamic_fraction", 0.1)
            rank = np.zeros(g.mask.shape)
            rank[g.mask] = np.argsort(np.argsort(-score[g.mask] + 1e-9 * m.rng.random(g.mask.sum()))) / g.mask.sum()
            new |= g.mask & ((rank < frac) | m.ibm & (rank < 1.5 * frac))
        m.ibm = new & g.mask
        m.sync()
        m.disaggregate(m.ibm)
        m.next_reg = m.t + hy.get("update_days", 10)

    # ------------------------------------------------ feeding, growth and mortality
    def units(m):
        """Flatten agents and biomass (class, column, level) groups into one set of interacting units."""
        S, g, a = m.sp, m.g, m.ag
        Nk = m.EN[:, :, None] * m.Fv
        # Classes holding a negligible share of their species are skipped here; grow_euler carries them over
        # untouched, so mass is preserved exactly and the interacting set stays a manageable size.
        keep = m.EB > m.hy.get("min_class_frac", 1e-7) * m.EB.sum((1, 2, 3))[:, None, None, None]
        s, st, k, j, i = np.nonzero((Nk > 0) & keep[:, :, None])
        M, E = m.EB[s, st, j, i] / m.EN[s, st, j, i], m.ER[s, st, j, i] / m.EN[s, st, j, i]
        G = np.where(st == ADULT, m.EG[s, j, i] / np.maximum(m.EN[s, ADULT, j, i], 1e-300), 0)
        ja, ia, _ = m.cols()
        cat = lambda x, y: np.concatenate([x, y])
        u = dict(s=cat(a.sp, s), st=cat(a.st, st), k=cat(a.k, k), j=cat(ja, j), i=cat(ia, i), n=cat(a.n, Nk[s, st, k, j, i]),
                 M=cat(a.W + a.E, M), G=cat(a.G, G), W=cat(a.W, M - E), E=cat(a.E, E),
                 Z=np.concatenate([a.Z, m.Zb[s, st, :, j, i]]), na=len(a),
                 ecc=np.ravel_multi_index((s, st, j, i), m.EN.shape))
        u["L"] = S.length(u["W"], u["s"])
        u["cell"] = (u["k"] * g.ny + u["j"]) * g.nx + u["i"]
        return u

    def feed_and_die(m, u):
        """Size- and taxon-structured multi-prey feeding, respiration and mortality for every unit."""
        S, g, dt, rng = m.sp, m.g, m.dt, m.rng
        s, st, n, M, G, Z, cell, na, nt = u["s"], u["st"], u["n"], u["M"], u["G"], u["Z"], u["cell"], u["na"], S.ntype
        T = m.f["T"].ravel()[cell]
        phi = np.exp(-0.5 * ((T - Z[:, T_OPT]) / Z[:, T_WID]) ** 2)
        fQ = Q10 ** ((T - S.T_ref[s]) / 10)
        mi = metabolic_index(m.npzd.C[O2_IDX].ravel()[cell], T, M, S.A_o[s], S.E_o[s], S.eps_o[s])
        aer = np.clip((mi - 1) / np.maximum(S.phi_crit[s] - 1, 1e-6), 0, 1)      # aerobic scope for activity
        Cmax = S.cmax[s] * Z[:, METAB] * M ** (2 / 3) * fQ * phi * aer * (st > EGG)
        R = S.rmet[s] * Z[:, METAB] * M ** 0.8 * fQ * (ACTIVITY[st] + (st >= JUV) * S.c_swim[s] * (Z[:, SPEED] / S.mu[s, SPEED]) ** 2)
        ft, mass = N_PLANK + 4 * s + st, n * (M + G)
        flat, size = cell * nt + ft, g.nc3 * nt
        B = np.bincount(flat, mass, size).reshape(g.nc3, nt)
        B[:, :N_PLANK] = (m.npzd.C[1:1+N_PLANK] * g.vol).reshape(N_PLANK, -1).T / RHO_N
        food = m.cfg["npzd"].get("detritus_food", "pool")
        if food == "bed_flux":
            # detritus is the benthos proxy: only on the sea bed, and only what sinks onto it this step
            flux = (m.npzd.C[4] * m.npzd.c["w_sink"] * dt * g.area[None]).ravel() / RHO_N    # mmol N -> g fish
            B[:, 3] = np.where(g.bed.ravel(), np.minimum(B[:, 3], flux), 0.0)
        elif food == "bed":                                   # the benthos proxy lives on the sea bed only
            B[:, 3] = np.where(g.bed.ravel(), B[:, 3], 0.0)
        Lp = np.bincount(flat, mass * u["L"], size).reshape(g.nc3, nt) / np.maximum(B, 1e-300)
        Lp[:, :N_PLANK] = PLANKTON_L
        Lp[B <= 0] = 1.0
        p = np.nonzero(Cmax > 0)[0]
        cp, sp = cell[p], s[p]
        avail = S.taxo[sp, st[p]] * S.size_window(Lp[cp] / u["L"][p, None], sp[:, None]) * B[cp]
        avail[:, :N_PLANK - 1] *= S.plankton_taper(u["L"][p], sp)[:, None]
        tot = avail.sum(1)
        A = tot / g.volx.ravel()[cp]
        demand = avail * (n[p] * Cmax[p] * A / (S.K_food[sp] + A) * dt / np.maximum(tot, 1e-300))[:, None]
        D = np.bincount((cp[:, None] * nt + np.arange(nt)).ravel(), demand.ravel(), size).reshape(g.nc3, nt)
        scale = np.where(D > 0, np.minimum(1, 0.9 * B / np.where(D > 0, D, 1)), 1)
        lam = np.where(B > 0, D * scale / np.maximum(B, 1e-300), 0)
        sel = 1 / (1 + np.exp(-np.log(19) * (u["L"] - S.L50[s]) / (S.L95[s] - S.L50[s])))
        Fs = m.fishing_rate(u["j"], u["i"], s)
        mF = Fs / 365 * sel * (m.t >= m.cfg["run"].get("fishing_start_day", 0))
        mT = S.m_thermal[s] * np.maximum(np.abs(T - Z[:, T_OPT]) / Z[:, T_WID] - 1, 0) ** 2
        young = (st == LARVA) | (st == JUV)
        col = u["j"] * g.nx + u["i"]
        crowd = np.bincount((s * 4 + st) * g.ncol + col, mass, 4 * S.n * g.ncol)[(s * 4 + st) * g.ncol + col] / g.area.ravel()[col]
        fade = np.exp(-u["L"] / S.crowd_L[s])                          # crowding matters less as fish grow
        mN = (S.m0[s] * M ** -0.25 + S.m_early[s] * (st <= LARVA)) * (1 + young * crowd * fade / S.K_nursery[s]
                                                                       + (st == ADULT) * crowd / S.K_adult[s])
        mH = S.m_hypoxia[s] * np.clip(1 - mi, 0, 1) ** 2                          # suffocation below phi = 1
        mtot = mF + mT + mN + mH
        pdie, lam_u = 1 - np.exp(-mtot * dt), lam[cell, ft]
        dP = n * lam_u
        dO = (n - dP) * pdie
        dF, dT, dH = dO * mF / mtot, dO * mT / mtot, dO * mH / mtot
        if na:
            b = lambda N, q: rng.binomial(np.round(N).astype(np.int64), np.clip(q, 0, 1)).astype(float)
            dP[:na] = b(n[:na], lam_u[:na])
            dO[:na] = b(n[:na] - dP[:na], pdie[:na])
            dF[:na] = b(dO[:na], mF[:na] / mtot[:na])
            dT[:na] = b(dO[:na] - dF[:na], mT[:na] / np.maximum(mT[:na] + mN[:na] + mH[:na], 1e-300))
            dH[:na] = b(dO[:na] - dF[:na] - dT[:na], mH[:na] / np.maximum(mH[:na] + mN[:na], 1e-300))
        dN = dO - dF - dT - dH
        realised = np.bincount(flat, dP * (M + G), size).reshape(g.nc3, nt)
        ratio = np.ones_like(B)
        ratio[:, N_PLANK:] = np.where(lam[:, N_PLANK:] > 0, realised[:, N_PLANK:] / np.maximum(lam[:, N_PLANK:] * B[:, N_PLANK:], 1e-300), 0)
        eaten = demand * (scale * ratio)[cp]
        for q in range(N_PLANK):
            m.npzd.C[1 + q] -= (np.bincount(cp, eaten[:, q], g.nc3) * RHO_N / g.volx.ravel()).reshape(g.nz, g.ny, g.nx)
        intake = np.zeros(len(n))
        intake[p] = eaten.sum(1)
        ns = n - dP - dO
        ci = intake / np.maximum(n, 1e-300)
        dead = M + G
        m.toD += np.bincount(cell, ns * (1 - S.alpha[s]) * ci + (n - ns) * ci + (dN + dT + dH) * dead, g.nc3)
        m.toN += np.bincount(cell, ns * R * dt, g.nc3)
        for c, d in (("predation", dP), ("fishing", dF), ("thermal", dT), ("hypoxia", dH), ("natural", dN)):
            m.loss[:, CAUSES.index(c)] += np.bincount(s, d * dead, S.n)
        m.catch += np.bincount(s, dF * dead, S.n)
        m.ext -= RHO_N * (dF * dead).sum()
        if na:                                     # who is dying, of what, and with which genes
            m.tally_selection(s[:na], [dP[:na], dF[:na], dT[:na], dH[:na], dN[:na]],
                              ("predation", "fishing", "thermal", "hypoxia", "natural"))
        m.mi_sum += np.bincount(s, n * mi, S.n)
        m.mi_n += np.bincount(s, n, S.n)
        return dict(ns=ns, net=S.alpha[s] * ci - R * dt, T=T, deaths=np.stack([dP, dF, dT, dH, dN], 1))

    def grow_agents(m, u, r):
        a, S, na = m.ag, m.sp, u["na"]
        if not na:
            return
        net, cell = r["net"][:na], u["cell"][:na]
        pos = net >= 0
        dE = np.where(pos, np.minimum(net, np.maximum(S.reserve_max[a.sp] * a.W - a.E, 0)), net)
        rest = np.where(pos, net - dE, 0)
        soma = np.where(a.st == ADULT, (1 - S.kappa_R[a.sp]) * np.clip(1 - a.W / S.W_inf[a.sp], 0, 1), 1) * rest
        a.E, a.G, a.W = a.E + dE, a.G + rest - soma, a.W + soma
        a.n = r["ns"][:na]
        starve = (a.E < 0) & (a.n > 0)
        m.toN += np.bincount(cell, starve * a.n * a.E, m.g.nc3)
        m.toD += np.bincount(cell, starve * a.n * (a.W + a.G), m.g.nc3)
        m.loss[:, CAUSES.index("starvation")] += np.bincount(a.sp, starve * a.n * (a.W + a.G), S.n)
        m.tally_selection(a.sp, [starve * a.n], ("starvation",))
        order = [CAUSES.index(c) for c in ("predation", "fishing", "thermal", "hypoxia", "natural")]
        cause = np.where(starve, CAUSES.index("starvation"), np.take(order, np.argmax(r["deaths"][:na], 1)))
        a.E = np.where(starve, 0, a.E)
        a.n = np.where(starve, 0, a.n)
        m.remove(a.n <= 0, cause)

    def grow_euler(m, u, r):
        """Same allocation as agents (reserve first, then structure or gonad). A class only knows its mean reserve,
        so individual reserves are taken as uniform within +/- reserve_spread of it and exactly the fraction
        whose reserve would go negative starves; survivors and casualties keep their conditional means."""
        S, g, na = m.sp, m.g, u["na"]
        e = slice(na, None)
        net, ns, G, s, st, cell, ecc = r["net"][e], r["ns"][e], u["G"][e], u["s"][e], u["st"][e], u["cell"][e], u["ecc"]
        W, E = u["W"][e], u["E"][e]
        net = (np.bincount(ecc, ns * net, m.EN.size) / np.maximum(np.bincount(ecc, ns, m.EN.size), 1e-300))[ecc]  # fish move between levels
        pos = net >= 0
        dE = np.where(pos, np.minimum(net, np.maximum(S.reserve_max[s] * W - E, 0)), net)
        rest = np.where(pos, net - dE, 0)
        soma = np.where(st == ADULT, (1 - S.kappa_R[s]) * np.clip(1 - W / S.W_inf[s], 0, 1), 1) * rest
        E1, d = E + dE, m.hy.get("reserve_spread", 0.5) * np.maximum(E, 0)
        lo, hi = E1 - d, E1 + d
        f = np.where(d > 0, np.clip(-lo / np.maximum(hi - lo, 1e-300), 0, 1), (E1 < 0).astype(float))
        e_dead = np.where(d > 0, 0.5 * (lo + np.minimum(hi, 0)), E1)          # mean reserve of the starving (< 0)
        e_live = np.where(f < 1, np.where(d > 0, 0.5 * (np.maximum(lo, 0) + hi), E1), 0)
        dst = ns * f
        m.toN += np.bincount(cell, dst * np.minimum(e_dead, 0), g.nc3)         # respiration they could not pay
        m.toD += np.bincount(cell, dst * (W + G), g.nc3)
        m.loss[:, CAUSES.index("starvation")] += np.bincount(s, dst * (W + G), S.n)
        ns = ns - dst
        Wn, Gn = W + soma, G + rest - soma
        sz, sh = m.EN.size, m.EN.shape
        EN = np.bincount(ecc, ns, sz).reshape(sh)
        orphan = (np.bincount(ecc, None, sz).reshape(sh) == 0) & (m.EN > 0)   # class with no feeding unit
        if orphan.any():                                  # carry it over untouched rather than losing its mass
            EN = np.where(orphan, m.EN, EN)
        m.ES *= (EN / np.maximum(m.EN, 1e-300))[:, :, None, None]
        m.EGEN *= EN / np.maximum(m.EN, 1e-300)
        m.NJ *= (EN[:, JUV] / np.maximum(m.EN[:, JUV], 1e-300))[:, None]
        m.EN = EN
        m.EB = np.where(orphan, m.EB, np.bincount(ecc, ns * (Wn + e_live), sz).reshape(sh))
        m.ER = np.where(orphan, m.ER, np.bincount(ecc, ns * e_live, sz).reshape(sh))
        ad = st == ADULT
        gi = np.ravel_multi_index((s[ad], u["j"][e][ad], u["i"][e][ad]), m.EG.shape)
        has = np.bincount(gi, None, m.EG.size).reshape(m.EG.shape) > 0
        m.EG = np.where(has, np.bincount(gi, ns[ad] * Gn[ad], m.EG.size).reshape(m.EG.shape), m.EG)
        w = np.bincount(ecc, r["ns"][e], sz)
        m.Tclass = (np.bincount(ecc, r["ns"][e] * r["T"][e], sz) / np.maximum(w, 1e-300)).reshape(sh)
        gpos = (np.bincount(ecc, ns * soma, sz) / np.maximum(EN.ravel(), 1e-300)).reshape(sh)
        m.promote_euler(gpos)

    # ------------------------------------------------ life history
    def promote_euler(m, gpos):
        S, Zb = m.sp, m.zbar()
        W = (m.EB - m.ER) / np.maximum(m.EN, 1e-300)
        Wj, Wm = S.weight(S.L_juv)[:, None, None], S.weight(Zb[:, JUV, L_MAT])
        We = (S.mu[:, EGG_M] * (1 - YOLK))[:, None, None]
        esc = lambda st, lo, hi: np.where(W[:, st] >= hi, 1, np.clip(gpos[:, st] / np.maximum(hi - lo, 1e-12), 0, 1))
        dur = S.egg_days[:, None, None] * Q10 ** (-(m.Tclass[:, EGG] - S.T_ref[:, None, None]) / 10)
        # juveniles: advance through the progress bins at K times the one-class rate; only the last bin matures
        K = m.NJ.shape[1]
        tot = m.NJ.sum(1)                                     # keep the bins summing to the juvenile count
        m.NJ = np.where((tot > 0)[:, None], m.NJ * (m.EN[:, JUV] / np.maximum(tot, 1e-300))[:, None],
                        (np.arange(K) == 0)[None, :, None, None] * m.EN[:, JUV, None])
        r = np.where(W[:, JUV] >= Wm, 1, np.clip(K * gpos[:, JUV] / np.maximum(Wm - Wj, 1e-12), 0, 1))[:, None]
        flow = m.NJ * r
        fj = np.clip(flow[:, -1] / np.maximum(m.EN[:, JUV], 1e-300), 0, 1)
        m.NJ -= flow
        m.NJ[:, 1:] += flow[:, :-1]
        fl = esc(LARVA, We, Wj)
        m.NJ[:, 0] += m.EN[:, LARVA] * fl                     # new juveniles start at the first bin
        for st, f in ((JUV, fj), (LARVA, fl), (EGG, np.minimum(1, m.dt / dur))):
            for X, ex in ((m.EN, 0), (m.EB, 0), (m.ER, 0), (m.ES, 2), (m.EGEN, 0)):
                d = X[:, st] * f.reshape(f.shape[:1] + (1,) * ex + f.shape[1:])
                X[:, st] -= d
                X[:, st + 1] += d

    def develop_agents(m):
        a, S, g = m.ag, m.sp, m.g
        if not len(a):
            return
        _, _, cell = m.cols()
        T = m.f["T"].ravel()[cell]
        egg = a.st == EGG
        a.dev = a.dev + egg * m.dt / (S.egg_days[a.sp] * Q10 ** (-(T - S.T_ref[a.sp]) / 10))
        L = S.length(a.W, a.sp)
        up = np.select([egg, a.st == LARVA, a.st == JUV], [a.dev >= 1, L >= S.L_juv[a.sp], L >= a.Z[:, L_MAT]], False)
        a.st = a.st + up
        a.age = a.age + m.dt
        old = a.age > a.amax
        m.toD += np.bincount(cell, old * a.n * (a.W + a.E + a.G), g.nc3)
        m.loss[:, CAUSES.index("senescence")] += np.bincount(a.sp, old * a.n * (a.W + a.E + a.G), S.n)
        m.remove(old, CAUSES.index("senescence"))

    def mate(m, fem, males, col):
        """A mate for every ripe female, as an index into `males` (-1 = none, she does not spawn).

        `genetics.mate_search` sets how far she looks: "column" (default) is her own grid cell, so gene flow
        is limited by dispersal and populations can diverge; "global" is the old behaviour, any conspecific
        male in the domain, which homogenises the species. With `mate_choice_sd` she samples
        `mate_choice_sample` males and prefers those like herself in `mate_choice_traits` (assortative
        mating, in standardised breeding-value units), the ingredient that lets divergence close on itself."""
        a, g, rng, gc = m.ag, m.g, m.rng, m.cfg.get("genetics", {})
        sd, ns = float(gc.get("mate_choice_sd", 0.0)), int(gc.get("mate_choice_sample", 1))
        keys, query = a.sp[males] * g.ncol + col[males], a.sp[fem] * g.ncol + col[fem]
        if sd > 0 and ns > 1 and len(males):
            cand = pick(rng, keys, query, ns)                                  # (nfem, ns) indices into males
            ti = [TRAITS.index(t) for t in gc.get("mate_choice_traits", ["T_opt"])]
            dz = a.A[males[np.maximum(cand, 0)]][:, :, ti] - a.A[fem][:, None, ti]
            w = np.exp(-0.5 * (dz ** 2).sum(-1) / sd ** 2) * (cand >= 0)
            cw = np.cumsum(w, 1)
            hit = np.minimum((cw < rng.random(len(fem))[:, None] * cw[:, -1:]).sum(1), ns - 1)
            dad = np.where(cw[:, -1] > 0, cand[np.arange(len(fem)), hit], -1)
        else:
            dad = pick(rng, keys, query)
        if gc.get("mate_search", "column") == "global":
            dad = np.where(dad < 0, pick(rng, a.sp[males], a.sp[fem]), dad)
        return dad

    def reproduce_agents(m):
        a, S, g, rng, hy = m.ag, m.sp, m.g, m.rng, m.hy
        if not len(a):
            return
        j, i, cell = m.cols()
        col = j * g.nx + i
        T = m.f["T"].ravel()[cell]
        phi = np.exp(-0.5 * ((T - a.Z[:, T_OPT]) / a.Z[:, T_WID]) ** 2)
        ripe = (a.st == ADULT) & (a.G >= S.gonad_frac[a.sp] * a.W) \
            & in_season((m.f["doy"] - a.Z[:, SPAWN]) % 365, S.spawn_doy[a.sp]) \
            & (phi > 0.3) & (a.n >= 1)
        male = ripe & (a.sex == 1)
        m.toD += np.bincount(cell, male * a.n * a.G, g.nc3)
        a.G = np.where(male, 0, a.G)
        fem = np.nonzero(ripe & (a.sex == 0))[0]
        males = np.nonzero((a.st == ADULT) & (a.sex == 1))[0]
        dad = m.mate(fem, males, col)
        k = hy.get("offspring_agents", 2)
        em = a.Z[fem, EGG_M]
        each = np.floor(a.n[fem] * a.G[fem] / em / k) * (dad >= 0)
        ok = each >= 1
        fem, dad, em, each = fem[ok], males[dad[ok]], em[ok], each[ok]
        if not len(fem):
            return
        a.G[fem] -= each * k * em / a.n[fem]
        r = np.repeat(np.arange(len(fem)), k)
        f, d = fem[r], dad[r]
        s = a.sp[f]
        if m.gen.L:
            gl = m.gen.meiosis(a.gl[f], a.gl[d], s, rng)
            A = m.gen.breeding(gl, s)
        else:
            gl, A = None, inherit(a.A[f], a.A[d], S.h2[s], rng)
        n, e = each[r], em[r]
        m.eggs += np.bincount(s, n, S.n)
        room = max(0, hy.get("max_agents", 10 ** 5) - len(a))
        if len(r) > room:
            x = slice(room, None)
            m.to_euler(s[x], np.full(len(r) - room, EGG), j[f][x], i[f][x], n[x], e[x], 0 * e[x], A[x], YOLK * e[x],
                       (np.maximum(a.gen[f], a.gen[d]) + 1)[x])
        x = slice(0, room)
        f, d, s, A, n, e = f[x], d[x], s[x], A[x], n[x], e[x]
        gl = gl[x] if gl is not None else None
        if len(f):
            lon = np.clip(a.lon[f] + (rng.random(len(f)) - 0.5) * 0.2 * g.dlon, g.lon_e[i[f]], g.lon_e[i[f] + 1] - 1e-9)
            lat = np.clip(a.lat[f] + (rng.random(len(f)) - 0.5) * 0.2 * g.dlat, g.lat_e[j[f]], g.lat_e[j[f] + 1] - 1e-9)
            a.add(sp=s, st=EGG, sex=rng.integers(0, 2, len(f)), k=0, n=n, W=e * (1 - YOLK), E=e * YOLK,
                  amax=S.lifespan_days[s] * np.exp(0.15 * rng.standard_normal(len(f))), tb=m.t, lon=lon, lat=lat,
                  lon0=lon, lat0=lat, gen=np.maximum(a.gen[f], a.gen[d]) + 1, mother=a.id[f], father=a.id[d],
                  origin=2, r0=np.log(n * e), **(dict(A=A, Z=m.phenotype(s, A), gl=gl) if gl is not None
                                                  else dict(A=A, Z=m.phenotype(s, A))))
            m.born += len(f)

    def reproduce_euler(m):
        S, g = m.sp, m.g
        Zb = m.zbar()
        Wa = (m.EB - m.ER)[:, ADULT]
        phi = np.exp(-0.5 * ((m.Tclass[:, ADULT] - Zb[:, ADULT, T_OPT]) / Zb[:, ADULT, T_WID]) ** 2)
        season = in_season((m.f["doy"] - Zb[:, ADULT, SPAWN]) % 365, S.spawn_doy[:, None, None])  # local timing
        go = season & (phi > 0.3) & (m.EG >= S.gonad_frac[:, None, None] * Wa) & (m.EN[:, ADULT] > 0)
        if not go.any():
            return
        eggs = np.where(go, 0.5 * m.EG / Zb[:, ADULT, EGG_M], 0)
        off = offspring_moments(np.moveaxis(m.ES[:, ADULT], 2, -1), m.EN[:, ADULT][:, None], S.h2[:, :, None, None])
        m.ES[:, EGG] += np.moveaxis(off, -1, 2) * eggs[:, None, None]
        m.EGEN[:, EGG] += eggs * (m.EGEN[:, ADULT] / np.maximum(m.EN[:, ADULT], 1e-300) + 1)
        m.EN[:, EGG] += eggs
        m.EB[:, EGG] += 0.5 * m.EG * go
        m.ER[:, EGG] += 0.5 * YOLK * m.EG * go
        m.toD[:g.ncol] += (0.5 * m.EG * go).sum(0).ravel()
        m.eggs += eggs.sum((1, 2))
        m.EG = np.where(go, 0, m.EG)

    def learn(m):
        """Within-lifetime REINFORCE on a learned offset to the innate policy (reward = log biomass change)."""
        a = m.ag
        if m.be["mode"] != "learn" or not len(a):
            return
        r = np.log(np.maximum(a.n * (a.W + a.E + a.G), 1e-300)) - a.r0
        upd = np.any(a.grad != 0, 1)
        a.dl = a.dl + upd[:, None] * m.be.get("learning_rate", 0.05) * (r - a.base)[:, None] * a.grad
        a.base = a.base + upd * 0.1 * (r - a.base)

    def cleanup(m):
        """Biomass classes with less than one fish become detritus (keeps the fields well conditioned)."""
        tiny = (m.EN < m.hy.get("min_class_fish", 1.0)) & (m.EN > 0) | (m.EN <= 0) & (m.EB != 0)
        mass = np.where(tiny, m.EB, 0) + np.where(tiny[:, ADULT], m.EG, 0)[:, None] * (np.arange(4) == ADULT)[None, :, None, None]
        m.toD[:m.g.ncol] += mass.sum((0, 1)).ravel()
        for X in (m.EN, m.EB, m.ER, m.EGEN):
            X[tiny] = 0
        m.NJ[np.broadcast_to(tiny[:, JUV, None], m.NJ.shape)] = 0
        m.ES[np.broadcast_to(tiny[:, :, None, None], m.ES.shape)] = 0
        m.EG[tiny[:, ADULT]] = 0

    # ------------------------------------------------ ocean coupling and main loop
    def ocean_step(m):
        p, f, g = m.npzd, m.f, m.g
        excr = (m.toN * RHO_N).reshape(g.nz, g.ny, g.nx) / g.volx
        p.C[0] += excr
        p.C[4] += (m.toD * RHO_N).reshape(g.nz, g.ny, g.nx) / g.volx
        p.C[O2_IDX] = np.maximum(p.C[O2_IDX] - O2N * excr, 0)                  # fish respiration
        m.o2_ext += p.air_sea(f["T"], m.cfg["ocean"].get("wind_speed", 7.0), m.dt)
        for _ in range(p.c["substeps"]):
            p.biology(f["T"], m.light, m.dt / p.c["substeps"])
        m.ext += p.restore(f["upwell"], m.dt)
        m.o2_ext += getattr(p, "o2_in", 0.0)
        p.C = m.transport3d(p.C * g.vol, f["Qx"], f["Qy"], p.c["kh"]) / g.volx
        p.vdiff(f["Kz"], m.dt_s)
        m.ext -= p.sink(m.dt)
        kref = int(np.searchsorted(g.z_e[1:], p.c.get("export_depth", 100.0)))
        deep = np.zeros((g.nz, 1, 1), bool)
        deep[kref + 1:] = True                                      # fish inputs already below the reference depth
        toD, toN = m.toD.reshape(g.nz, g.ny, g.nx), m.toN.reshape(g.nz, g.ny, g.nx)
        m.exp_sink = p.export_ref * C2N                             # mmol C per step through the reference depth
        m.exp_fish = (toD * deep).sum() * RHO_N * C2N               # faeces and carcasses released below it
        m.exp_resp = (toN * deep).sum() * RHO_N * C2N               # respired below it (migrators' active transport)

    def step(m):
        m.reset_step_tallies()
        m.force()
        if m.t >= m.next_reg - 1e-9:
            m.regulate()
        m.fields()
        m.vertical()
        m.move_euler()
        m.move_agents()
        m.sync()
        m.Zb = m.zbar()
        u = m.units()
        r = m.feed_and_die(u)
        m.grow_agents(u, r)
        m.grow_euler(u, r)
        m.learn()
        m.develop_agents()
        m.reproduce_agents()
        m.reproduce_euler()
        m.sync()
        m.cleanup()
        m.ocean_step()
        m.t += m.dt
        m.step_i += 1

    def diagnostics(m):
        S, a, g = m.sp, m.ag, m.g
        adult = (np.arange(4) == ADULT)[None, :]
        bio = (m.EB + m.EG[:, None] * adult[..., None, None]).sum((2, 3))
        num = m.EN.sum((2, 3))
        PS = m.ES.sum((1, 4, 5))
        abio, acnt = np.zeros(S.n), np.bincount(a.sp, None, S.n)
        if len(a):
            cls = a.sp * 4 + a.st
            abio = np.bincount(a.sp, a.n * (a.W + a.E + a.G), S.n)
            bio = bio + np.bincount(cls, a.n * (a.W + a.E + a.G), 4 * S.n).reshape(S.n, 4)
            num = num + np.bincount(cls, a.n, 4 * S.n).reshape(S.n, 4)
            PS = PS + np.stack([np.bincount(a.sp, a.n * a.A[:, t] ** p, S.n) for t in range(len(TRAITS)) for p in range(1, 5)],
                               1).reshape(S.n, len(TRAITS), 4)
        mean, var = central_moments(PS, num.sum(1)[:, None])[:2]
        acting = (a.beh >= 0) & (a.st >= LARVA)
        beh = np.bincount(a.sp[acting] * len(BEHAVIORS) + a.beh[acting], a.n[acting], S.n * len(BEHAVIORS)).reshape(S.n, -1)
        tot = m.total_N()
        return dict(time=m.t, biomass=bio, numbers=num, agent_biomass=abio, agents=acnt, losses=m.loss / m.dt,
                    catch=m.catch / m.dt, eggs=m.eggs / m.dt, trait_mean=S.mu + S.sd * mean, trait_sd=S.sd * np.sqrt(var),
                    plankton=(m.npzd.C[:N_NUTR] * g.vol).sum((1, 2, 3)), oxygen=(m.npzd.C[O2_IDX] * g.vol).sum(),
                    o2_min=m.npzd.C[O2_IDX][g.wet].min(), o2_ext=m.o2_ext,
                    metabolic_index=m.mi_sum / np.maximum(m.mi_n, 1e-300), **m.adaptation(),
                    export_sinking=m.exp_sink / m.dt, export_fish=m.exp_fish / m.dt, export_respired=m.exp_resp / m.dt, total_N=tot, ext_N=m.ext,
                    budget_error=(tot - m.total0 - m.ext) / m.total0, ibm_cells=m.ibm.sum(),
                    behavior=beh / np.maximum(beh.sum(1, keepdims=True), 1e-300),
                    max_gen=m.max_generation())

    def max_generation(m):
        """Furthest generation per species: agents' own counts, and the (rounded) mean generation of any
        biomass class holding at least one fish."""
        a, S = m.ag, m.sp
        ag = np.array([a.gen[a.sp == q].max() if (a.sp == q).any() else 0 for q in range(S.n)])
        live = m.EN >= m.hy.get("min_class_fish", 1.0)
        eu = np.where(live, m.EGEN / np.maximum(m.EN, 1e-300), 0).reshape(S.n, -1).max(1)
        return np.maximum(ag, np.round(eu)).astype(int)

    def log(m, d, wall):
        S, eta = m.sp, wall * (m.nsteps - m.step_i)
        now = (m.start + dtm.timedelta(days=m.t)).strftime("%Y-%m-%d %H:%M")
        P, Zp, Kp = m.npzd.C[1, 0][m.g.mask].mean(), m.npzd.C[2, 0][m.g.mask].mean(), m.npzd.C[3, 0][m.g.mask].mean()
        print(f"[{m.step_i:5d}/{m.nsteps}] {now}  day {m.t:7.2f} | agents {len(m.ag):7,d} (+{m.born} -{m.died}) "
              f"| IBM cells {d['ibm_cells']:4d} | sfc P {P:4.2f} Z {Zp:4.2f} K {Kp:4.2f} | N err {d['budget_error']:+.1e} "
              f"| {wall:5.2f} s/step  ETA {eta / 60:5.1f} min", flush=True)
        if m.verbose >= 2 or m.step_i % max(1, m.cfg["run"].get("species_log_every", 10)) == 0:
            for q, name in enumerate(S.names):
                b = d["biomass"][q]
                print(f"      {name:10s} B {b.sum() / 1e6:10.3e} t [" + " ".join(f"{x / max(b.sum(), 1e-30):4.0%}" for x in b)
                      + f"] agents {d['agents'][q]:6d} ({d['agent_biomass'][q] / max(b.sum(), 1e-30):4.0%} B) "
                      f"eggs/d {d['eggs'][q]:9.2e} catch/d {d['catch'][q] / 1e6:8.2e} t | T_opt {d['trait_mean'][q, T_OPT]:5.2f}"
                      f"±{d['trait_sd'][q, T_OPT]:4.2f} L_mat {d['trait_mean'][q, L_MAT]:6.2f} gen≤{d['max_gen'][q]}", flush=True)

    def banner(m):
        S, g, c = m.sp, m.g, m.cfg
        print("=" * 100)
        print(f"fishNET  run '{c['run']['name']}'  start {m.start:%Y-%m-%d %H:%M}  {m.nsteps} steps x {c['run']['dt_hours']} h"
              f"  ({c['run']['days']} days)")
        print(f"grid {g.nx}x{g.ny}x{g.nz}  ({g.mask.sum()} ocean columns, depth {g.z_e[-1]:.0f} m)  ocean: {c['ocean']['source']}"
              f"  hybrid: {m.hy['mode']}  behaviour: {m.be['mode']}  heritable: {c.get('genetics', {}).get('heritable', True)}")
        for q, n in enumerate(S.names):
            print(f"  {n:10s} B0 {S.biomass[q]:5.2f} g/m2  L_inf {S.L_inf[q]:6.1f} cm  W_inf {S.W_inf[q]:9.1f} g  "
                  f"T_opt {S.mu[q, T_OPT]:4.1f}  F {S.F[q]:.2f}/yr  diet: " +
                  ", ".join(t for t, v in zip(PLANKTON + S.names, np.r_[S.taxo[q, ADULT, :N_PLANK], S.taxo[q, ADULT, N_PLANK::4]]) if v > 0))
        for who, key, val in S.diet_unmatched:
            print(f"  note: {who}'s fish_diet lists '{key}' = {val}, which is not a species in this run; it is ignored")
        print(f"initial agents {len(m.ag):,d} in {m.ibm.sum()} IBM cells | total N {m.total0:.4e} mmol")
        print("=" * 100, flush=True)

    def spin_up(m):
        """Settle the fish and plankton before day 0: step for run.spinup_days while the forcing replays its first
        run.spinup_cycle_days, writing nothing. The clock runs from -spinup_days to 0 (so nothing is fished, and
        fish born now have negative birth days); budgets and life histories then start from the settled state."""
        r = m.cfg["run"]
        n = int(round(r.get("spinup_days", 30) / m.dt))
        if n <= 0:
            return
        cycle = max(float(r.get("spinup_cycle_days", 7.0)), m.dt)
        sw = getattr(m.ocean, "sw", None)
        keep = [x.copy() for x in (sw.h, sw.u, sw.v)] + [sw.t] if sw is not None else None
        if m.verbose:
            print(f"spin-up: {n * m.dt:g} days replaying the first {cycle:g} days of forcing (not written to output)",
                  flush=True)
        m.t, m.spin, m.next_reg = -n * m.dt, (n * m.dt, cycle), -n * m.dt
        tic, every = time.time(), max(1, int(round(5 / m.dt)))
        for k in range(n):
            m.step()
            if m.verbose and ((k + 1) % every == 0 or k + 1 == n):
                print(f"  spin-up day {m.t:+7.2f} | agents {len(m.ag):7,d} | fish {m.fish_mass() / 1e12:.4g} Mt "
                      f"| {time.time() - tic:.0f} s", flush=True)
        if sw is not None:                                  # the physics restarts where it was at day 0
            sw.h, sw.u, sw.v, sw.t = keep
        m.t, m.step_i, m.spin = 0.0, 0, None
        m.next_reg = min(m.next_reg, 0.0)
        m.obits, m.ext, m.o2_ext = [], 0.0, 0.0
        m.total0 = m.total_N()

    def run(m):
        resumed = bool(m.out) and m.out.resumed
        if m.verbose and not resumed:
            m.banner()
        if not resumed:
            m.spin_up()
        out_every = max(1, int(round(m.cfg["run"].get("output_every_hours", 24) / m.cfg["run"]["dt_hours"])))
        log_every = max(1, m.cfg["run"].get("log_every", 1))
        ag_h = m.cfg["run"].get("agent_output_every_hours", 120)
        ag_every = max(1, int(round(ag_h / m.cfg["run"]["dt_hours"]))) if ag_h > 0 else 0
        ck_days = m.cfg["run"].get("restart_every_days", 30)
        ck_every = max(1, int(round(ck_days / m.dt))) if ck_days > 0 and m.out else 0
        d = m.diagnostics()
        if m.out and not resumed:                           # a resumed run's current step is already on file
            m.out.write(m, d, True, True)
        if resumed and m.step_i >= m.nsteps and m.verbose:
            print(f"nothing to run: the restart is at day {m.t:g} and run.days = {m.cfg['run']['days']}; raise run.days")
        t0 = time.time()
        while m.step_i < m.nsteps:
            tic = time.time()
            m.step()
            d = m.diagnostics()
            if m.out:
                m.out.write(m, d, m.step_i % out_every == 0, bool(ag_every) and m.step_i % ag_every == 0)
            if m.verbose and m.step_i % log_every == 0:
                m.log(d, time.time() - tic)
            if ck_every and m.step_i % ck_every == 0 and m.step_i < m.nsteps:
                m.checkpoint()
        if ck_every:
            m.checkpoint()                                  # the end state, so a finished run can be extended
        if m.out:
            m.out.close(m)
        if m.verbose:
            print(f"done: {m.nsteps} steps in {(time.time() - t0) / 60:.2f} min; final N budget error {d['budget_error']:+.2e}")
        return m

    # ---------------------------------------------------------------- restarts
    def checkpoint(m):
        """Write the full model state to <out>/restart/ for `--resume`. The ocean is not stored (it is rebuilt from
        the namelist) apart from the shallow-water layers, which carry state. The outputs are flushed first and their
        lengths recorded, so a resume drops whatever was written after this point."""
        r, d = m.cfg["run"], m.out.dir / "restart"
        d.mkdir(exist_ok=True)
        state = dict(version=1, day=m.t, step=m.step_i, cfg=m.cfg, out=m.out.snapshot(),
                     model={k: v for k, v in m.__dict__.items() if k not in ("out", "ocean", "cfg")})
        if hasattr(m.ocean, "sw"):
            state["sw"] = dict(h=m.ocean.sw.h, u=m.ocean.sw.u, v=m.ocean.sw.v, t=m.ocean.sw.t)
        path = d / f"restart_day{m.t:07.1f}.pkl"
        tmp = path.with_suffix(".tmp")                      # a crash mid-write never leaves a truncated restart
        with open(tmp, "wb") as f:
            pickle.dump(state, f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, path)
        keep = r.get("restart_keep", 2)
        for old in (sorted(d.glob("restart_day*.pkl"))[:-keep] if keep > 0 else []):
            old.unlink()
        if m.verbose:
            print(f"restart file written: {path} ({path.stat().st_size / 1e6:.0f} MB)", flush=True)

    @classmethod
    def resume(cls, cfg, out_dir, path, verbose=None):
        """Continue from a restart file. `cfg` is the current namelist: [run] days is the new end day (raise it to
        extend a finished run) and output/restart settings and the ocean come from it, while fish, plankton,
        species parameters and the random stream come from the restart. dt, start and grid must match."""
        st = load_restart(path)
        old, r = st["cfg"], cfg["run"]
        if (old["run"]["dt_hours"], old["run"].get("start"), old["grid"]) != (r["dt_hours"], r.get("start"), cfg["grid"]):
            sys.exit(f"--resume: {path} was written with a different run.dt_hours, run.start or [grid] than this namelist")
        m = cls.__new__(cls)
        m.__dict__.update(st["model"])
        if "EGEN" not in m.__dict__:                        # restarts written before generations were tracked in biomass
            m.EGEN = np.zeros_like(m.EN)
        if "NJ" not in m.__dict__:                          # ...or before juvenile progress bins
            m.NJ = m.EN[:, JUV, None] * np.ones((1, max(1, int(cfg["hybrid"].get("juvenile_bins", 1))), 1, 1)) \
                / max(1, int(cfg["hybrid"].get("juvenile_bins", 1)))
        m.cfg, m.hy, m.be = cfg, cfg["hybrid"], cfg["behavior"]
        m.verbose = r.get("verbose", 1) if verbose is None else verbose
        m.nsteps = int(round(r["days"] / m.dt))
        if m.verbose:
            print(f"resuming '{r['name']}' from {path} at day {m.t:g} (step {m.step_i}), running to day {r['days']}",
                  flush=True)
            if old["species"] != cfg["species"]:
                print("  note: [[species]] differs from the restart; species parameters come from the restart", flush=True)
        m.ocean = Ocean(cfg["ocean"], m.g, m.start, m.verbose)
        if "sw" in st and hasattr(m.ocean, "sw"):
            m.ocean.sw.h, m.ocean.sw.u, m.ocean.sw.v, m.ocean.sw.t = (st["sw"][k] for k in "huvt")
        m.out = Output(m, out_dir, resume=st["out"]) if out_dir else None
        return m


class Output:
    """NetCDF output: gridded fields, agent snapshots, per-step series and completed life histories."""

    def __init__(o, m, out_dir, resume=None):
        o.dir = Path(out_dir)
        o.dir.mkdir(parents=True, exist_ok=True)
        o.resumed = resume is not None
        if o.resumed:                                        # reopen for appending, cut back to the restart point
            o.F, o.Aa, o.Sr = [nc.Dataset(trim_nc(o.dir / n, *resume[n]), "a") for n in ("fields.nc", "agents.nc", "series.nc")]
            return
        g, S, c = m.g, m.sp, m.cfg
        attrs = dict(species=",".join(S.names), colors=",".join(S.color), traits=",".join(TRAITS),
                     stages=",".join(STAGES), causes=",".join(CAUSES), behaviors=",".join(BEHAVIORS), start=str(m.start),
                     dt_hours=c["run"]["dt_hours"], hybrid=m.hy["mode"], engine=m.be.get("engine", "softmax"),
                     behavior_mode=m.be["mode"], ocean_area=float((g.area * g.mask).sum()), config=repr(c))
        dims = dict(species=S.n, stage=4, trait=len(TRAITS), depth=g.nz, lat=g.ny, lon=g.nx, cause=len(CAUSES), pool=N_NUTR,
                    behavior=len(BEHAVIORS))
        o.F, o.Aa, o.Sr = [o._open(n, dims, attrs, g) for n in ("fields.nc", "agents.nc", "series.nc")]
        o.Aa.createDimension("rec", None)
        var = lambda ds, name, dims, units="", dtype="f4": ds.createVariable(name, dtype, dims, zlib=True).setncattr("units", units)
        F4 = ("time", "depth", "lat", "lon")
        for v, u, dims in (("mask", "1", ("lat", "lon")), ("bottom_depth", "m", ("lat", "lon"))):
            static = o.F.createVariable(v, "f4", dims, zlib=True)
            static.units = u
            static[:] = g.mask.astype(float) if v == "mask" else g.depth
        for v, u in (("temp", "degC"), ("N", "mmol N m-3"), ("P", "mmol N m-3"), ("Z", "mmol N m-3"),
                     ("krill", "mmol N m-3"), ("D", "mmol N m-3"), ("O2", "mmol O2 m-3"), ("par", "W m-2")):
            var(o.F, v, F4, u)
        var(o.F, "u", F4, "m s-1"); var(o.F, "v", F4, "m s-1"); var(o.F, "w", F4, "m s-1")
        var(o.F, "fish_biomass", ("time", "species", "stage", "lat", "lon"), "g m-2")
        var(o.F, "fish_numbers", ("time", "species", "stage", "lat", "lon"), "m-2")
        var(o.F, "agent_biomass", ("time", "species", "lat", "lon"), "g m-2")
        var(o.F, "trait_mean", ("time", "species", "trait", "lat", "lon"), "breeding value")
        var(o.F, "ibm", ("time", "lat", "lon"), "", "i1")
        for k, t in (("time", "f8"), ("id", "i8"), ("species", "i2"), ("stage", "i1"), ("sex", "i1"), ("gen", "i4"),
                     ("mother", "i8"), ("father", "i8"), ("origin", "i1"), ("behavior", "i1")):
            var(o.Aa, k, ("rec",), "", t)
        for k in ("lon", "lat", "depth", "n", "W", "E", "G", "L", "age"):
            var(o.Aa, k, ("rec",), "", "f8" if k == "n" else "f4")
        var(o.Aa, "A", ("rec", "trait"), "breeding value"); var(o.Aa, "Z", ("rec", "trait"), "phenotype")
        S3 = ("time", "species")
        for k, dm, u in (("biomass", S3 + ("stage",), "g"), ("numbers", S3 + ("stage",), ""), ("agent_biomass", S3, "g"),
                         ("agents", S3, ""), ("losses", S3 + ("cause",), "g d-1"), ("catch", S3, "g d-1"),
                         ("eggs", S3, "d-1"), ("trait_mean", S3 + ("trait",), ""), ("trait_sd", S3 + ("trait",), ""),
                         ("plankton", ("time", "pool"), "mmol N"), ("total_N", ("time",), "mmol N"),
                         ("ext_N", ("time",), "mmol N"), ("budget_error", ("time",), ""), ("ibm_cells", ("time",), ""),
                         ("max_gen", S3, ""), ("behavior", S3 + ("behavior",), "fraction of agent fish"),
                         ("oxygen", ("time",), "mmol O2"), ("o2_min", ("time",), "mmol O2 m-3"),
                         ("o2_ext", ("time",), "mmol O2"), ("metabolic_index", S3, ""),
                         ("export_sinking", ("time",), "mmol C d-1"), ("export_fish", ("time",), "mmol C d-1"),
                         ("export_respired", ("time",), "mmol C d-1"), ("fst", S3, ""),
                         ("g_var", S3 + ("trait",), "genetic variance"),
                         ("sel", S3 + ("cause", "trait"), "selection differential")):
            var(o.Sr, k, dm, u, "f8")
        for ds in (o.F, o.Aa, o.Sr):          # put the headers on disk now: readers (the dashboard) may open the
            ds.sync()                         # files during a spin-up, long before the first record is written

    def _open(o, name, dims, attrs, g):
        ds = nc.Dataset(o.dir / name, "w")
        ds.setncatts(attrs)
        for k, v in dims.items():
            ds.createDimension(k, v)
        if name != "agents.nc":
            ds.createDimension("time", None)
            ds.createVariable("time", "f8", ("time",)).units = f"days since {attrs['start']}"
            for k, v in (("lon", g.lon), ("lat", g.lat), ("depth", g.z)):
                ds.createVariable(k, "f8", (k,))[:] = v
        return ds

    def write(o, m, d, fields, agents):
        n = len(o.Sr.dimensions["time"])
        for k in o.Sr.variables:
            if k in d:
                o.Sr[k][n] = d[k]
        g, S, a = m.g, m.sp, m.ag
        if agents and len(a):
            r0 = len(o.Aa.dimensions["rec"])
            vals = dict(time=np.full(len(a), m.t), id=a.id, species=a.sp, stage=a.st, sex=a.sex, gen=a.gen, mother=a.mother,
                        father=a.father, origin=a.origin, behavior=a.beh, lon=a.lon, lat=a.lat, depth=g.z[a.k], n=a.n, W=a.W, E=a.E,
                        G=a.G, L=S.length(a.W, a.sp), age=a.age, A=S.mu[a.sp] + S.sd[a.sp] * a.A, Z=a.Z)
            for k, v in vals.items():
                o.Aa[k][r0:r0 + len(a)] = v
        if not fields:
            return
        for ds in (o.Sr, o.Aa):
            ds.sync()
        n, f = len(o.F.dimensions["time"]), m.f
        o.F["time"][n] = m.t
        for k, v in (("temp", f["T"]), ("N", m.npzd.C[0]), ("P", m.npzd.C[1]), ("Z", m.npzd.C[2]),
                     ("krill", m.npzd.C[3]), ("D", m.npzd.C[4]), ("O2", m.npzd.C[O2_IDX]),
                     ("par", m.light), ("u", m.uc), ("v", m.vc), ("w", m.wc)):
            o.F[k][n] = np.where(g.mask, v, np.nan)
        B, N = m.EB + m.EG[:, None] * (np.arange(4) == ADULT)[None, :, None, None], m.EN.copy()
        PS1, AB = m.ES[:, :, :, 0].sum(1), np.zeros((S.n, g.ny, g.nx))
        if len(a):
            j, i, _ = m.cols()
            idx, ix = np.ravel_multi_index((a.sp, a.st, j, i), B.shape), np.ravel_multi_index((a.sp, j, i), AB.shape)
            mass = a.n * (a.W + a.E + a.G)
            B = B + np.bincount(idx, mass, B.size).reshape(B.shape)
            N = N + np.bincount(idx, a.n, B.size).reshape(B.shape)
            AB = np.bincount(ix, mass, AB.size).reshape(AB.shape)
            PS1 = PS1 + np.stack([np.bincount(ix, a.n * a.A[:, t], AB.size).reshape(AB.shape) for t in range(len(TRAITS))], 1)
        Nt = N.sum(1)
        o.F["fish_biomass"][n] = np.where(g.mask, B / g.area, np.nan)
        o.F["fish_numbers"][n] = np.where(g.mask, N / g.area, np.nan)
        o.F["agent_biomass"][n] = np.where(g.mask, AB / g.area, np.nan)
        tm = S.mu[:, :, None, None] + S.sd[:, :, None, None] * PS1 / np.maximum(Nt, 1e-300)[:, None]
        o.F["trait_mean"][n] = np.where((Nt[:, None] > 0) & g.mask, tm, np.nan)
        o.F["ibm"][n] = m.ibm
        o.F.sync()

    def snapshot(o):
        """For a restart: each file's length along its record dimension, plus every variable's rows since its last
        complete HDF chunk. A crash while that chunk is being rewritten can corrupt it (that is how the 1e36 junk
        in killed runs arises), whereas complete chunks are never touched again, so trim_nc takes those from the file
        and the partial one from here."""
        snap = {}
        for name, ds, dim in (("fields.nc", o.F, "time"), ("agents.nc", o.Aa, "rec"), ("series.nc", o.Sr, "time")):
            ds.sync()
            n, tail = len(ds.dimensions[dim]), {}
            for k, v in ds.variables.items():
                if v.dimensions[:1] == (dim,):
                    c0 = n // v.chunking()[0] * v.chunking()[0]
                    tail[k] = (c0, np.ma.getdata(v[c0:n]))
            snap[name] = (dim, n, tail)
        return snap

    def close(o, m):
        m.remove(np.ones(len(m.ag), bool), CAUSES.index("alive"))
        L = nc.Dataset(o.dir / "lifehist.nc", "w")
        L.setncatts(dict(species=",".join(m.sp.names), traits=",".join(TRAITS), causes=",".join(CAUSES),
                         stages=",".join(STAGES)))
        rec = {{"sp": "species", "st": "stage"}.get(k, k): np.concatenate([r[k] for r in m.obits]) for k in m.obits[0]} if m.obits else {}
        L.createDimension("agent", len(rec.get("id", [])))
        L.createDimension("trait", len(TRAITS))
        for k, v in rec.items():
            x = L.createVariable(k, v.dtype if v.dtype != np.float64 else "f4", ("agent",) + (("trait",) if v.ndim == 2 else ()),
                                 zlib=True)
            x[:] = v
        for ds in (o.F, o.Aa, o.Sr, L):
            ds.close()


def trim_nc(path, dim, n, tail):
    """Rewrite `path` with the first `n` records along its unlimited dimension `dim`: whole chunks from the file,
    the partial last chunk from the restart (`tail`, see Output.snapshot). Anything written after the restart goes."""
    tmp = path.with_suffix(".resume.nc")
    with nc.Dataset(path) as src, nc.Dataset(tmp, "w", format=src.data_model) as dst:
        src.set_auto_mask(False)
        if len(src.dimensions[dim]) < max((c0 for c0, _ in tail.values()), default=0):
            sys.exit(f"--resume: {path} is shorter than when the restart was written; it cannot be continued")
        dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})
        for k, d in src.dimensions.items():
            dst.createDimension(k, None if d.isunlimited() else len(d))
        for k, v in src.variables.items():
            f, ch = v.filters() or {}, v.chunking()
            x = dst.createVariable(k, v.dtype, v.dimensions, zlib=bool(f.get("zlib")), complevel=f.get("complevel") or 4,
                                   shuffle=bool(f.get("shuffle")), chunksizes=None if ch == "contiguous" else ch,
                                   fill_value=v.getncattr("_FillValue") if "_FillValue" in v.ncattrs() else None)
            x.setncatts({a: v.getncattr(a) for a in v.ncattrs() if a != "_FillValue"})
            if v.dimensions[:1] != (dim,):
                x[:] = v[:]
                continue
            c0, rows = tail[k]
            blk = max(1, 2 ** 22 // max(1, int(np.prod(v.shape[1:]))))    # agents.nc and fields.nc run to gigabytes
            for i in range(0, c0, blk):
                x[i:min(i + blk, c0)] = v[i:min(i + blk, c0)]
            if n > c0:
                x[c0:n] = rows
    os.replace(tmp, path)
    return path


class _Unpickler(pickle.Unpickler):
    """Restart files name fishNET's classes after whichever module wrote them (`__main__` when run as a script)."""
    def find_class(self, module, name):
        if module in ("__main__", "fishnet") and name in globals():
            return globals()[name]
        return super().find_class(module, name)


def literature_biomass(cfg, verbose=1):
    """Replace each species' starting `biomass` with the literature estimate in run.literature_biomass
    (default observations/literature_biomass.toml, written by `python obsval.py --tables`). The table holds
    global tonnes, spread as tonnes x 1e6 / 3.6e14 m2 like the namelists. Species it lacks keep theirs."""
    path = find_file(cfg["run"].get("literature_biomass", "observations/literature_biomass.toml"))
    table = tomllib.loads(Path(path).read_text())
    kept = []
    for s in cfg["species"]:
        e = table.get(s["name"])
        if e and "biomass_t" in e:
            s["biomass"] = e["biomass_t"] * 1e6 / 3.6e14
        else:
            kept.append(s["name"])
    if verbose:
        print(f"initial biomass from {path}" + (f" (namelist value kept for: {', '.join(kept)})" if kept else ""),
              flush=True)


def load_restart(path):
    """The dict a checkpoint wrote: model state, config, day/step and the output bookkeeping."""
    with open(path, "rb") as f:
        return _Unpickler(f).load()


def latest_restart(out):
    files = sorted((out / "restart").glob("restart_day*.pkl"))
    if not files:
        sys.exit(f"--resume: no restart files in {out / 'restart'} (set run.restart_every_days > 0 and rerun)")
    return files[-1]


# ---------------------------------------------------------------- ensembles
def get_param(cfg, path):
    """Value at a dotted path: section.key, species.<name>.key, species.<name>.traits.T_opt.0 (list index)."""
    keys = path.split(".")
    match = lambda s: keys[1] == "*" or key_name(keys[1]) == key_name(s["name"])
    t = next(s for s in cfg["species"] if match(s)) if keys[0] == "species" else cfg
    for k in keys[2:] if keys[0] == "species" else keys:
        t = t[int(k)] if isinstance(t, list) else t[k]
    return t


def set_param(cfg, path, value, scale=False):
    """Set a dotted path (species.* sets every species); scale=True multiplies each target's own value."""
    keys = path.split(".")
    match = lambda s: keys[1] == "*" or key_name(keys[1]) == key_name(s["name"])
    targets = [s for s in cfg["species"] if match(s)] if keys[0] == "species" else [cfg]
    if not targets:
        sys.exit(f"ensemble: no species matches '{path}'")
    for t in targets:
        for k in (keys[2:] if keys[0] == "species" else keys)[:-1]:
            t = t[int(k)] if isinstance(t, list) else t.setdefault(k, {})
        k = int(keys[-1]) if isinstance(t, list) else keys[-1]
        t[k] = t[k] * value if scale else value


def ensemble_members(cfg):
    """Parameter sets from [ensemble.sweep]: lists, {min, max, n, log} ranges or {scale = [...]} of the base value."""
    e, sweep = cfg["ensemble"], cfg["ensemble"]["sweep"]
    names, method = list(sweep), e.get("method", "grid")

    def discrete(n):
        sp = sweep[n]
        if isinstance(sp, list):
            return sp
        if "scale" in sp:
            return [float(x) for x in sp["scale"]]     # factors, applied to each target's own value
        return [float(x) for x in (np.geomspace if sp.get("log") else np.linspace)(sp["min"], sp["max"], sp.get("n", 3))]

    if method == "grid":
        members = [dict(zip(names, combo)) for combo in itertools.product(*map(discrete, names))]
    elif method == "oat":
        members = [{}] + [{n: v} for n in names for v in discrete(n)]
    else:                                      # "lhs" (stratified) or "random" in the unit cube
        N, rng = e.get("members", 10), np.random.default_rng(e.get("seed", 0))
        u = (np.argsort(rng.random((N, len(names))), 0) + rng.random((N, len(names)))) / N if method == "lhs" \
            else rng.random((N, len(names)))
        members = []
        for row in u:
            mem = {}
            for q, n in enumerate(names):
                sp = sweep[n]
                if isinstance(sp, dict) and "min" in sp:
                    lo, hi = (np.log(sp["min"]), np.log(sp["max"])) if sp.get("log") else (sp["min"], sp["max"])
                    mem[n] = float(np.exp(lo + row[q] * (hi - lo)) if sp.get("log") else lo + row[q] * (hi - lo))
                else:
                    vals = discrete(n)
                    mem[n] = vals[min(int(row[q] * len(vals)), len(vals) - 1)]
            members.append(mem)
    seed = cfg["run"].get("seed", 0)
    return [{**m, "run.seed": seed + r} for m in members for r in range(e.get("replicates", 1))]


def run_member(job):
    """Run one member; failures are reported, not fatal, so one bad parameter set cannot sink the ensemble."""
    i, cfg, out = job
    t0 = time.time()
    try:
        cfg["run"]["restart_every_days"] = 0                  # members are rerun, not resumed
        d = Model(cfg, out).run().diagnostics()
        return i, time.time() - t0, d["biomass"].sum(1), d["budget_error"], ""
    except (Exception, SystemExit) as err:
        return i, time.time() - t0, None, np.nan, f"{type(err).__name__}: {err}"


def ensemble(cfg, out):
    import copy, multiprocessing
    e = cfg["ensemble"]
    members = ensemble_members(cfg)
    names = sorted({k for m in members for k in m})
    print(f"fishNET ensemble '{cfg['run']['name']}': {len(members)} members, method '{e.get('method', 'grid')}', "
          f"{e.get('processes', 1)} processes -> {out}")
    print("  swept: " + ", ".join(n for n in names if n != "run.seed"))
    jobs, oceans = [], {}
    for i, mem in enumerate(members):
        c = copy.deepcopy({k: v for k, v in cfg.items() if k != "ensemble"})
        for p, v in mem.items():
            set_param(c, p, v, isinstance(e["sweep"].get(p), dict) and "scale" in e["sweep"][p])
        c["run"].update(verbose=e.get("member_verbose", 0), days=e.get("days", c["run"]["days"]), name=f"member_{i:03d}")
        if not e.get("keep_agents", False):
            c["run"]["agent_output_every_hours"] = 0
        if c["ocean"]["source"] == "shallow_water":   # spin up each distinct ocean once, before the workers start
            oceans.setdefault(repr((c["ocean"], c["grid"], c["run"].get("start"))), c)
        jobs.append((i, c, out / f"member_{i:03d}"))
        print(f"  member {i:3d}: " + ", ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in mem.items()))
    for c in oceans.values():
        Ocean(c["ocean"], Grid(c["grid"]), dtm.datetime.fromisoformat(c["run"].get("start", "2000-01-01T00:00")))
    out.mkdir(parents=True, exist_ok=True)
    t0, done, results = time.time(), 0, {}
    with multiprocessing.Pool(e.get("processes", 1)) as pool:
        for i, wall, B, err, fail in pool.imap_unordered(run_member, jobs):
            done += 1
            results[i] = fail
            eta = (time.time() - t0) / done * (len(jobs) - done) / 60
            msg = f"FAILED ({fail})" if fail else "final biomass (t) " + " ".join(f"{b / 1e6:.3g}" for b in B) + f" | N err {err:+.1e}"
            print(f"  [{done:3d}/{len(jobs)}] member {i:3d} {wall / 60:5.1f} min | {msg} | ETA {eta:.1f} min", flush=True)
    collect_ensemble(out, members, names, [i for i, f in results.items() if not f])


def collect_ensemble(out, members, names, ok):
    """Stack every member's time series into ensemble.nc with the parameter values per member."""
    import xarray as xr
    series = {i: xr.open_dataset(out / f"member_{i:03d}" / "series.nc", decode_times=False) for i in sorted(ok)}
    first = next(iter(series.values()))
    nt = min(len(s.time) for s in series.values())
    keep = ["biomass", "numbers", "agents", "trait_mean", "trait_sd", "eggs", "catch", "loss", "plankton", "budget_error", "behavior"]
    ds = xr.Dataset({k: xr.concat([series[i][k].isel(time=slice(0, nt)) if i in series else
                                   xr.full_like(first[k].isel(time=slice(0, nt)), np.nan, dtype=float) for i in range(len(members))],
                                  "member") for k in keep if k in first})
    ds = ds.assign_coords(time=first.time.values[:nt], member=np.arange(len(members)))
    for q, n in enumerate(names):
        vals = [m.get(n) for m in members]
        num = all(isinstance(v, (int, float)) or v is None for v in vals)
        ds[f"param_{q}"] = ("member", np.array([np.nan if v is None else v for v in vals], float) if num
                            else np.array([str(v) for v in vals]))
        ds[f"param_{q}"].attrs["path"] = n
    ds.attrs = {**first.attrs, "params": ",".join(names), "failed": ",".join(str(i) for i in range(len(members)) if i not in ok)}
    ds.to_netcdf(out / "ensemble.nc")
    print(f"ensemble summary: {out / 'ensemble.nc'} ({len(ok)}/{len(members)} members succeeded)")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    cfg = load_config(sys.argv[1])
    cfg["ocean"].setdefault("sw_cache_dir", str(Path(sys.argv[1]).parent / "sw_cache"))
    out = Path(cfg["run"].get("out_dir", "runs")) / cfg["run"]["name"]
    if not out.is_absolute():
        out = Path(sys.argv[1]).parent / out
    if "--ensemble" in sys.argv:
        k = sys.argv.index("--ensemble") + 1                       # optional member count: --ensemble 40
        if k < len(sys.argv) and sys.argv[k].isdigit():
            cfg.setdefault("ensemble", {})["members"] = int(sys.argv[k])
            if cfg["ensemble"].get("method", "grid") in ("grid", "oat"):
                print(f"note: method '{cfg['ensemble']['method']}' sets its own member count; {sys.argv[k]} ignored")
        ensemble(cfg, out.with_name(out.name + "_ensemble"))
    elif "--resume" in sys.argv:
        k = sys.argv.index("--resume") + 1                         # optional file: --resume runs/x/restart/restart_day00365.0.pkl
        path = Path(sys.argv[k]) if k < len(sys.argv) and not sys.argv[k].startswith("--") else latest_restart(out)
        Model.resume(cfg, out, path).run()
        print(f"output in {out}")
    else:
        Model(cfg, out).run()
        print(f"output in {out}")


if __name__ == "__main__":
    main()
