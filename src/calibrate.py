#!/usr/bin/env python3
"""Fit fishNET species parameters to observations.

    python calibrate.py namelist.toml targets.toml --method abc --members 60 --processes 4
    python calibrate.py namelist.toml targets.toml --method smc --members 40 --rounds 4
    python calibrate.py namelist.toml targets.toml --method emulator --members 60 --chain 20000

Priors are declared in the namelist under [calibration.prior], using the same syntax as an ensemble
sweep, and observations in a separate TOML (see targets.toml). Every method writes posterior.nc with
the accepted parameter sets, their losses and predictions, plus a short report.

  abc       run the prior, keep the best `--keep` fraction (rejection ABC)
  smc       several rounds, each resampling and perturbing the survivors of the last (ABC-SMC)
  emulator  run the prior, fit a Gaussian process to the loss, then MCMC on the emulator (cheap,
            but only as good as the surface the design covers)
"""
import argparse, copy, sys, tomllib, multiprocessing
from pathlib import Path
import numpy as np
import fishnet as fn


# ---------------------------------------------------------------- observations
def load_targets(path):
    """Read observation targets: name, what to compare, the value, and its uncertainty."""
    t = tomllib.loads(Path(path).read_text())
    out = []
    for i, d in enumerate(t.get("target", [])):
        d = dict(d)
        d.setdefault("name", f"target_{i}")
        d.setdefault("weight", 1.0)
        d.setdefault("log", d["kind"] in ("biomass", "catch", "chlorophyll", "export"))
        if "sigma" not in d:
            sys.exit(f"target '{d['name']}' needs sigma (observation + representativeness uncertainty)")
        out.append(d)
    if not out:
        sys.exit(f"no [[target]] entries in {path}")
    return out


def predict(ds, targets, names):
    """Model equivalent of each target, from a run's series.nc."""
    t = ds["time"].values
    out = np.full(len(targets), np.nan)
    for i, d in enumerate(targets):
        lo, hi = d.get("day_range", [t.max() - 365, t.max()])
        w = (t >= lo) & (t <= hi)
        if not w.any():
            continue
        sp = names.index(d["species"]) if "species" in d else slice(None)
        k = d["kind"]
        if k == "biomass":
            x = ds["biomass"][w, sp].sum("stage")
        elif k == "catch":
            x = ds["catch"][w, sp]
        elif k == "chlorophyll":
            x = ds["plankton"][w, 1]
        elif k == "export":
            x = ds["export_sinking"][w]
        elif k == "metabolic_index":
            x = ds["metabolic_index"][w, sp]
        elif k in ds:
            x = ds[k][w, sp] if ds[k].ndim > 1 else ds[k][w]
        else:
            sys.exit(f"target '{d['name']}': unknown kind '{k}'")
        x = np.asarray(x, float)
        out[i] = np.nanmean(x) * d.get("scale", 1.0)
    return out


def loss(pred, targets):
    """-2 log likelihood, Gaussian in the target's own units (log units where the target says so)."""
    tot = 0.0
    for p, d in zip(pred, targets):
        obs, sig = float(d["value"]), float(d["sigma"])
        if not np.isfinite(p) or (d["log"] and p <= 0):
            return np.inf
        if d["log"]:                       # in log space a relative sigma is already a log-sd
            p, obs = np.log(p), np.log(obs)
            sig = sig if d.get("sigma_relative", True) else sig / np.exp(obs)
        tot += d["weight"] * ((p - obs) / sig) ** 2
    return float(tot)


# ---------------------------------------------------------------- priors and sampling
def prior_spec(cfg):
    p = cfg.get("calibration", {}).get("prior")
    if not p:
        sys.exit("no [calibration.prior] section in the namelist")
    names = list(p)
    bounds, logs = [], []
    for n in names:
        d = p[n]
        if not isinstance(d, dict) or "min" not in d:
            sys.exit(f"prior '{n}' must be a table with min and max")
        bounds.append((float(d["min"]), float(d["max"])))
        logs.append(bool(d.get("log", False)))
    return names, np.array(bounds), np.array(logs)


def to_unit(x, bounds, logs):
    lo, hi = np.where(logs, np.log(bounds[:, 0]), bounds[:, 0]), np.where(logs, np.log(bounds[:, 1]), bounds[:, 1])
    return (np.where(logs, np.log(np.maximum(x, 1e-300)), x) - lo) / (hi - lo)


def from_unit(u, bounds, logs):
    lo, hi = np.where(logs, np.log(bounds[:, 0]), bounds[:, 0]), np.where(logs, np.log(bounds[:, 1]), bounds[:, 1])
    x = lo + np.clip(u, 0, 1) * (hi - lo)
    return np.where(logs, np.exp(x), x)


def lhs(n, d, rng):
    return (np.argsort(rng.random((n, d)), 0) + rng.random((n, d))) / n


# ---------------------------------------------------------------- running members
def _run(job):
    i, cfg, out, targets, names = job
    try:
        import xarray as xr
        fn.Model(cfg, out).run()
        ds = xr.open_dataset(out / "series.nc", decode_times=False)
        pred = predict(ds, targets, names)
        return i, pred, loss(pred, targets), ""
    except (Exception, SystemExit) as err:
        return i, np.full(len(targets), np.nan), np.inf, f"{type(err).__name__}: {err}"


def evaluate(theta, cfg, names, targets, out, a, tag, start=0):
    """Run one batch of parameter sets and return their predictions and losses."""
    jobs = []
    for i, row in enumerate(theta):
        c = copy.deepcopy({k: v for k, v in cfg.items() if k not in ("ensemble", "calibration")})
        for n, v in zip(names, row):
            fn.set_param(c, n, float(v))
        c["run"].update(verbose=0, days=a.days or c["run"]["days"], name=f"{tag}_{start + i:03d}",
                        agent_output_every_hours=0, seed=c["run"].get("seed", 0) + (start + i))
        jobs.append((i, c, out / f"{tag}_{start + i:03d}", targets, [s["name"] for s in c["species"]]))
    pred = np.full((len(theta), len(targets)), np.nan)
    L = np.full(len(theta), np.inf)
    with multiprocessing.Pool(a.processes) as pool:
        for k, (i, p, l, err) in enumerate(pool.imap_unordered(_run, jobs)):
            pred[i], L[i] = p, l
            print(f"  [{k + 1:3d}/{len(jobs)}] {tag} {start + i:3d}  loss {l:11.3f}" + (f"  FAILED {err[:60]}" if err else ""),
                  flush=True)
    return pred, L


# ---------------------------------------------------------------- Gaussian process emulator
class GP:
    """Zero-mean GP on the unit cube with an RBF kernel; length scale and noise picked by marginal likelihood."""

    def __init__(gp, X, y):
        gp.X, gp.mu, gp.sd = X, y.mean(), y.std() + 1e-12
        gp.y = (y - gp.mu) / gp.sd
        d2 = ((X[:, None] - X[None]) ** 2).sum(-1)
        best = None
        for ell in np.geomspace(0.1, 3.0, 12):
            for noise in (1e-3, 1e-2, 1e-1):
                K = np.exp(-0.5 * d2 / ell ** 2) + noise * np.eye(len(X))
                try:
                    Lc = np.linalg.cholesky(K)
                except np.linalg.LinAlgError:
                    continue
                al = np.linalg.solve(Lc.T, np.linalg.solve(Lc, gp.y))
                ml = -0.5 * gp.y @ al - np.log(np.diag(Lc)).sum()
                if best is None or ml > best[0]:
                    best = (ml, ell, noise, Lc, al)
        gp.ml, gp.ell, gp.noise, gp.L, gp.alpha = best

    def __call__(gp, Xs, var=False):
        k = np.exp(-0.5 * ((Xs[:, None] - gp.X[None]) ** 2).sum(-1) / gp.ell ** 2)
        mean = k @ gp.alpha * gp.sd + gp.mu
        if not var:
            return mean
        v = np.linalg.solve(gp.L, k.T)
        return mean, np.maximum(1 + gp.noise - (v ** 2).sum(0), 0) * gp.sd ** 2

    def loo_r2(gp):
        """Leave-one-out R^2 from the standard GP shortcut, to show whether the surface is trustworthy."""
        Ki = np.linalg.inv(gp.L.T) @ np.linalg.inv(gp.L)
        res = gp.alpha / np.diag(Ki)
        return float(1 - (res ** 2).sum() / ((gp.y - gp.y.mean()) ** 2).sum())


def mcmc(logp, d, n, rng, step=0.1):
    """Random-walk Metropolis on the unit cube."""
    x = np.full(d, 0.5)
    lp, keep, acc = logp(x[None])[0], [], 0
    for i in range(n):
        y = x + step * rng.standard_normal(d)
        if (y < 0).any() or (y > 1).any():
            keep.append(x.copy())
            continue
        ly = logp(y[None])[0]
        if np.log(rng.random()) < ly - lp:
            x, lp, acc = y, ly, acc + 1
        keep.append(x.copy())
    return np.array(keep[n // 5:]), acc / n


# ---------------------------------------------------------------- methods
def run_abc(cfg, names, bounds, logs, targets, out, a, rng):
    U = lhs(a.members, len(names), rng)
    theta = from_unit(U, bounds, logs)
    pred, L = evaluate(theta, cfg, names, targets, out, a, "prior")
    ok = np.isfinite(L)
    if ok.sum() < 2:
        sys.exit("every prior member failed; check the prior ranges")
    eps = np.quantile(L[ok], a.keep)
    return theta, pred, L, L <= eps, dict(epsilon=float(eps))


def run_smc(cfg, names, bounds, logs, targets, out, a, rng):
    U = lhs(a.members, len(names), rng)
    theta = from_unit(U, bounds, logs)
    pred, L = evaluate(theta, cfg, names, targets, out, a, "round0")
    keep = L <= np.quantile(L[np.isfinite(L)], a.keep)
    hist = [dict(round=0, alive=int(keep.sum()), best=float(np.nanmin(L)), epsilon=float(np.quantile(L[np.isfinite(L)], a.keep)))]
    for r in range(1, a.rounds):
        U0 = to_unit(theta[keep], bounds, logs)
        cov = np.cov(U0.T) if U0.shape[0] > len(names) else np.eye(len(names)) * 0.01
        cov = np.atleast_2d(cov) * 2.0 / max(len(names), 1) + 1e-4 * np.eye(len(names))
        idx = rng.integers(0, len(U0), a.members)
        Un = np.clip(U0[idx] + rng.multivariate_normal(np.zeros(len(names)), cov, a.members), 0, 1)
        tn = from_unit(Un, bounds, logs)
        pn, Ln = evaluate(tn, cfg, names, targets, out, a, f"round{r}", start=r * a.members)
        theta, pred, L = np.vstack([theta, tn]), np.vstack([pred, pn]), np.r_[L, Ln]
        eps = np.quantile(Ln[np.isfinite(Ln)], a.keep) if np.isfinite(Ln).any() else np.inf
        keep = L <= max(eps, np.nanmin(L) * (1 + 1e-9))
        hist.append(dict(round=r, alive=int(keep.sum()), best=float(np.nanmin(L)), epsilon=float(eps)))
        print(f"  round {r}: epsilon {eps:.3f}, best {np.nanmin(L):.3f}, {keep.sum()} kept", flush=True)
    return theta, pred, L, keep, dict(rounds=hist)


def run_emulator(cfg, names, bounds, logs, targets, out, a, rng):
    U = lhs(a.members, len(names), rng)
    theta = from_unit(U, bounds, logs)
    pred, L = evaluate(theta, cfg, names, targets, out, a, "design")
    ok = np.isfinite(L)
    if ok.sum() < len(names) + 3:
        sys.exit("too few successful members to fit an emulator")
    gp = GP(U[ok], np.log(L[ok] + 1e-9))
    r2 = gp.loo_r2()
    print(f"  emulator: length scale {gp.ell:.2f}, noise {gp.noise:g}, leave-one-out R2 {r2:.2f}")
    if r2 < 0.3:
        print("  warning: the emulator fits poorly, so treat the posterior as indicative only")
    chain, acc = mcmc(lambda x: -0.5 * np.exp(gp(x)), len(names), a.chain, rng)
    print(f"  MCMC: {len(chain)} samples kept, acceptance {acc:.2f}")
    post = from_unit(chain[:: max(1, len(chain) // 2000)], bounds, logs)
    return theta, pred, L, ok & (L <= np.quantile(L[ok], a.keep)), dict(posterior=post, loo_r2=r2)


# ---------------------------------------------------------------- reporting
def report(names, bounds, theta, pred, L, keep, targets, extra, out):
    import netCDF4 as nc
    post = extra.get("posterior", theta[keep])
    print(f"\nposterior from {keep.sum()} of {len(theta)} members (best loss {np.nanmin(L):.3f})")
    print(f"{'parameter':<38}{'prior range':>22}{'posterior median':>18}{'   90% interval':>22}")
    for j, n in enumerate(names):
        q = np.nanquantile(post[:, j], [0.05, 0.5, 0.95])
        print(f"{n:<38}{bounds[j, 0]:>10.4g}..{bounds[j, 1]:<11.4g}{q[1]:>18.4g}   {q[0]:>9.4g}..{q[2]:<9.4g}")
    best = np.nanargmin(L)
    print(f"\n{'target':<22}{'observed':>12}{'sigma':>10}{'best member':>14}{'posterior mean':>16}")
    for i, d in enumerate(targets):
        pm = np.nanmean(pred[keep, i]) if keep.any() else np.nan
        print(f"{d['name']:<22}{float(d['value']):>12.4g}{float(d['sigma']):>10.3g}{pred[best, i]:>14.4g}{pm:>16.4g}")
    with nc.Dataset(out / "posterior.nc", "w") as ds:
        ds.createDimension("member", len(theta))
        ds.createDimension("param", len(names))
        ds.createDimension("target", len(targets))
        ds.createDimension("sample", len(post))
        ds.createVariable("theta", "f8", ("member", "param"))[:] = theta
        ds.createVariable("loss", "f8", ("member",))[:] = np.where(np.isfinite(L), L, 1e30)
        ds.createVariable("accepted", "i1", ("member",))[:] = keep.astype(np.int8)
        ds.createVariable("prediction", "f8", ("member", "target"))[:] = np.nan_to_num(pred, nan=-9999)
        ds.createVariable("posterior", "f8", ("sample", "param"))[:] = post
        ds.createVariable("observed", "f8", ("target",))[:] = [float(d["value"]) for d in targets]
        ds.createVariable("sigma", "f8", ("target",))[:] = [float(d["sigma"]) for d in targets]
        ds.setncatts(dict(params=",".join(names), targets=",".join(d["name"] for d in targets),
                          best_loss=float(np.nanmin(L)), **{k: str(v) for k, v in extra.items() if k != "posterior"}))
    print(f"\nwritten to {out / 'posterior.nc'}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("namelist")
    p.add_argument("targets")
    p.add_argument("--method", choices=("abc", "smc", "emulator"), default="abc")
    p.add_argument("--members", type=int, default=40, help="members per round")
    p.add_argument("--rounds", type=int, default=3, help="smc only")
    p.add_argument("--keep", type=float, default=0.2, help="accepted fraction")
    p.add_argument("--chain", type=int, default=20000, help="emulator MCMC steps")
    p.add_argument("--days", type=float, default=0, help="override run length")
    p.add_argument("--processes", type=int, default=2)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default=None)
    a = p.parse_args()
    cfg = fn.load_config(a.namelist)
    cfg["ocean"].setdefault("sw_cache_dir", str(Path(a.namelist).parent / "sw_cache"))
    targets = load_targets(a.targets)
    names, bounds, logs = prior_spec(cfg)
    out = Path(a.out or Path(a.namelist).parent / "runs" / f"{cfg['run']['name']}_calibration")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    print(f"fishNET calibration '{a.method}': {len(names)} parameters, {len(targets)} targets, "
          f"{a.members} members{f' x {a.rounds} rounds' if a.method == 'smc' else ''}, {a.processes} processes")
    print("  parameters: " + ", ".join(names))
    print("  targets:    " + ", ".join(d["name"] for d in targets))
    if cfg["ocean"]["source"] == "shallow_water":          # spin the ocean up once, before the workers start
        fn.Ocean(cfg["ocean"], fn.Grid(cfg["grid"]), fn.dtm.datetime.fromisoformat(cfg["run"].get("start", "2000-01-01T00:00")))
    method = dict(abc=run_abc, smc=run_smc, emulator=run_emulator)[a.method]
    theta, pred, L, keep, extra = method(cfg, names, bounds, logs, targets, out, a, rng)
    report(names, bounds, theta, pred, L, keep, targets, extra, out)


if __name__ == "__main__":
    main()
