#!/usr/bin/env python3
"""Score a fishNET run against observations.

    python skill.py runs/atlantic_hybrid observations.toml
    python skill.py runs/atlantic_hybrid observations.toml --ensemble runs/x_ensemble/ensemble.nc

Reads `[[target]]` entries (a single number, as calibration uses) and `[[series]]` entries (a time
series of observations) from the observations file, and reports:

  point targets    bias, relative error, and whether the model sits inside the stated uncertainty
  time series      bias, RMSE, centred RMS, correlation, standard-deviation ratio, Nash-Sutcliffe
  size spectrum    slope and intercept of the community abundance spectrum against an observed slope
  events           marine heatwaves and population collapses: did the model produce them, and when
  reliability      rank histogram and CRPS, when an ensemble is supplied

It writes skill.md (the scorecard), skill.nc (the numbers) and a Taylor diagram.
"""
import argparse, sys, tomllib
from pathlib import Path
import numpy as np
import xarray as xr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------- metrics
def metrics(model, obs):
    """Standard deterministic skill metrics for a paired series."""
    m, o = np.asarray(model, float), np.asarray(obs, float)
    w = np.isfinite(m) & np.isfinite(o)
    m, o = m[w], o[w]
    if len(m) < 2:
        return {k: np.nan for k in ("n", "bias", "rmse", "crms", "corr", "sd_ratio", "nse")} | {"n": len(m)}
    bias = float(m.mean() - o.mean())
    rmse = float(np.sqrt(((m - o) ** 2).mean()))
    crms = float(np.sqrt((((m - m.mean()) - (o - o.mean())) ** 2).mean()))
    sd_m, sd_o = m.std(), o.std()
    corr = float(np.corrcoef(m, o)[0, 1]) if sd_m > 0 and sd_o > 0 else np.nan
    nse = float(1 - ((m - o) ** 2).sum() / max(((o - o.mean()) ** 2).sum(), 1e-300))
    return dict(n=len(m), bias=bias, rmse=rmse, crms=crms, corr=corr,
                sd_ratio=float(sd_m / sd_o) if sd_o > 0 else np.nan, nse=nse)


def crps(ens, obs):
    """Continuous ranked probability score from an ensemble sample (lower is better)."""
    e = np.sort(np.asarray(ens, float))
    n = len(e)
    if n == 0 or not np.isfinite(obs):
        return np.nan
    a = np.abs(e - obs).mean()
    b = np.abs(e[:, None] - e[None]).sum() / (2 * n * n)
    return float(a - b)


def rank_histogram(ens, obs):
    """Where the observation falls among the ensemble members (flat = reliable)."""
    return int((np.asarray(ens) < obs).sum())


# ---------------------------------------------------------------- model equivalents
def series_of(ds, d, names):
    """Model time series matching an observation entry, as (time, values)."""
    t = ds["time"].values
    k = d["kind"]
    sp = names.index(d["species"]) if "species" in d else None
    if k == "biomass":
        x = ds["biomass"][:, sp].sum("stage") if sp is not None else ds["biomass"].sum(("species", "stage"))
    elif k == "catch":
        x = ds["catch"][:, sp] if sp is not None else ds["catch"].sum("species")
    elif k == "chlorophyll":
        x = ds["plankton"][:, 1]
    elif k == "export":
        x = ds["export_sinking"]
    elif k in ds:
        x = ds[k][:, sp] if sp is not None and ds[k].ndim > 1 else ds[k]
    else:
        sys.exit(f"observation '{d['name']}': unknown kind '{k}'")
    return t, np.asarray(x, float) * d.get("scale", 1.0)


def size_spectrum(run, ds, names, nbins=12):
    """Community abundance size spectrum: log10 numbers per log10 mass bin, and its slope.
    Uses agent snapshots when present (individual masses), otherwise the stage structure."""
    W, N = [], []
    ag = run / "agents.nc"
    if ag.exists() and xr.open_dataset(ag, decode_times=False).sizes.get("rec", 0):
        a = xr.open_dataset(ag, decode_times=False)
        last = float(a["time"].max())
        w = a["time"].values == last
        W = (a["W"].values + a["E"].values)[w]
        N = a["n"].values[w]
    source = "agents"
    if len(W) < 20:                              # no agent snapshots: fall back to the stage structure
        source = "stage structure (coarse)"
        B, Nc = ds["biomass"][-30:].mean("time").values, ds["numbers"][-30:].mean("time").values
        W, N = (B / np.maximum(Nc, 1e-300)).ravel(), Nc.ravel()
    W, N = np.asarray(W, float), np.asarray(N, float)
    ok = (W > 0) & (N > 0)
    W, N = W[ok], N[ok]
    if len(W) < 5:
        return None
    edges = np.linspace(np.log10(W.min()), np.log10(W.max()) + 1e-9, nbins + 1)
    idx = np.clip(np.digitize(np.log10(W), edges) - 1, 0, nbins - 1)
    tot = np.bincount(idx, N, nbins)
    mid = 0.5 * (edges[1:] + edges[:-1])
    width = np.diff(edges)
    use = tot > 0
    if use.sum() < 3:
        return None
    y = np.log10(tot[use] / width[use])                     # normalised: numbers per unit log mass
    slope, intercept = np.polyfit(mid[use], y, 1)
    return dict(slope=float(slope), intercept=float(intercept), mid=mid[use], y=y, bins=int(use.sum()),
                source=source)


def find_heatwaves(F, threshold=1.5, min_days=5, area_fraction=0.02, control=None):
    """Marine heatwaves. The seasonal cycle is removed per cell with an annual + semiannual harmonic fit
    (so a warm event does not hide inside its own climatology), and an event is declared when more than
    `area_fraction` of the ocean exceeds `threshold` degC of anomaly for at least `min_days`."""
    t = F["time"].values
    if len(t) < 8:
        return []
    sst = F["temp"].isel(depth=0).values                                  # (time, lat, lon)
    wet = np.isfinite(sst[0]) & (F["mask"].values > 0 if "mask" in F else True)
    x = sst[:, wet]
    w = 2 * np.pi * t / 365
    A = np.column_stack([np.ones_like(t), np.cos(w), np.sin(w), np.cos(2 * w), np.sin(2 * w)])
    if control is not None:                                               # climatology from a control run
        ct = control["time"].values
        cs = control["temp"].isel(depth=0).values[:, wet]
        Ac = np.column_stack([np.ones_like(ct), np.cos(2 * np.pi * ct / 365), np.sin(2 * np.pi * ct / 365),
                              np.cos(4 * np.pi * ct / 365), np.sin(4 * np.pi * ct / 365)])
        anom = x - A @ np.linalg.lstsq(Ac, cs, rcond=None)[0]
    else:
        anom = x - A @ np.linalg.lstsq(A, x, rcond=None)[0]
    frac = (anom > threshold).mean(1)
    out, start = [], None
    for i, fr in enumerate(frac):
        if fr > area_fraction and start is None:
            start = i
        elif fr <= area_fraction and start is not None:
            if t[i - 1] - t[start] >= min_days:
                out.append(dict(start=float(t[start]), end=float(t[i - 1]), peak=float(anom[start:i].max()),
                                area=float(frac[start:i].max())))
            start = None
    if start is not None and t[-1] - t[start] >= min_days:
        out.append(dict(start=float(t[start]), end=float(t[-1]), peak=float(anom[start:].max()),
                        area=float(frac[start:].max())))
    return out


def find_collapses(ds, names, drop=0.5, window=730):
    """Population collapses: biomass falling below `drop` of its running peak within `window` days."""
    t = ds["time"].values
    B = ds["biomass"].sum("stage").values
    out = []
    for q, n in enumerate(names):
        b = B[:, q]
        peak, peak_t, done = b[0], t[0], False
        for i in range(1, len(t)):
            if b[i] > peak:
                peak, peak_t, done = b[i], t[i], False
            elif not done and b[i] < drop * peak and t[i] - peak_t <= window:
                out.append(dict(species=n, start=float(peak_t), time=float(t[i]),
                                ratio=float(b[i] / max(peak, 1e-300))))
                done = True
    return out


# ---------------------------------------------------------------- report
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("run")
    p.add_argument("observations", nargs="?")
    p.add_argument("--ensemble", help="ensemble.nc for reliability scores")
    p.add_argument("--heatwave-threshold", type=float, default=1.5)
    p.add_argument("--climatology", help="control run whose seasonal cycle defines the anomalies")
    p.add_argument("--collapse-drop", type=float, default=0.5)
    a = p.parse_args()
    run = Path(a.run)
    ds = xr.open_dataset(run / "series.nc", decode_times=False)
    names = ds.attrs["species"].split(",")
    obs = tomllib.loads(Path(a.observations).read_text()) if a.observations else {}
    lines = [f"# fishNET skill report: `{run}`", ""]
    store = {}

    point = obs.get("target", [])
    if point:
        lines += ["## Point targets", "", "| target | observed | sigma | model | bias | within sigma |",
                  "|---|---|---|---|---|---|"]
        for d in point:
            t, x = series_of(ds, d, names)
            lo, hi = d.get("day_range", [t.max() - 365, t.max()])
            mod = float(np.nanmean(x[(t >= lo) & (t <= hi)]))
            o, sig = float(d["value"]), float(d["sigma"])
            rel = mod / o - 1
            ok = abs(np.log(max(mod, 1e-300) / o)) <= sig if d.get("log", True) else abs(mod - o) <= sig
            lines.append(f"| {d['name']} | {o:.4g} | {sig:.3g} | {mod:.4g} | {rel:+.1%} | {'yes' if ok else 'no'} |")
            store[f"point_{d['name']}"] = [o, mod, rel]
        lines.append("")

    taylor = []
    ser = obs.get("series", [])
    if ser:
        lines += ["## Time series", "", "| series | n | bias | RMSE | centred RMS | corr | sd ratio | NSE |",
                  "|---|---|---|---|---|---|---|---|"]
        for d in ser:
            t, x = series_of(ds, d, names)
            od, ov = np.asarray(d["days"], float), np.asarray(d["values"], float)
            mod = np.interp(od, t, x, left=np.nan, right=np.nan) * 1.0
            if d.get("log", True):
                mod, ov = np.log(np.maximum(mod, 1e-300)), np.log(np.maximum(ov, 1e-300))
            mt = metrics(mod, ov)
            taylor.append((d["name"], mt))
            store[f"series_{d['name']}"] = [mt["bias"], mt["rmse"], mt["corr"], mt["sd_ratio"], mt["nse"]]
            lines.append(f"| {d['name']} | {mt['n']} | {mt['bias']:+.3g} | {mt['rmse']:.3g} | {mt['crms']:.3g} | "
                         f"{mt['corr']:.2f} | {mt['sd_ratio']:.2f} | {mt['nse']:.2f} |")
        lines.append("")

    spec = size_spectrum(run, ds, names)
    if spec:
        want = obs.get("size_spectrum", {}).get("slope")
        lines += ["## Community size spectrum", "",
                  f"- modelled slope **{spec['slope']:.2f}** (numbers per unit log mass), intercept "
                  f"{spec['intercept']:.2f}, from {spec['bins']} populated bins ({spec['source']})"]
        if want is not None:
            lines.append(f"- observed slope {want:.2f}, difference {spec['slope'] - float(want):+.2f}")
        store["size_spectrum"] = [spec["slope"], spec["intercept"]]
        lines.append("")

    fields = run / "fields.nc"
    if fields.exists():
        F = xr.open_dataset(fields, decode_times=False)
        ctrl = xr.open_dataset(Path(a.climatology) / "fields.nc", decode_times=False) if a.climatology else None
        hw = find_heatwaves(F, a.heatwave_threshold, control=ctrl)
        lines += ["## Events", "", f"- marine heatwaves (SST anomaly > {a.heatwave_threshold} degC): {len(hw)}"]
        if ctrl is None and float(F["time"].max() - F["time"].min()) < 730:
            lines.append("  - note: the run is shorter than two years, so the fitted seasonal cycle absorbs part of"
                         " any event; pass --climatology with a control run for a clean anomaly")
        for h in hw[:8]:
            lines.append(f"  - day {h['start']:.0f}-{h['end']:.0f}, peak +{h['peak']:.2f} degC over up to "
                         f"{h['area']:.0%} of the ocean")
        store["n_heatwaves"] = [len(hw)]
    else:
        lines += ["## Events", ""]
    col = find_collapses(ds, names, a.collapse_drop)
    lines.append(f"- collapses (biomass below {a.collapse_drop:.0%} of a recent peak): {len(col)}")
    for c in col[:8]:
        lines.append(f"  - {c['species']}: peak day {c['start']:.0f} to {c['ratio']:.0%} by day {c['time']:.0f}")
    store["n_collapses"] = [len(col)]
    for d in obs.get("event", []):                       # observed events to match against
        kind, when, tol = d.get("kind", "collapse"), float(d["day"]), float(d.get("tolerance_days", 180))
        cand = [c["time"] for c in col if c.get("species") == d.get("species")] if kind == "collapse" \
            else [h["start"] for h in (hw if fields.exists() else [])]
        hit = [c for c in cand if abs(c - when) <= tol]
        lines.append(f"- observed {kind} {d.get('species', '')} near day {when:.0f}: "
                     f"{'matched at day ' + f'{hit[0]:.0f}' if hit else 'not reproduced'}")
        store[f"event_{d.get('name', kind)}"] = [1.0 if hit else 0.0]
    lines.append("")

    if a.ensemble:
        E = xr.open_dataset(a.ensemble, decode_times=False)
        lines += ["## Ensemble reliability", "", "| target | members | observed | ensemble mean | CRPS | rank |",
                  "|---|---|---|---|---|---|"]
        for d in point:
            sp = names.index(d["species"]) if "species" in d else slice(None)
            x = E["biomass"][:, :, sp].sum("stage") if d["kind"] == "biomass" else E["plankton"][:, :, 1]
            vals = np.asarray(x[:, -1], float)
            o = float(d["value"])
            lines.append(f"| {d['name']} | {len(vals)} | {o:.4g} | {np.nanmean(vals):.4g} | "
                         f"{crps(vals, o):.4g} | {rank_histogram(vals, o)}/{len(vals)} |")
            store[f"crps_{d['name']}"] = [crps(vals, o), rank_histogram(vals, o)]
        lines.append("")

    out = run / "skill.md"
    out.write_text("\n".join(lines) + "\n")
    if taylor:
        taylor_diagram(taylor, run / "skill_taylor.png")
    sk = xr.Dataset({k: (f"{k}_value", np.asarray(v, float)) for k, v in store.items()})
    sk.to_netcdf(run / "skill.nc")
    print("\n".join(lines))
    print(f"written: {out}, {run / 'skill.nc'}" + (f", {run / 'skill_taylor.png'}" if taylor else ""))


def taylor_diagram(entries, path):
    """Correlation as angle, standard-deviation ratio as radius; the reference sits at (1, 0)."""
    fig = plt.figure(figsize=(6, 5.5))
    ax = fig.add_subplot(111, polar=True)
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    for r in (0.5, 1.0, 1.5):
        ax.plot(np.linspace(0, np.pi / 2, 100), np.full(100, r), color="0.85", lw=0.8)
    for name, mt in entries:
        if not np.isfinite(mt["corr"]) or not np.isfinite(mt["sd_ratio"]):
            continue
        ax.plot(np.arccos(np.clip(mt["corr"], -1, 1)), mt["sd_ratio"], "o", label=name)
    ax.plot(0, 1, "k*", markersize=13, label="observations")
    ax.set_thetagrids(np.degrees(np.arccos([0.0, 0.4, 0.7, 0.9, 0.95, 0.99])),
                      labels=["0", "0.4", "0.7", "0.9", "0.95", "0.99"])
    ax.set_title("Taylor diagram: correlation and variability")
    ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.1), fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
