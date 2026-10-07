#!/usr/bin/env python3
"""Figures and animations from a fishNET run directory.

    python plot.py runs/<name> [--no-anim]
"""
import sys
from pathlib import Path
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3, "font.size": 9})


def load(run):
    F, S, A = (xr.open_dataset(run / f, decode_times=False) for f in ("fields.nc", "series.nc", "agents.nc"))
    L = xr.open_dataset(run / "lifehist.nc", decode_times=False) if (run / "lifehist.nc").exists() else xr.Dataset()
    meta = {k: F.attrs.get(k, "").split(",") for k in ("species", "colors", "traits", "stages", "causes", "behaviors")}
    return F, S, A, L, meta


def mapax(ax, F, field, title, cmap="viridis", **kw):
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad("0.75")
    im = ax.pcolormesh(F.lon, F.lat, field, cmap=cm, shading="auto", **kw)
    ax.set_title(title)
    ax.set_aspect(1 / np.cos(np.radians(float(F.lat.mean()))))
    ax.grid(False)
    return im


def timeseries(F, S, M, out):
    sp, col, t = M["species"], M["colors"], S.time.values
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    B = S.biomass.sum("stage") / 1e6
    for q, n in enumerate(sp):
        ax[0, 0].semilogy(t, B[:, q], color=col[q], label=n)
        ax[0, 1].plot(t, S.agents[:, q], color=col[q])
        ax[0, 2].semilogy(t, np.maximum(S.eggs[:, q], 1), color=col[q])
        ax[1, 2].plot(t, S.max_gen[:, q], color=col[q])
    ax[0, 0].set(title="total biomass (t)", xlabel="day"); ax[0, 0].legend()
    ax[0, 1].set(title="super-individuals (agents)", xlabel="day")
    ax[0, 2].set(title="egg production (d$^{-1}$)", xlabel="day")
    ax[1, 2].set(title="max generation among agents", xlabel="day")
    for p, name, c in zip(range(5), ["N", "P", "Z", "K", "D"], ["k", "g", "orange", "magenta", "brown"]):
        ax[1, 0].semilogy(t, S.plankton[:, p], color=c, label=name)
    ax[1, 0].set(title="NPZKD inventories (mmol N)", xlabel="day"); ax[1, 0].legend()
    ax[1, 1].plot(t, S.budget_error, "k")
    ax[1, 1].set(title="relative N budget error", xlabel="day")
    ax2 = ax[1, 1].twinx()
    ax2.plot(t, S.ibm_cells, color="tab:blue", alpha=0.5)
    ax2.set_ylabel("IBM cells", color="tab:blue")
    fig.tight_layout(); fig.savefig(out / "timeseries.png"); plt.close(fig)


def stages(S, M, out):
    sp, st = M["species"], M["stages"]
    fig, ax = plt.subplots(1, len(sp), figsize=(4 * len(sp), 3.2), squeeze=False)
    for q, n in enumerate(sp):
        b = S.biomass[:, q].values / 1e6
        ax[0, q].stackplot(S.time, b.T, labels=st, colors=plt.cm.viridis(np.linspace(0.1, 0.9, 4)))
        ax[0, q].set(title=f"{n}: biomass by stage (t)", xlabel="day")
    ax[0, 0].legend(loc="upper left")
    fig.tight_layout(); fig.savefig(out / "stages.png"); plt.close(fig)


def traits(S, M, out):
    sp, col, tr = M["species"], M["colors"], M["traits"]
    fig, ax = plt.subplots(4, 4, figsize=(15, 9), sharex=True)
    print(len(tr))
    for k, name in enumerate(tr):
        a = ax.flat[k]
        for q, n in enumerate(sp):
            m, s = S.trait_mean[:, q, k].values, S.trait_sd[:, q, k].values
            a.plot(S.time, m - m[0], color=col[q], label=n)
            a.fill_between(S.time, m - m[0] - s, m - m[0] + s, color=col[q], alpha=0.12)
        a.axhline(0, color="k", lw=0.5)
        a.set_title(f"{name}: change in mean")
    ax.flat[0].legend(); [a.set_xlabel("day") for a in ax[-1]]
    fig.suptitle("Evolution of population-mean breeding values (shading: +/- genetic SD)")
    fig.tight_layout(); fig.savefig(out / "traits.png"); plt.close(fig)


def behaviour(S, M, out):
    if "behavior" not in S or float(S.behavior.sum()) == 0:
        return
    sp, names = M["species"], M["behaviors"]
    fig, ax = plt.subplots(1, len(sp), figsize=(4 * len(sp), 3.4), squeeze=False)
    for q, n in enumerate(sp):
        ax[0, q].stackplot(S.time, S.behavior[:, q].values.T, labels=names, colors=plt.cm.tab10(np.arange(len(names))))
        ax[0, q].set(title=f"{n}: behaviour budget (agents)", xlabel="day", ylim=(0, 1))
    ax[0, 0].legend(loc="lower left", fontsize=7)
    fig.tight_layout(); fig.savefig(out / "behaviour.png"); plt.close(fig)


def biogeochem(F, S, M, out):
    """Oxygen, the metabolic index of each species, and carbon export split by pathway."""
    if "oxygen" not in S:
        return
    sp, col = M["species"], M["colors"]
    fig, ax = plt.subplots(1, 4, figsize=(18, 3.8))
    ax[0].plot(S.time, S.o2_min, color="crimson", label="basin minimum")
    ax[0].plot(S.time, S.oxygen / float(S.oxygen[0]) * float(S.o2_min[0]), color="0.5", lw=1, label="inventory (scaled)")
    ax[0].set(title="dissolved oxygen (mmol m$^{-3}$)", xlabel="day"); ax[0].legend(fontsize=7)
    if "O2" in F:
        o = F.O2.isel(time=-1).mean("lon")
        im = ax[1].pcolormesh(F.lat, F.depth, o, cmap="viridis", shading="auto")
        ax[1].invert_yaxis(); plt.colorbar(im, ax=ax[1]); ax[1].set(title="zonal-mean O$_2$, final day", xlabel="lat", ylabel="depth (m)")
    for q, n in enumerate(sp):
        ax[2].plot(S.time, S.metabolic_index[:, q], color=col[q], label=n)
    ax[2].axhline(1, color="k", ls="--", lw=0.8)
    ax[2].set(title="metabolic index (1 = suffocation)", xlabel="day", yscale="log"); ax[2].legend(fontsize=7)
    A = float(S.attrs.get("ocean_area", 0)) or 1.0
    gc = lambda x: x * 12.011 / 1000 / A * 365
    ax[3].stackplot(S.time, gc(S.export_sinking), gc(S.export_fish), gc(S.export_respired),
                    labels=["sinking detritus", "fish faeces/carcasses", "fish respiration"], colors=["0.4", "#c46", "#48c"])
    ax[3].set(title="carbon export below the reference depth", xlabel="day",
              ylabel="g C m$^{-2}$ yr$^{-1}$" if A > 1 else "mmol C d$^{-1}$")
    ax[3].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out / "biogeochem.png"); plt.close(fig)


def losses(S, M, out):
    sp, causes = M["species"], M["causes"]
    dt = np.gradient(S.time.values)
    tot = (S.losses * dt[:, None, None]).sum("time").values
    use = [c for c in range(len(causes)) if tot[:, c].sum() > 0]
    frac = tot[:, use] / tot[:, use].sum(1, keepdims=True)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    left = np.zeros(len(sp))
    for c, f in zip(use, frac.T):
        ax[0].barh(sp, f, left=left, label=causes[c]); left += f
    ax[0].set(title="share of biomass lost by cause", xlim=(0, 1)); ax[0].legend(fontsize=8)
    for q, n in enumerate(sp):
        ax[1].plot(S.time, S.catch[:, q] / 1e6, color=M["colors"][q], label=n)
    ax[1].set(title="catch (t d$^{-1}$)", xlabel="day"); ax[1].legend()
    fig.tight_layout(); fig.savefig(out / "losses.png"); plt.close(fig)


def arrows(ax, F, f, k=2):
    """Surface-current arrows (every k-th cell) so eddies and gyres are visible."""
    u, v = (x.isel(depth=0) if "depth" in x.dims else x for x in (f.u, f.v))
    return ax.quiver(F.lon[::k], F.lat[::k], u[::k, ::k], v[::k, ::k], color="k", alpha=0.6, scale=12, width=0.003)


def maps(F, A, M, out, it=-1):
    sp, col = M["species"], M["colors"]
    ns = len(sp)
    fig, ax = plt.subplots(2, max(3, ns), figsize=(4.2 * max(3, ns), 7.5))
    f = F.isel(time=it)
    plt.colorbar(mapax(ax[0, 0], F, f.temp[0], "SST (°C) + surface currents", "RdYlBu_r"), ax=ax[0, 0], shrink=0.7)
    arrows(ax[0, 0], F, f)
    plt.colorbar(mapax(ax[0, 1], F, np.log10(f.P[0]), "log10 surface P (mmol N m$^{-3}$)", "YlGn"), ax=ax[0, 1], shrink=0.7)
    plt.colorbar(mapax(ax[0, 2], F, np.log10(f.Z[:2].mean("depth")), "log10 Z 0-75 m", "PuBu"), ax=ax[0, 2], shrink=0.7)
    for a in ax[0, 3:]:
        a.axis("off")
    snap = A.where(A.time == A.time.max(), drop=True) if A.sizes.get("rec", 0) else None
    for q, n in enumerate(sp):
        b = np.log10(np.maximum(f.fish_biomass[q].sum("stage"), 1e-3))
        plt.colorbar(mapax(ax[1, q], F, b, f"log10 {n} biomass (g m$^{{-2}}$)", "magma", vmin=-2), ax=ax[1, q], shrink=0.7)
        ax[1, q].contour(F.lon, F.lat, f.ibm, [0.5], colors="c", linewidths=0.6)
        if snap is not None:
            s = snap.species == q
            ax[1, q].scatter(snap.lon[s], snap.lat[s], s=1, c="w", alpha=0.4)
    fig.suptitle(f"day {float(f.time):.1f} (cyan = IBM cells, dots = agents at last snapshot)")
    fig.tight_layout(); fig.savefig(out / "maps.png"); plt.close(fig)


def trait_maps(F, M, out, trait="T_opt"):
    sp = M["species"]
    k = M["traits"].index(trait)
    fig, ax = plt.subplots(2, len(sp), figsize=(4.2 * len(sp), 7), squeeze=False)
    for row, it in enumerate((0, -1)):
        f = F.isel(time=it)
        for q, n in enumerate(sp):
            tm = f.trait_mean[q, k].where(f.fish_biomass[q].sum("stage") > 1e-2)
            plt.colorbar(mapax(ax[row, q], F, tm, f"{n} {trait} day {float(f.time):.0f}", "coolwarm"), ax=ax[row, q], shrink=0.7)
    fig.suptitle(f"Local mean breeding value for {trait} (masked where biomass < 0.01 g m$^{{-2}}$)")
    fig.tight_layout(); fig.savefig(out / f"trait_map_{trait}.png"); plt.close(fig)


def hovmoller(F, M, out):
    sp = M["species"]
    fig, ax = plt.subplots(1, len(sp) + 1, figsize=(4 * (len(sp) + 1), 3.6))
    sst = F.temp[:, 0].mean("lon")
    ax[0].pcolormesh(F.time, F.lat, sst.T, cmap="RdYlBu_r", shading="auto"); ax[0].set(title="zonal-mean SST", xlabel="day")
    for q, n in enumerate(sp):
        b = np.log10(np.maximum(F.fish_biomass[:, q].sum("stage").mean("lon"), 1e-3))
        ax[q + 1].pcolormesh(F.time, F.lat, b.T, cmap="magma", vmin=-2, shading="auto")
        ax[q + 1].set(title=f"{n}: zonal-mean log10 biomass", xlabel="day")
    ax[0].set_ylabel("latitude")
    fig.tight_layout(); fig.savefig(out / "hovmoller.png"); plt.close(fig)


def life_history(A, L, M, out):
    sp, col, causes = M["species"], M["colors"], M["causes"]
    fig, ax = plt.subplots(1, 4, figsize=(17, 3.8))
    if L.sizes.get("agent", 0):
        born = L.origin == 2
        for q, n in enumerate(sp):
            s = born & (L.species == q)
            if s.sum():
                ax[0].hist((L.t_end - L.tb)[s], bins=40, histtype="step", color=col[q], label=n)
            ax[1].bar(np.arange(len(causes)) + 0.2 * q - 0.3, [(s & (L.cause == c)).sum() for c in range(len(causes))],
                      0.2, color=col[q])
            g = L.gen[L.species == q]
            ax[2].hist(g, bins=np.arange(g.max() + 2) - 0.5 if len(g) else 1, histtype="step", color=col[q])
        ax[0].set(title="agent lifetimes, born in model (days)", yscale="log"); ax[0].legend()
        ax[1].set_xticks(range(len(causes)), causes, rotation=40, ha="right")
        ax[1].set(title="fate of agents born in model", yscale="log")
        ax[2].set(title="generation of all agents", yscale="log")
    if A.sizes.get("rec", 0):
        idx = np.random.default_rng(0).choice(A.sizes["rec"], min(20000, A.sizes["rec"]), replace=False)
        a = A.isel(rec=np.sort(idx))
        for q, n in enumerate(sp):
            s = a.species == q
            ax[3].scatter(a.age[s] / 365, a.L[s], s=2, color=col[q], alpha=0.3, label=n)
        ax[3].set(title="length at age (agent snapshots)", xlabel="age (yr)", ylabel="L (cm)", yscale="log")
    fig.tight_layout(); fig.savefig(out / "life_history.png"); plt.close(fig)


def tracks(F, A, M, out, n=40):
    if not A.sizes.get("rec", 0):
        return
    ids, cnt = np.unique(A.id.values, return_counts=True)
    keep = ids[np.argsort(cnt)[-n:]]
    fig, ax = plt.subplots(figsize=(8, 6))
    mapax(ax, F, F.temp[-1, 0], "longest-tracked agents over SST", "RdYlBu_r")
    for i in keep:
        s = A.where(A.id == i, drop=True).sortby("time")
        q = int(s.species[0])
        ax.plot(s.lon, s.lat, "-", color=M["colors"][q], lw=1)
        ax.plot(s.lon[-1], s.lat[-1], "o", color=M["colors"][q], ms=3)
    for q, nm in enumerate(M["species"]):
        ax.plot([], [], color=M["colors"][q], label=nm)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / "tracks.png"); plt.close(fig)


def depth_use(A, M, out):
    if not A.sizes.get("rec", 0):
        return
    fig, ax = plt.subplots(1, len(M["species"]), figsize=(3.4 * len(M["species"]), 3.2), sharey=True)
    for q, n in enumerate(M["species"]):
        s = (A.species == q).values
        for st, name in enumerate(M["stages"]):
            w = s & (A.stage == st).values
            if w.any():
                d, c = np.unique(A.depth.values[w], return_counts=True)
                ax[q].plot(c / c.sum(), d, "o-", label=name)
        ax[q].set(title=n, xlabel="fraction of agents"); ax[q].invert_yaxis() if q == 0 else None
    ax[0].set_ylabel("depth (m)"); ax[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out / "depth_use.png"); plt.close(fig)


def save_anim(fig, update, frames, path):
    writer = animation.FFMpegWriter(fps=8, bitrate=2400) if animation.writers.is_available("ffmpeg") \
        else animation.PillowWriter(fps=8)
    path = path.with_suffix(".mp4" if isinstance(writer, animation.FFMpegWriter) else ".gif")
    animation.FuncAnimation(fig, update, frames=frames).save(path, writer=writer)
    plt.close(fig)
    print("  wrote", path)


def anim_fields(F, M, out):
    sp = M["species"]
    fig, ax = plt.subplots(2, max(2, (len(sp) + 3) // 2 + 1), figsize=(15, 7.5))
    ax = ax.flat
    ims = [mapax(ax[0], F, F.temp[0, 0], "SST", "RdYlBu_r", vmin=float(F.temp[:, 0].min()), vmax=float(F.temp[:, 0].max())),
           mapax(ax[1], F, np.log10(F.P[0, 0]), "log10 surface P", "YlGn", vmin=-2, vmax=0.7)]
    for q, n in enumerate(sp):
        ims.append(mapax(ax[q + 2], F, np.log10(np.maximum(F.fish_biomass[0, q].sum("stage"), 1e-3)),
                         f"log10 {n} (g m$^{{-2}}$)", "magma", vmin=-2, vmax=2))
    for a in ax[len(sp) + 2:]:
        a.axis("off")
    for a, im in zip(ax, ims):
        plt.colorbar(im, ax=a, shrink=0.7)
    title = fig.suptitle("")
    cont = []
    q0 = arrows(ax[0], F, F.isel(time=0))

    def update(i):
        f = F.isel(time=i)
        q0.set_UVC(*(x.isel(depth=0)[::2, ::2] for x in (f.u, f.v)))
        ims[0].set_array(f.temp[0].values.ravel())
        ims[1].set_array(np.log10(f.P[0]).values.ravel())
        for q in range(len(sp)):
            ims[q + 2].set_array(np.log10(np.maximum(f.fish_biomass[q].sum("stage"), 1e-3)).values.ravel())
        for c in cont:
            c.remove()
        cont.clear()
        for q in range(len(sp)):
            cont.append(ax[q + 2].contour(F.lon, F.lat, f.ibm, [0.5], colors="c", linewidths=0.5))
        title.set_text(f"day {float(f.time):.1f}   (cyan outline = IBM cells)")
    fig.tight_layout()
    save_anim(fig, update, F.sizes["time"], out / "anim_fields")


def anim_agents(F, A, M, out):
    if not A.sizes.get("rec", 0):
        return
    times = np.unique(A.time.values)
    ft = F.time.values
    fig, ax = plt.subplots(figsize=(9, 6.5))
    im = mapax(ax, F, np.log10(F.P[0, 0]), "", "Greens", vmin=-2, vmax=0.7)
    plt.colorbar(im, ax=ax, shrink=0.7, label="log10 surface P")
    sc = [ax.scatter([], [], s=4, color=M["colors"][q], label=n, alpha=0.7) for q, n in enumerate(M["species"])]
    ax.legend(loc="lower left", fontsize=8, markerscale=3)
    qv = arrows(ax, F, F.isel(time=0))

    def update(k):
        a = A.where(A.time == times[k], drop=True)
        f = F.isel(time=int(np.argmin(np.abs(ft - times[k]))))
        im.set_array(np.log10(f.P[0]).values.ravel())
        qv.set_UVC(*(x.isel(depth=0)[::2, ::2] for x in (f.u, f.v)))
        for q, s in enumerate(sc):
            w = a.species == q
            s.set_offsets(np.c_[a.lon[w], a.lat[w]] if w.any() else np.empty((0, 2)))
            s.set_sizes(1 + 0.4 * np.log10(np.maximum(a.n[w].values, 1)))
        ax.set_title(f"agents on day {times[k]:.1f}  (n = {a.sizes['rec']:,d}; size ~ log fish represented)")
    save_anim(fig, update, len(times), out / "anim_agents")


def ensemble_plots(run, out):
    E = xr.open_dataset(run / "ensemble.nc", decode_times=False)
    sp, col = E.attrs["species"].split(","), E.attrs["colors"].split(",")
    ok = np.isfinite(E.biomass.values).all(axis=tuple(range(1, E.biomass.ndim)))
    E = E.isel(member=np.nonzero(ok)[0])
    B = E.biomass.sum("stage") / 1e6
    params = [E[v] for v in E.data_vars if v.startswith("param_") and E[v].attrs["path"] != "run.seed"]
    print(f"ensemble: {E.sizes['member']} members, {E.sizes['time']} steps, parameters: "
          + ", ".join(p.attrs["path"] for p in params))
    fig, ax = plt.subplots(1, len(sp), figsize=(4.2 * len(sp), 3.6), squeeze=False)
    for q, n in enumerate(sp):
        a, b = ax[0, q], B[:, :, q]
        a.plot(E.time, b.T, color=col[q], alpha=0.35, lw=0.8)
        a.fill_between(E.time, b.quantile(0.1, "member"), b.quantile(0.9, "member"), color=col[q], alpha=0.2)
        a.plot(E.time, b.median("member"), color="k", lw=1.8)
        a.set(yscale="log", title=f"{n}: biomass (t), median and 10-90%", xlabel="day")
    fig.tight_layout(); fig.savefig(out / "ensemble_spread.png"); plt.close(fig)
    y = np.log10(B.isel(time=slice(-max(1, E.sizes["time"] // 4), None)).mean("time").values + 1e-9)
    code = lambda p: p.values.astype(float) if p.dtype.kind == "f" else np.unique(p.values, return_inverse=True)[1].astype(float)
    X = [(p.attrs["path"], code(p)) for p in params if np.nanstd(code(p)) > 0]
    rank = lambda x: np.argsort(np.argsort(x))
    R = np.array([[np.corrcoef(rank(x), rank(y[:, q]))[0, 1] for q in range(len(sp))] for _, x in X]).reshape(len(X), len(sp))
    fig, ax = plt.subplots(figsize=(1.6 * len(sp) + 4, 0.45 * len(X) + 1.8))
    im = ax.imshow(R, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set(xticks=range(len(sp)), xticklabels=sp, yticks=range(len(X)), yticklabels=[n for n, _ in X],
           title="Spearman correlation: parameter vs final-quarter log biomass")
    for (i, j), r in np.ndenumerate(R):
        ax.text(j, i, f"{r:+.2f}", ha="center", va="center", fontsize=8)
    plt.colorbar(im, ax=ax); ax.grid(False)
    fig.tight_layout(); fig.savefig(out / "ensemble_sensitivity.png"); plt.close(fig)
    fig, ax = plt.subplots(len(X), len(sp), figsize=(3.2 * len(sp), 2.5 * len(X)), squeeze=False)
    for i, (n, x) in enumerate(X):
        for q in range(len(sp)):
            ax[i, q].scatter(x, 10 ** y[:, q], color=col[q], s=18)
            ax[i, q].set(yscale="log", xlabel=n if i == len(X) - 1 else "", ylabel=f"{sp[q]} (t)" if True else "")
            if i == 0:
                ax[i, q].set_title(sp[q])
        ax[i, 0].set_ylabel(f"{n}\nbiomass (t)", fontsize=8)
    fig.tight_layout(); fig.savefig(out / "ensemble_scatter.png"); plt.close(fig)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    run = Path(sys.argv[1])
    out = run / "figures"
    out.mkdir(exist_ok=True)
    if (run / "ensemble.nc").exists():
        ensemble_plots(run, out)
        print("figures in", out, "(plot single members with: python plot.py", run / "member_000)")
        return
    F, S, A, L, M = load(run)
    print(f"plotting {run}: {F.sizes['time']} field frames, {A.sizes.get('rec', 0):,d} agent records, "
          f"{L.sizes.get('agent', 0):,d} life histories")
    for name, fn in (("timeseries", lambda: timeseries(F, S, M, out)), ("stages", lambda: stages(S, M, out)),
                     ("traits", lambda: traits(S, M, out)), ("behaviour", lambda: behaviour(S, M, out)),
                     ("biogeochem", lambda: biogeochem(F, S, M, out)), ("losses", lambda: losses(S, M, out)),
                     ("maps", lambda: maps(F, A, M, out)), ("trait maps", lambda: trait_maps(F, M, out)),
                     ("hovmoller", lambda: hovmoller(F, M, out)), ("life history", lambda: life_history(A, L, M, out)),
                     ("tracks", lambda: tracks(F, A, M, out)), ("depth use", lambda: depth_use(A, M, out))):
        print(" ", name)
        fn()
    if "--no-anim" not in sys.argv:
        anim_fields(F, M, out)
        anim_agents(F, A, M, out)
    print("figures in", out)


if __name__ == "__main__":
    main()
