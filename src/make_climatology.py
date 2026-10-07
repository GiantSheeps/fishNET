#!/usr/bin/env python3
"""Build a smooth, loopable annual climatology from a multi-year monthly forcing file.

Monthly means are averaged by calendar month, then a mean + annual + semiannual harmonic is
fitted per cell and evaluated on a finer time axis, so linear interpolation between the stored
steps is visually smooth and the cycle closes on itself.

    python make_climatology.py --in forcing/forcing.nc --out forcing_clim/forcing.nc
    python make_climatology.py --in forcing/forcing.nc --slice-year 2001 --out forcing_1yr/forcing.nc

fishNET wraps netcdf forcing time modulo the file's span, so a one-year file loops forever.
"""
import argparse
from pathlib import Path

import numpy as np
import xarray as xr

YEAR = 365.0
MONTH_MID = np.array([15.5, 45.0, 74.5, 105.0, 135.5, 166.0,      # day-of-year centre of each
                      196.5, 227.5, 258.0, 288.5, 319.0, 349.5])  # month, the mean's true phase
VARS = ("thetao", "uo", "vo")


def design(t, nharm):
    """[1, cos wt, sin wt, cos 2wt, sin 2wt, ...] for days-of-year t."""
    w = 2 * np.pi / YEAR
    cols = [np.ones_like(t)]
    for k in range(1, nharm + 1):
        cols += [np.cos(k * w * t), np.sin(k * w * t)]
    return np.stack(cols, -1)


def fit_evaluate(clim, t_out, nharm):
    """clim (12, ...) monthly means -> (len(t_out), ...) harmonic reconstruction, land kept NaN."""
    shape = clim.shape[1:]
    Y = clim.reshape(12, -1)
    land = ~np.isfinite(Y).all(0)              # the mask is fixed in time (land, below the bed)
    coef = np.linalg.lstsq(design(MONTH_MID, nharm), np.nan_to_num(Y), rcond=None)[0]
    out = design(t_out, nharm) @ coef
    out[:, land] = np.nan
    return out.reshape((len(t_out),) + shape)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in", dest="src", default="forcing/forcing.nc")
    p.add_argument("--out", default="forcing_clim/forcing.nc")
    p.add_argument("--steps", type=int, default=73, help="time steps in the output year (73 = pentads)")
    p.add_argument("--harmonics", type=int, default=2, help="annual + semiannual by default")
    p.add_argument("--slice-year", type=int, help="write this calendar year's monthly means instead")
    p.add_argument("--coarsen", type=int, default=1, help="average NxN native cells (1 = native 1/12 deg)")
    p.add_argument("--year", type=int, default=2000, help="calendar year to stamp the output with")
    p.add_argument("--complevel", type=int, default=4)
    a = p.parse_args()

    ds = xr.open_dataset(a.src)
    if a.coarsen > 1:
        ds = ds.coarsen(longitude=a.coarsen, latitude=a.coarsen, boundary="trim").mean(skipna=True)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    if a.slice_year:
        sub = ds.sel(time=str(a.slice_year))
        print(f"slice {a.slice_year}: {sub.sizes['time']} steps")
        # Re-stamp month-start times onto the month centre so the seasonal phase is right. The centres
        # are spaced uniformly (365/12) rather than at each month's true midpoint: fishNET infers the
        # cycle length from the first interval, so uneven steps would make the loop drift ~1.5 d/year.
        cen = (sub.time.dt.month.values - 0.5) * YEAR / 12
        mid = np.datetime64(f"{a.slice_year}-01-01") + (cen * 86400).astype("timedelta64[s]")
        sub = sub.assign_coords(time=mid)
        write(sub, out, a)
        return

    t_out = (np.arange(a.steps) + 0.5) * YEAR / a.steps
    stamp = np.datetime64(f"{a.year}-01-01") + (t_out * 86400).astype("timedelta64[s]")
    print(f"climatology from {ds.sizes['time']} monthly steps "
          f"({str(ds.time.values[0])[:7]}..{str(ds.time.values[-1])[:7]}) "
          f"-> {a.steps} steps, {a.harmonics} harmonics")
    counts = np.bincount(ds.time.dt.month.values - 1, minlength=12)
    print("  samples per month:", " ".join(f"{m}:{c}" for m, c in zip("JFMAMJJASOND", counts)))

    res = xr.Dataset(coords=dict(time=stamp, depth=ds.depth, latitude=ds.latitude, longitude=ds.longitude))
    for v in VARS:
        clim = ds[v].groupby("time.month").mean("time").transpose("month", "depth", "latitude", "longitude")
        vals = fit_evaluate(clim.values.astype("float64"), t_out, a.harmonics).astype("float32")
        res[v] = (("time", "depth", "latitude", "longitude"), vals)
        res[v].attrs = dict(ds[v].attrs)
        rms = np.sqrt(np.nanmean((clim.values - fit_evaluate(clim.values.astype("float64"), MONTH_MID, a.harmonics)) ** 2))
        print(f"  {v}: {vals.shape}, fit residual rms {rms:.4f} {ds[v].attrs.get('units', '')}")
        del clim, vals
    write(res, out, a)


def write(res, out, a):
    enc = {v: dict(zlib=True, complevel=a.complevel, dtype="float32") for v in VARS if v in res}
    res.attrs["note"] = "annual climatology for cyclic forcing; time is a nominal year and loops"
    res.to_netcdf(out, encoding=enc)
    print(f"\nwrote {out} ({out.stat().st_size / 1e6:.0f} MB, {res.sizes['time']} steps)")
    print(f"  time {str(res.time.values[0])[:13]} .. {str(res.time.values[-1])[:13]}")


if __name__ == "__main__":
    main()
