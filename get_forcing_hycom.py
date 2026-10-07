#!/usr/bin/env python3
"""Download HYCOM + NCODA ocean reanalysis (GOFS 3.1) as fishNET netCDF forcing. No account needed.

HYCOM is served over OPeNDAP by the Naval Research Laboratory, free and without registration, which
makes it the quickest route to realistic forcing. GLORYS (`get_forcing.py`) is the better product for
published work, but needs a Copernicus Marine account.

    python get_forcing_hycom.py --preset north_atlantic --start 2010-01-01 --end 2011-12-31
    python get_forcing_hycom.py --lon -70 -10 --lat 20 60 --start 2012-01-01 --end 2012-12-31 --every 1

Writes <out>/forcing.nc (water_temp, water_u, water_v on time, depth, lat, lon; fishNET auto-detects
these names) and prints a namelist fragment. The file is written incrementally, so an interrupted
download leaves a usable file covering the dates fetched so far.

Coverage: GLBv0.08/expt_53.X, 1994-2015, 1/12 degree, 40 levels, 3-hourly (sampled with --every).
Later years live in other experiments; pass --url to point at one.
"""
import argparse, sys, time
from pathlib import Path

import numpy as np

URL = "https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/{year}"
PRESETS = dict(north_sea=(-4, 12, 51, 62), benguela=(8, 20, -35, -15), california=(-130, -115, 30, 45),
               nw_atlantic=(-75, -50, 35, 50), peru=(-85, -70, -20, -5), north_atlantic=(-70, -10, 20, 60))
VARS = ("water_temp", "water_u", "water_v")


def open_year(url, year, tries=4):
    import xarray as xr
    for k in range(tries):
        try:
            return xr.open_dataset(url.format(year=year), drop_variables=["tau"])
        except Exception as e:
            print(f"  open {year}: {type(e).__name__}, retry {k + 1}/{tries}", flush=True)
            time.sleep(5 * (k + 1))
    sys.exit(f"cannot open {url.format(year=year)}")


def fetch(da, tries=4):
    for k in range(tries):
        try:
            return da.values
        except Exception as e:
            print(f"  fetch: {type(e).__name__}: {str(e)[:80]}, retry {k + 1}/{tries}", flush=True)
            time.sleep(5 * (k + 1))
    raise RuntimeError("OPeNDAP fetch failed")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--preset", choices=sorted(PRESETS))
    p.add_argument("--lon", nargs=2, type=float, metavar=("W", "E"))
    p.add_argument("--lat", nargs=2, type=float, metavar=("S", "N"))
    p.add_argument("--start", default="2010-01-01")
    p.add_argument("--end", default="2010-12-31")
    p.add_argument("--every", type=float, default=1.0, help="days between kept snapshots (0.125 = every 3 h)")
    p.add_argument("--stride", type=int, default=8, help="keep every Nth grid point (8 = ~0.64 deg)")
    p.add_argument("--depth", type=float, default=600.0, help="deepest level to keep (m)")
    p.add_argument("--levels", type=int, default=12, help="keep about this many levels, spread by log depth "
                                                          "(0 = every level HYCOM has; fewer = much faster)")
    p.add_argument("--hour", type=int, default=12, help="hour of day to sample")
    p.add_argument("--chunk", type=int, default=8, help="time steps per OPeNDAP request")
    p.add_argument("--out", default="forcing_hycom")
    p.add_argument("--url", default=URL, help="OPeNDAP template with {year}")
    a = p.parse_args()
    if a.preset:
        a.lon, a.lat = PRESETS[a.preset][:2], PRESETS[a.preset][2:]
    if not (a.lon and a.lat):
        p.error("give --preset or both --lon and --lat")

    import netCDF4 as nc
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "forcing.nc"
    t0, t1 = np.datetime64(a.start), np.datetime64(a.end) + np.timedelta64(1, "D")
    years = range(int(str(t0)[:4]), int(str(t1 - np.timedelta64(1, "s"))[:4]) + 1)

    ds0 = open_year(a.url, years[0])
    lat = ds0.lat.sel(lat=slice(*a.lat)).values[::a.stride]
    lon = ds0.lon.sel(lon=slice(*a.lon)).values[::a.stride]
    zall = ds0.depth.sel(depth=slice(0, a.depth)).values
    if a.levels and a.levels < len(zall):                  # thin by log depth: dense near the surface
        want = np.expm1(np.linspace(0, np.log1p(zall[-1]), a.levels))
        kz = sorted({int(np.argmin(np.abs(zall - w))) for w in want} | {0, len(zall) - 1})
    else:
        kz = list(range(len(zall)))
    dep = zall[kz]
    if not (len(lat) and len(lon)):
        sys.exit("the box selects no grid points; check --lon/--lat (HYCOM lon is -180..180 here)")
    print(f"box lon {lon[0]:.2f}..{lon[-1]:.2f} ({len(lon)}), lat {lat[0]:.2f}..{lat[-1]:.2f} ({len(lat)}), "
          f"{len(dep)} levels to {dep[-1]:.0f} m, {a.start}..{a.end} every {a.every} d")
    per = len(dep) * len(lat) * len(lon) * 4 * 3 / 1e6
    print(f"~{per:.2f} MB per snapshot (3 variables)")

    f = nc.Dataset(path, "w")
    f.setncatts(dict(source=a.url.format(year="<year>"), title="HYCOM + NCODA GOFS 3.1 reanalysis subset",
                     note="fetched by get_forcing_hycom.py for fishNET"))
    f.createDimension("time", None)
    for k, v in (("depth", dep), ("lat", lat), ("lon", lon)):
        f.createDimension(k, len(v))
        f.createVariable(k, "f8", (k,))[:] = v
    f["depth"].units, f["lat"].units, f["lon"].units = "m", "degrees_north", "degrees_east"
    tv = f.createVariable("time", "f8", ("time",))
    tv.units, tv.calendar = "days since 2000-01-01 00:00:00", "gregorian"
    for v, u in zip(VARS, ("degC", "m s-1", "m s-1")):
        f.createVariable(v, "f4", ("time", "depth", "lat", "lon"), zlib=True, complevel=1,
                         fill_value=np.float32(np.nan)).units = u

    n, t_start = 0, time.time()
    for year in years:
        ds = open_year(a.url, year)
        want = np.arange(max(t0, np.datetime64(f"{year}-01-01")), min(t1, np.datetime64(f"{year + 1}-01-01")),
                         np.timedelta64(int(round(a.every * 24)), "h")) + np.timedelta64(a.hour, "h")
        have = ds.time.values
        idx = sorted({int(np.argmin(np.abs(have - w))) for w in want if abs(have - w).min() < np.timedelta64(6, "h")})
        print(f"{year}: {len(idx)} snapshots of {len(have)} available", flush=True)
        sel = dict(lat=slice(*a.lat), lon=slice(*a.lon))
        for c0 in range(0, len(idx), a.chunk):
            part = idx[c0:c0 + a.chunk]
            for v in VARS:
                da = ds[v].isel(time=part, depth=kz).sel(lat=sel["lat"], lon=sel["lon"])[:, :, ::a.stride, ::a.stride]
                f[v][n:n + len(part)] = fetch(da)
            tv[n:n + len(part)] = (have[part] - np.datetime64("2000-01-01")) / np.timedelta64(1, "D")
            n += len(part)
            f.sync()
            done = (time.time() - t_start)
            print(f"  {n} snapshots, {path.stat().st_size / 1e6:6.1f} MB, {done / 60:5.1f} min elapsed, "
                  f"{str(have[part[-1]])[:13]}", flush=True)
        ds.close()
    f.close()
    describe(path, dep, lat, lon)


def describe(path, dep, lat, lon):
    import netCDF4 as nc
    f = nc.Dataset(path)
    f.set_auto_mask(False)
    t = f["time"][:]
    T = f["water_temp"][0]
    edges = [z for z in (0, 25, 75, 150, 300, 600, 1200, 2500, 4000) if z <= dep[-1]]
    print(f"\n{path}: {f['water_temp'].shape} (time, depth, lat, lon), {path.stat().st_size / 1e6:.1f} MB")
    print(f"  temp {np.nanmin(T):.2f} .. {np.nanmax(T):.2f} degC, {100 * np.isnan(T).mean():.0f}% missing (land/below bed)")
    start = (np.datetime64("2000-01-01") + np.timedelta64(int(round(t[0] * 24)), "h")).astype(str)[:16]
    print(f"""
[run]
start = "{start}"

[grid]
nx = {max(8, int((lon[-1] - lon[0]) / 1.0))}
ny = {max(8, int((lat[-1] - lat[0]) / 1.0))}
lon = [{lon[0]:.2f}, {lon[-1]:.2f}]
lat = [{lat[0]:.2f}, {lat[-1]:.2f}]
depth_edges = [{", ".join(f"{x:.0f}" for x in edges)}]

[ocean]
source = "netcdf"
file = "{path}"
""")
    f.close()


if __name__ == "__main__":
    main()
