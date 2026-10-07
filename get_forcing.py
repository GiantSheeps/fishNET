#!/usr/bin/env python3
"""Download GLORYS/CMEMS ocean reanalysis and write forcing fishNET can read.

Needs a free Copernicus Marine account (https://data.marine.copernicus.eu) and the toolbox:

    pip install copernicusmarine
    copernicusmarine login                    # stores credentials once

Then, for example:

    python get_forcing.py --lon -70 -10 --lat 20 60 --start 2015-01-01 --end 2019-12-31
    python get_forcing.py --preset north_sea --start 2000-01-01 --end 2020-12-31 --daily

It writes forcing.nc (thetao, uo, vo) and bathymetry.nc (deptho) into --out, plus a namelist
fragment to paste into your [ocean] section. Point the run at them with:

    [ocean]
    source = "netcdf"
    file = "forcing/forcing.nc"
    bathymetry = "forcing/bathymetry.nc"
"""
import argparse, sys, textwrap
from pathlib import Path

# Monthly and daily physics reanalysis (1993-present) and the matching static file.
MONTHLY = "cmems_mod_glo_phy_my_0.083deg_P1M-m"
DAILY = "cmems_mod_glo_phy_my_0.083deg_P1D-m"
STATIC = "cmems_mod_glo_phy_my_0.083deg_static"
PRESETS = dict(north_sea=(-4, 12, 51, 62), benguela=(8, 20, -35, -15), california=(-130, -115, 30, 45),
               nw_atlantic=(-75, -50, 35, 50), peru=(-85, -70, -20, -5), north_atlantic=(-70, -10, 20, 60))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--preset", choices=sorted(PRESETS), help="named region instead of --lon/--lat")
    p.add_argument("--lon", nargs=2, type=float, metavar=("W", "E"))
    p.add_argument("--lat", nargs=2, type=float, metavar=("S", "N"))
    p.add_argument("--start", default="2010-01-01")
    p.add_argument("--end", default="2019-12-31")
    p.add_argument("--depth", type=float, default=1000.0, help="deepest level to fetch (m)")
    p.add_argument("--daily", action="store_true", help="daily means instead of monthly (much larger)")
    p.add_argument("--out", default="forcing")
    p.add_argument("--dataset", help="override the dataset id")
    a = p.parse_args()
    if a.preset:
        a.lon, a.lat = PRESETS[a.preset][:2], PRESETS[a.preset][2:]
    if not (a.lon and a.lat):
        p.error("give --preset or both --lon and --lat")
    try:
        import copernicusmarine as cm
    except ImportError:
        sys.exit("pip install copernicusmarine, then run 'copernicusmarine login'")

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    box = dict(minimum_longitude=a.lon[0], maximum_longitude=a.lon[1],
               minimum_latitude=a.lat[0], maximum_latitude=a.lat[1])
    print(f"region lon {a.lon[0]}..{a.lon[1]}, lat {a.lat[0]}..{a.lat[1]}, {a.start}..{a.end}, 0..{a.depth:.0f} m")
    cm.subset(dataset_id=a.dataset or (DAILY if a.daily else MONTHLY), variables=["thetao", "uo", "vo"],
              start_datetime=a.start, end_datetime=a.end, minimum_depth=0.0, maximum_depth=a.depth,
              output_filename="forcing.nc", output_directory=str(out), overwrite=True, **box)
    cm.subset(dataset_id=STATIC, variables=["deptho"], output_filename="bathymetry.nc",
              output_directory=str(out), overwrite=True, **box)
    describe(out / "forcing.nc", out / "bathymetry.nc", a)


def describe(forcing, bathy, a):
    """Report what arrived and print a namelist fragment matching it."""
    import numpy as np
    import xarray as xr
    ds = xr.open_dataset(forcing)
    T = ds["thetao"]
    lon, lat, dep = ds["longitude"].values, ds["latitude"].values, ds["depth"].values
    print(f"\n{forcing}: {T.shape} (time, depth, lat, lon), {forcing.stat().st_size / 1e6:.1f} MB")
    print(f"  time  {str(ds.time.values[0])[:10]} .. {str(ds.time.values[-1])[-0:] and str(ds.time.values[-1])[:10]}"
          f"  ({ds.sizes['time']} steps)")
    print(f"  depth {dep[0]:.1f} .. {dep[-1]:.1f} m ({len(dep)} levels)")
    print(f"  temp  {float(T.min()):.2f} .. {float(T.max()):.2f} degC, "
          f"{100 * float(np.isnan(T.isel(time=0)).mean()):.0f}% missing (land and below the bed)")
    edges = np.r_[0, 0.5 * (dep[1:] + dep[:-1]), dep[-1] + (dep[-1] - dep[-2]) / 2]
    keep = edges <= min(a.depth, edges[-1])
    print(textwrap.dedent(f"""
        [grid]                       # a grid inside the downloaded box (coarsen as you like)
        nx = {max(8, int((lon[-1] - lon[0]) / 0.25))}
        ny = {max(8, int((lat[-1] - lat[0]) / 0.25))}
        lon = [{lon[0]:.3f}, {lon[-1]:.3f}]
        lat = [{lat[0]:.3f}, {lat[-1]:.3f}]
        depth_edges = [{", ".join(f"{x:.1f}" for x in edges[keep][:8])}]

        [ocean]
        source = "netcdf"
        file = "{forcing}"
        bathymetry = "{bathy}"
        start = "{str(ds.time.values[0])[:10]}T00:00"
        """))


if __name__ == "__main__":
    main()
