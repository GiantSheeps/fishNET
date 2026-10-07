#!/usr/bin/env python3
"""Build a global monthly forcing file by streaming GLORYS one month at a time.

A global 1/12 deg subset is ~2.3 GB per month and ~129 GB in memory if processed whole, so each
month is downloaded, coarsened, written to a small staging file and the raw copy deleted. Resumable:
months already staged are skipped, so it can be interrupted and restarted.

    python get_forcing_global.py                       # 1999-12 .. 2002-01, 0-4000 m, coarsen to 1 deg
    python get_forcing_global.py --start 1993-01 --end 2020-12

Then feed the result to make_climatology.py to fit the smooth annual cycle.
"""
import argparse, datetime as dtm, shutil
from pathlib import Path

import numpy as np
import xarray as xr

DATASET = "cmems_mod_glo_phy_my_0.083deg_P1M-m"
VARS = ["thetao", "uo", "vo"]


def months(start, end):
    y0, m0 = (int(x) for x in start.split("-"))
    y1, m1 = (int(x) for x in end.split("-"))
    y, m = y0, m0
    while (y, m) <= (y1, m1):
        yield y, m
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--start", default="1999-12")
    p.add_argument("--end", default="2002-01")
    p.add_argument("--depth", type=float, default=4000.0, help="deepest level (m)")
    p.add_argument("--coarsen", type=int, default=12, help="average NxN native cells (12 = 1/12 deg -> 1 deg)")
    p.add_argument("--out", default="forcing_global")
    p.add_argument("--keep-raw", action="store_true")
    a = p.parse_args()

    import copernicusmarine as cm
    out = Path(a.out); stage = out / "monthly"; raw = out / "raw"
    stage.mkdir(parents=True, exist_ok=True); raw.mkdir(parents=True, exist_ok=True)
    todo = list(months(a.start, a.end))
    print(f"{len(todo)} months, global, 0-{a.depth:.0f} m, coarsen {a.coarsen}x "
          f"(1/12 deg -> {a.coarsen / 12:.2f} deg)... staging in {stage}")

    for k, (y, m) in enumerate(todo, 1):
        small = stage / f"{y}{m:02d}.nc"
        if small.exists():
            print(f"[{k:3d}/{len(todo)}] {y}-{m:02d} staged already")
            continue
        last = (dtm.date(y + (m == 12), m % 12 + 1, 1) - dtm.timedelta(days=1)).day
        big = raw / f"{y}{m:02d}.nc"
        if not big.exists():
            cm.subset(dataset_id=DATASET, variables=VARS,
                      minimum_longitude=-180.0, maximum_longitude=179.92,
                      minimum_latitude=-80.0, maximum_latitude=90.0,
                      minimum_depth=0.0, maximum_depth=a.depth,
                      start_datetime=f"{y}-{m:02d}-01", end_datetime=f"{y}-{m:02d}-{last}",
                      output_filename=big.name, output_directory=str(raw),
                      overwrite=True, disable_progress_bar=True)
        ds = xr.open_dataset(big)
        c = ds.coarsen(longitude=a.coarsen, latitude=a.coarsen, boundary="trim").mean(skipna=True)
        enc = {v: dict(zlib=True, complevel=4, dtype="float32") for v in VARS}
        c.to_netcdf(small, encoding=enc)
        c.close(); ds.close()
        if not a.keep_raw:
            big.unlink()
        print(f"[{k:3d}/{len(todo)}] {y}-{m:02d}  {small.stat().st_size / 1e6:6.1f} MB staged"
              f"  ({c.sizes['latitude']}x{c.sizes['longitude']}, {c.sizes['depth']} levels)")

    files = sorted(stage.glob("*.nc"))
    print(f"\nconcatenating {len(files)} months")
    ds = xr.open_mfdataset(files, combine="by_coords")
    enc = {v: dict(zlib=True, complevel=4, dtype="float32") for v in VARS}
    dest = out / "monthly_all.nc"
    ds.to_netcdf(dest, encoding=enc)
    print(f"wrote {dest} ({dest.stat().st_size / 1e6:.0f} MB, {ds.sizes['time']} steps, "
          f"{ds.sizes['latitude']}x{ds.sizes['longitude']}, {ds.sizes['depth']} levels)")
    if not a.keep_raw:
        shutil.rmtree(raw, ignore_errors=True)
    print(f"\nnext:  python make_climatology.py --in {dest} --out {out}/forcing.nc --steps 73 --harmonics 2")


if __name__ == "__main__":
    main()
