#!/usr/bin/env python3
"""Join staged GLORYS months (from get_forcing_global.py) into one dated forcing file for a hindcast.

    python make_forcing_sequence.py --in forcing_glorys_1993_2012/monthly --out forcing_glorys_1993_2012/forcing.nc

Unlike make_climatology.py this keeps the real sequence: interannual variability, El Nino, trends. Each monthly
mean is stamped at the middle of its month (GLORYS stamps the first day), so fishNET's linear interpolation
between steps has the right phase. fishNET reads the dates relative to run.start, so a run starting at the
file's first month follows the calendar; a run longer than the file loops back to its start (with a warning).
"""
import argparse, sys
from pathlib import Path
import numpy as np
import xarray as xr

VARS = ["thetao", "uo", "vo"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="src", required=True, help="folder of staged YYYYMM.nc files")
    ap.add_argument("--out", required=True)
    ap.add_argument("--start", help="first month to include, YYYY-MM")
    ap.add_argument("--end", help="last month to include, YYYY-MM")
    a = ap.parse_args()
    files = sorted(Path(a.src).glob("[12][0-9][0-9][0-9][01][0-9].nc"))
    key = lambda p: p.stem[:4] + "-" + p.stem[4:]
    files = [p for p in files if (not a.start or key(p) >= a.start) and (not a.end or key(p) <= a.end)]
    if not files:
        sys.exit("no staged months found")
    months = [np.datetime64(key(p)) for p in files]
    gaps = [str(m) for m, n in zip(months, months[1:]) if (n - m).astype(int) != 1]
    if gaps:
        sys.exit(f"missing months after {', '.join(gaps)}: download them first")
    ds = xr.open_mfdataset(files, combine="by_coords")[VARS]
    start = ds.time.values.astype("datetime64[M]")
    mid = start.astype("datetime64[s]") + ((start + 1).astype("datetime64[s]") - start.astype("datetime64[s]")) // 2
    ds = ds.assign_coords(time=mid)
    ds.attrs = dict(note=f"GLORYS12 monthly means {key(files[0])}..{key(files[-1])}, coarsened to 1 deg, stamped at "
                         "mid-month; a dated sequence (not a climatology). Built by make_forcing_sequence.py.")
    enc = {v: dict(zlib=True, complevel=4, dtype="float32") for v in VARS}
    out = Path(a.out)
    ds.to_netcdf(out, encoding=enc)
    print(f"wrote {out} ({out.stat().st_size / 1e6:.0f} MB): {ds.sizes['time']} months "
          f"{str(ds.time.values[0])[:10]} .. {str(ds.time.values[-1])[:10]}, "
          f"{ds.sizes['latitude']}x{ds.sizes['longitude']}, {ds.sizes['depth']} levels")


if __name__ == "__main__":
    main()
