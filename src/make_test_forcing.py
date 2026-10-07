#!/usr/bin/env python3
"""Write a GLORYS-lookalike forcing file for testing the netCDF path without a Copernicus account.

The data are invented, but the packaging copies CMEMS GLORYS as closely as it can: variable names
(thetao, uo, vo, deptho), descending latitude, packed short integers with scale_factor/add_offset,
_FillValue land masking that deepens with depth, time as "hours since 1950-01-01" on a monthly axis,
depth positive downward on the GLORYS level set, and bathymetry in a separate static file.

    python make_test_forcing.py --out test_forcing            # then point [ocean] at it
    python make_test_forcing.py --out x --lon 0 360           # 0..360 longitude convention
"""
import argparse
import numpy as np
import netCDF4 as nc
from pathlib import Path

GLORYS_LEVELS = np.array([0.494, 1.541, 2.646, 3.819, 5.078, 6.441, 7.930, 9.573, 11.405, 13.467,
                          15.810, 18.496, 21.599, 25.211, 29.445, 34.434, 40.344, 47.374, 55.764,
                          65.807, 77.854, 92.326, 109.729, 130.666, 155.851, 186.126, 222.475,
                          266.040, 318.127, 380.213, 453.938, 541.089, 643.567, 763.333, 902.339])


def pack(var, data, fill=-32767):
    """Store as short integers with scale_factor/add_offset, exactly as CMEMS does."""
    good = np.isfinite(data)
    lo, hi = float(data[good].min()), float(data[good].max())
    scale = max((hi - lo) / 65000.0, 1e-9)
    var.setncatts(dict(scale_factor=scale, add_offset=(hi + lo) / 2, _FillValue=np.int16(fill),
                       missing_value=np.int16(fill)))
    var[:] = np.ma.masked_array(np.where(good, data, 0.0), mask=~good)      # land/below-bed -> _FillValue


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default="test_forcing")
    p.add_argument("--lon", nargs=2, type=float, default=[-70.0, -10.0])
    p.add_argument("--lat", nargs=2, type=float, default=[20.0, 60.0])
    p.add_argument("--res", type=float, default=0.5, help="degrees (GLORYS itself is 1/12)")
    p.add_argument("--months", type=int, default=12)
    p.add_argument("--levels", type=int, default=20)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    lon = np.arange(a.lon[0], a.lon[1], a.res)
    lat = np.arange(a.lat[0], a.lat[1], a.res)[::-1]            # CMEMS ships north to south
    dep = GLORYS_LEVELS[np.unique(np.linspace(0, len(GLORYS_LEVELS) - 1, a.levels).astype(int))]
    t0 = (np.datetime64("2010-01-15") - np.datetime64("1950-01-01")) / np.timedelta64(1, "h")
    time = t0 + np.arange(a.months) * 730.5                     # monthly, hours since 1950
    LON, LAT = np.meshgrid(lon, lat)

    # A continent in the east with a shelf, plus a deep basin: bathymetry first, fields masked by it.
    coast = a.lon[1] - 6 - 2 * np.sin(np.radians(3 * LAT))
    off = np.clip(coast - LON, 0, None)
    depth = np.clip(60 + 3000 * (1 - np.exp(-off / 4.0)), 0, 4000)
    depth[LON >= coast] = np.nan                                # land
    seamount = 3000 * np.exp(-(((LON + 40) ** 2 + (LAT - 38) ** 2) / 8))
    depth = np.where(np.isnan(depth), np.nan, np.maximum(depth - seamount, 50))

    yr = 2 * np.pi * (np.arange(a.months) % 12) / 12
    sst = 26 - 0.42 * (LAT - 20) + 3.5 * np.sin(yr - 1.9)[:, None, None]
    thermo = 12 * np.exp(-dep[:, None, None] / 220)
    temp = 3.5 + (sst[:, None] - 3.5) * np.exp(-dep[None, :, None, None] / 180) + thermo
    temp = temp + 0.25 * rng.standard_normal(temp.shape)

    # A wind-driven gyre with a western boundary current, decaying with depth.
    xs = (LON - lon[0]) / (lon[-1] - lon[0])
    ys = (LAT - lat[-1]) / (lat[0] - lat[-1])
    psi = 4e4 * np.sin(2 * np.pi * ys) * (1 - xs) * (1 - np.exp(-xs / 0.05))
    decay = np.exp(-dep / 600)[:, None, None]
    dy_m = a.res * 111e3
    dx_m = dy_m * np.cos(np.radians(LAT))
    u = -(np.gradient(psi, axis=0) / dy_m)[None] * decay                    # m/s from the streamfunction
    v = (np.gradient(psi, axis=1) / dx_m)[None] * decay
    u, v = [np.repeat(x[None], a.months, 0) * (1 + 0.15 * np.sin(yr)[:, None, None, None]) for x in (u, v)]

    wet3 = np.isfinite(depth)[None] & (dep[:, None, None] < np.nan_to_num(depth, nan=-1)[None])
    for x in (temp, u, v):
        x[:, ~wet3] = np.nan

    f = nc.Dataset(out / "forcing.nc", "w", format="NETCDF4")
    f.setncatts(dict(title="synthetic GLORYS-lookalike forcing (invented data, CMEMS packaging)",
                     institution="fishNET make_test_forcing.py", source="not real reanalysis output"))
    for name, size in (("time", None), ("depth", len(dep)), ("latitude", len(lat)), ("longitude", len(lon))):
        f.createDimension(name, size)
    tv = f.createVariable("time", "f8", ("time",))
    tv.setncatts(dict(units="hours since 1950-01-01 00:00:00", calendar="gregorian", standard_name="time"))
    tv[:] = time
    dv = f.createVariable("depth", "f4", ("depth",))
    dv.setncatts(dict(units="m", positive="down", standard_name="depth"))
    dv[:] = dep
    for name, vals, unit, std in (("latitude", lat, "degrees_north", "latitude"),
                                  ("longitude", lon, "degrees_east", "longitude")):
        cv = f.createVariable(name, "f4", (name,))
        cv.setncatts(dict(units=unit, standard_name=std))
        cv[:] = vals
    for name, data, unit, std in (("thetao", temp, "degrees_C", "sea_water_potential_temperature"),
                                  ("uo", u, "m s-1", "eastward_sea_water_velocity"),
                                  ("vo", v, "m s-1", "northward_sea_water_velocity")):
        var = f.createVariable(name, "i2", ("time", "depth", "latitude", "longitude"), zlib=True)
        pack(var, data)
        var.setncatts(dict(units=unit, standard_name=std))
    f.close()

    b = nc.Dataset(out / "bathymetry.nc", "w", format="NETCDF4")
    b.setncatts(dict(title="synthetic GLORYS-lookalike static file"))
    for name, size in (("latitude", len(lat)), ("longitude", len(lon))):
        b.createDimension(name, size)
        cv = b.createVariable(name, "f4", (name,))
        cv.units = "degrees_north" if name == "latitude" else "degrees_east"
        cv[:] = lat if name == "latitude" else lon
    var = b.createVariable("deptho", "i2", ("latitude", "longitude"), zlib=True)
    pack(var, depth)
    var.setncatts(dict(units="m", standard_name="sea_floor_depth_below_geoid", positive="down"))
    b.close()

    print(f"wrote {out / 'forcing.nc'} ({(out / 'forcing.nc').stat().st_size / 1e6:.2f} MB) and {out / 'bathymetry.nc'}")
    print(f"  {a.months} monthly steps, {len(dep)} levels to {dep[-1]:.0f} m, {len(lat)}x{len(lon)} cells, "
          f"lat descending, packed int16, {100 * np.mean(~wet3):.0f}% masked")
    print(f"\n[ocean]\nsource = \"netcdf\"\nfile = \"{out / 'forcing.nc'}\"\n"
          f"bathymetry = \"{out / 'bathymetry.nc'}\"\nstart = \"2010-01-15T00:00\"")


if __name__ == "__main__":
    main()
