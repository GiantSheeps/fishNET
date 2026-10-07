#!/usr/bin/env python3
"""fishNET dashboard: a browser GUI to launch runs, watch them live and animate the output.

    python3 dashboard.py                 # serve on http://127.0.0.1:8765 and open a browser
    python3 dashboard.py --port 9000 --no-browser

Standard library only (plus numpy/netCDF4, which fishNET already needs). The server
launches `fishnet.py` as a subprocess, streams its log, and reads the netCDF output while
it is still being written, so maps, animations and time series fill in as the run goes.
"""
import os

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")   # read output files while the model writes them

import sys, json, re, ast, time, threading, subprocess, webbrowser, argparse, shutil, mimetypes, tomllib
from collections import deque
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from pathlib import Path
from urllib.parse import urlparse, parse_qs, unquote

import numpy as np
import netCDF4 as nc

ROOT = Path(__file__).resolve().parent
WEB = ROOT / "webui"
RUNS_DEFAULT = ROOT / "runs"
CONFIG_DIR = RUNS_DEFAULT / "_configs"

LEVEL_VARS = ("temp", "N", "P", "Z", "krill", "D", "O2", "par", "u", "v", "w")       # (time, depth, lat, lon)
SPECIES_STAGE_VARS = ("fish_biomass", "fish_numbers")                   # (time, species, stage, lat, lon)
SPECIES_VARS = ("agent_biomass",)                                       # (time, species, lat, lon)


# ------------------------------------------------------------------ namelist handling
def split_sections(text):
    """[(header or None, [lines]), ...] preserving the file verbatim."""
    out, cur, head = [], [], None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("[") and not s.startswith("#"):
            out.append((head, cur))
            head, cur = s.split("#")[0].strip(), [line]
        else:
            cur.append(line)
    out.append((head, cur))
    return out


def toml_value(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    return json.dumps(str(v))


def set_key(text, section, key, value, occurrence=0):
    """Set `key = value` inside `section` ([run], [[species]], ...), preserving comments."""
    secs, seen, done = split_sections(text), -1, False
    pat = re.compile(rf"^(\s*){re.escape(key)}(\s*)=\s*([^#]*)(#.*)?$")
    for si, (head, lines) in enumerate(secs):
        if head != section:
            continue
        seen += 1
        if seen != occurrence:
            continue
        for li, line in enumerate(lines):
            m = pat.match(line)
            if m:
                lines[li] = f"{m.group(1)}{key}{m.group(2) or ' '}= {toml_value(value)}" + (
                    f"   {m.group(4)}" if m.group(4) else "")
                done = True
                break
        if not done:                                        # key absent: add it under the header
            lines.insert(1, f"{key} = {toml_value(value)}")
            done = True
        break
    if not done:                                            # section absent: append one
        secs.append((section, [section, f"{key} = {toml_value(value)}"]))
    return "\n".join("\n".join(l for l in lines) for _, lines in secs)


def absolutise(text, base):
    """Rewrite species `file = "..."` paths and `out_dir` so the config can live anywhere."""
    out = []
    for line in text.splitlines():
        m = re.match(r"^(\s*file\s*=\s*)\"([^\"]+)\"(.*)$", line)
        if m and not Path(m.group(2)).is_absolute():
            p = base / m.group(2)
            p = p if p.exists() else base / Path(m.group(2)).name
            line = f'{m.group(1)}"{p.resolve()}"{m.group(3)}'
        out.append(line)
    text = "\n".join(out)
    text = set_key(text, "[run]", "out_dir", str(RUNS_DEFAULT))
    text = set_key(text, "[ocean]", "sw_cache_dir", str(ROOT / "sw_cache"))
    return text


def apply_overrides(text, ov):
    """ov: {"run.days": 30, "species.0.fishing_F": 0.2, ...} from the UI's quick controls."""
    for path, value in (ov or {}).items():
        keys = path.split(".")
        if keys[0] == "species" and len(keys) == 3:
            text = set_key(text, "[[species]]", keys[2], value, occurrence=int(keys[1]))
        elif len(keys) == 2:
            text = set_key(text, f"[{keys[0]}]", keys[1], value)
    return text


# ------------------------------------------------------------------ run process
STEP_RE = re.compile(
    r"\[\s*(\d+)/(\d+)\]\s+(\d{4}-\d\d-\d\d \d\d:\d\d)\s+day\s+([-\d.]+)\s*\|\s*agents\s+([\d,]+)"
    r"\s*\(\+(\d+)\s*-(\d+)\)\s*\|\s*IBM cells\s+(\d+).*?\|\s*([\d.]+)\s*s/step\s+ETA\s+([-\d.]+)")
SPECIES_RE = re.compile(r"^\s{6}(\S+)\s+B\s+([\d.e+-]+)\s+t\s+\[(.*?)\]\s+agents\s+(\d+)")


class Runner:
    """One model subprocess at a time, with its log and parsed progress."""

    def __init__(self):
        self.lock = threading.Lock()
        self.proc = None
        self.lines = deque(maxlen=6000)
        self.first = 0                      # index of lines[0] in the whole stream
        self.count = 0
        self.progress = {}
        self.species = {}
        self.run_name = None
        self.cfg_path = None
        self.started = None
        self.finished = None
        self.returncode = None
        self.command = None

    def running(self):
        return self.proc is not None and self.proc.poll() is None

    def start(self, cfg_text, run_name, ensemble=0, resume=False):
        if self.running():
            raise RuntimeError("a run is already in progress")
        CONFIG_DIR.mkdir(parents=True, exist_ok=True)
        cfg_path = CONFIG_DIR / f"{run_name}.toml"
        cfg_path.write_text(cfg_text)
        cmd = [sys.executable, "-u", str(ROOT / "fishnet.py"), str(cfg_path)]
        if ensemble:
            cmd += ["--ensemble", str(ensemble)]
        if resume:
            cmd.append("--resume")
        with self.lock:
            self.lines.clear()
            self.first = self.count = 0
            self.progress, self.species = {}, {}
            self.run_name, self.cfg_path = run_name, cfg_path
            self.started, self.finished, self.returncode = time.time(), None, None
            self.command = " ".join(cmd)
            self.proc = subprocess.Popen(cmd, cwd=str(ROOT), stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT, text=True, bufsize=1)
        threading.Thread(target=self._pump, args=(self.proc,), daemon=True).start()
        return {"run": run_name, "config": str(cfg_path), "command": self.command}

    def stop(self):
        if self.running():
            self.proc.terminate()
            try:
                self.proc.wait(timeout=8)
            except subprocess.TimeoutExpired:
                self.proc.kill()
            return True
        return False

    def _pump(self, proc):
        for line in proc.stdout:
            line = line.rstrip("\n")
            with self.lock:
                self.lines.append(line)
                self.count += 1
                self.first = self.count - len(self.lines)
                m = STEP_RE.search(line)
                if m:
                    self.progress = dict(step=int(m[1]), nsteps=int(m[2]), date=m[3], day=float(m[4]),
                                         agents=int(m[5].replace(",", "")), born=int(m[6]), died=int(m[7]),
                                         ibm_cells=int(m[8]), s_per_step=float(m[9]), eta_min=float(m[10]))
                s = SPECIES_RE.match(line)
                if s:
                    self.species[s[1]] = dict(biomass_t=float(s[2]),
                                              stages=[float(x.strip("% ")) / 100 for x in s[3].split()],
                                              agents=int(s[4]))
        proc.wait()
        with self.lock:
            self.finished, self.returncode = time.time(), proc.returncode

    def status(self, since=0):
        with self.lock:
            start = max(since - self.first, 0)
            lines = list(self.lines)[start:] if since < self.count else []
            return dict(running=self.running(), run=self.run_name, progress=dict(self.progress),
                        species={k: dict(v) for k, v in self.species.items()},
                        lines=lines, next=self.count, first=self.first, command=self.command,
                        started=self.started, finished=self.finished, returncode=self.returncode,
                        elapsed=(time.time() - self.started) if self.started else None)


RUNNER = Runner()


# ------------------------------------------------------------------ netCDF access
class Store:
    """Opens run output, re-opening while a run is live so new records appear.

    HDF5 is not thread-safe, so `lock` guards every read as well as every open: the
    threaded server would otherwise corrupt the library's heap under concurrent requests.
    """

    def __init__(self):
        self.lock = threading.RLock()
        self.ds = {}            # path -> (Dataset, opened_at, size)
        self.agents = {}        # path -> {"n": records scanned, "times": [...], "start": [...]}
        self.ranges = {}

    def _live(self, run):
        return RUNNER.running() and RUNNER.run_name == run

    def open(self, run, fname, force=False):
        path = run_dir(run) / fname
        if not path.exists():
            return None
        with self.lock:
            key = str(path)
            ds, opened, size = self.ds.get(key, (None, 0, -1))
            now, cur = time.time(), path.stat().st_size
            stale = force or ds is None or (self._live(run) and now - opened > 0.4) or cur != size
            if stale:
                if ds is not None:
                    try:
                        ds.close()
                    except Exception:
                        pass
                for attempt in range(4):
                    try:
                        ds = nc.Dataset(str(path), "r")
                        break
                    except OSError:
                        time.sleep(0.15 * (attempt + 1))
                else:
                    return None
                ds.set_auto_mask(False)
                self.ds[key] = (ds, now, cur)
            return ds

    def read(self, run, fname, fn, default=None):
        """Run `fn(dataset)`, re-opening once if the writer was mid-flush (HDF errors on a live file)."""
        for attempt in range(3):
            ds = self.open(run, fname, force=attempt > 0)
            if ds is None:
                return default
            try:
                return fn(ds)
            except (RuntimeError, OSError, IndexError, KeyError) as e:
                last = e
                time.sleep(0.15 * (attempt + 1))
        print(f"read {run}/{fname}: {type(last).__name__}: {last}", file=sys.stderr, flush=True)
        return default

    def close_run(self, run):
        with self.lock:
            for key in [k for k in self.ds if k.startswith(str(run_dir(run)))]:
                try:
                    self.ds.pop(key)[0].close()
                except Exception:
                    pass

    def agent_index(self, run):
        """Snapshot times in agents.nc and where each starts, extended incrementally."""
        return self.read(run, "agents.nc", self._agent_index, {"times": [], "start": [], "n": 0})

    def _agent_index(self, ds):
        if "rec" not in ds.dimensions:
            return {"times": [], "start": [], "n": 0}
        key, n = ds.filepath(), len(ds.dimensions["rec"])
        with self.lock:
            idx = self.agents.setdefault(key, {"n": 0, "times": [], "start": []})
            if n < idx["n"]:                                   # file was rewritten
                idx = self.agents[key] = {"n": 0, "times": [], "start": []}
            if n > idx["n"]:
                t = np.asarray(ds["time"][idx["n"]:n], dtype="f8")
                if t.size:
                    cuts = np.flatnonzero(np.diff(t)) + 1
                    starts = np.r_[0, cuts] + idx["n"]
                    times = t[np.r_[0, cuts]]
                    if idx["times"] and abs(times[0] - idx["times"][-1]) < 1e-9:
                        starts, times = starts[1:], times[1:]   # same snapshot continued
                    idx["times"] += [float(x) for x in times]
                    idx["start"] += [int(x) for x in starts]
                idx["n"] = n
            return {"times": list(idx["times"]), "start": list(idx["start"]), "n": idx["n"]}


STORE = Store()


def run_dir(run):
    return (RUNS_DEFAULT / run).resolve()


def latest_restart_day(run):
    """Day of the newest restart file for a run, or None."""
    files = sorted((run_dir(run) / "restart").glob("restart_day*.pkl"))
    return float(files[-1].stem.removeprefix("restart_day")) if files else None


def valid_run(run):
    p = run_dir(run)
    return p.is_dir() and str(p).startswith(str(RUNS_DEFAULT.resolve())) and not run.startswith("_")


def all_species_groups():
    groups = {}
    species_dir = ROOT / "species"
    if species_dir.is_dir():
        for p in sorted(species_dir.glob("*.toml")):
            try:
                with open(p, "rb") as f:
                    data = tomllib.load(f)
                if "group" in data:
                    groups[p.stem] = str(data["group"])
            except Exception:
                pass
    return groups


def species_groups(ds, names):
    """Functional group per species, read from the `group` key in each species file. fishNET stores the
    whole resolved config in the output's `config` attribute, so no extra model output is needed.
    Species without a group fall back to reading species/<name>.toml on disk, or "other"."""
    groups = {}
    if "config" in ds.ncattrs():
        try:
            cfg = ast.literal_eval(ds.getncattr("config"))
            for spec in cfg.get("species", []):
                if spec.get("name"):
                    groups[str(spec["name"])] = str(spec.get("group", "other"))
        except Exception:
            pass
    disk_groups = all_species_groups()
    for n in names:
        if n not in groups or groups[n] in ("other", "") or n in ("squid", "shark", "whale", "dolphin"):
            if n in disk_groups:
                groups[n] = disk_groups[n]
    return [groups.get(n, "other") for n in names]


def meta_attrs(ds):
    keys = ("species", "colors", "traits", "stages", "causes", "behaviors")
    out = {k: [s for s in ds.getncattr(k).split(",") if s] for k in keys if k in ds.ncattrs()}
    out["groups"] = species_groups(ds, out.get("species", []))
    out["start"] = ds.getncattr("start") if "start" in ds.ncattrs() else "2000-01-01 00:00:00"
    for k in ("dt_hours", "hybrid", "engine", "behavior_mode", "ocean_area"):
        if k in ds.ncattrs():
            v = ds.getncattr(k)
            out[k] = float(v) if isinstance(v, (np.floating, np.integer, float, int)) else str(v)
    return out


def list_runs():
    out = []
    if not RUNS_DEFAULT.exists():
        return out
    for p in sorted(RUNS_DEFAULT.iterdir()):
        if not p.is_dir() or p.name.startswith("_"):
            continue
        series, fields = p / "series.nc", p / "fields.nc"
        if not series.exists() and not fields.exists():
            continue
        info = dict(name=p.name, mtime=max((f.stat().st_mtime for f in p.glob("*.nc")), default=0),
                    size=sum(f.stat().st_size for f in p.glob("*.nc")),
                    running=RUNNER.running() and RUNNER.run_name == p.name,
                    ensemble=(p / "ensemble.nc").exists(),
                    figures=sorted(f.name for f in (p / "figures").glob("*")) if (p / "figures").is_dir() else [])
        info["restart_day"] = latest_restart_day(p.name)
        try:
            info["cfg_days"] = tomllib.loads((CONFIG_DIR / f"{p.name}.toml").read_text())["run"].get("days")
        except (OSError, tomllib.TOMLDecodeError, KeyError):
            pass
        ds = STORE.open(p.name, "series.nc")
        if ds is not None:
            info["steps"] = len(ds.dimensions["time"])
            info["day"] = float(ds["time"][-1]) if info["steps"] else 0.0
            info.update(meta_attrs(ds))
        out.append(info)
    out.sort(key=lambda r: -r["mtime"])
    return out


def parse_species(v):
    """The `species` query value: an index, -1 for every species, or "0,3,7" for a subset
    (what the dashboard's species filter sends when it is summing only the species on show)."""
    txt = str(v if v is not None else -1)
    if "," in txt:
        idx = [int(x) for x in txt.split(",") if x.strip() != ""]
        return idx or -1
    return int(txt)


def pick_species(a, species):
    """Select or sum axis 0 of a (species, ...) array: an index, a list of indices, or -1 for all."""
    if isinstance(species, list):
        return np.nansum(a[species], 0)
    return a[species] if species >= 0 else np.nansum(a, 0)


def field_slice(ds, var, t, species=-1, stage=-1, depth=0, trait=0):
    """A (lat, lon) float32 frame, NaN over land, with -1 meaning 'sum over that dimension'.
    `species` may also be a list of indices, summed together."""
    mask = np.asarray(ds["mask"][:]) > 0.5
    if var == "speed":
        u, v = (np.asarray(ds[k][t, depth], dtype="f4") for k in ("u", "v"))
        a = np.hypot(u, v)
    elif var == "bottom_depth":
        a = np.asarray(ds["bottom_depth"][:], dtype="f4")
    elif var in LEVEL_VARS:
        a = np.asarray(ds[var][t, depth], dtype="f4")
    elif var in SPECIES_STAGE_VARS:
        a = np.asarray(ds[var][t], dtype="f4")                  # (species, stage, lat, lon)
        a = pick_species(a, species)
        a = a[stage] if stage >= 0 else np.nansum(a, 0)
    elif var in SPECIES_VARS:
        a = pick_species(np.asarray(ds[var][t], dtype="f4"), species)
    elif var == "trait_mean":
        one = species[0] if isinstance(species, list) and species else (species if not isinstance(species, list) else 0)
        a = np.asarray(ds[var][t, max(one, 0), trait], dtype="f4")   # a breeding value cannot be summed
    elif var == "ibm":
        a = np.asarray(ds[var][t], dtype="f4")
    else:
        raise KeyError(var)
    a = np.where(mask, a, np.nan).astype("f4")
    return a


def percentiles(a, log=False):
    v = a[np.isfinite(a)]
    if log:
        v = v[v > 0]
    if v.size == 0:
        return 0.0, 1.0
    lo, hi = (float(x) for x in np.percentile(v, [1.0, 99.0]))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.min(v)), float(np.max(v)) if v.size else (0.0, 1.0)
    return lo, (hi if hi > lo else lo + 1e-12)


# ------------------------------------------------------------------ HTTP
class Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    server_version = "fishNET-dashboard"

    def log_message(self, fmt, *args):
        pass

    # ---- helpers
    def send_json(self, obj, code=200):
        body = json.dumps(obj, allow_nan=False, default=jsonable).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def send_bin(self, data, meta):
        self.send_response(200)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("X-Meta", json.dumps(meta, allow_nan=False, default=jsonable))
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def send_error_json(self, code, msg):
        self.send_json({"error": msg}, code)

    def body(self):
        n = int(self.headers.get("Content-Length", 0))
        return json.loads(self.rfile.read(n) or b"{}")

    # ---- routing
    def do_GET(self):
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        try:
            if u.path.startswith("/api/"):
                with STORE.lock:                      # HDF5 reads are serialised across request threads
                    return self.api_get(u.path[5:], q)
            return self.static(u.path)
        except BrokenPipeError:
            pass
        except Exception as e:
            print(f"GET {self.path} -> {type(e).__name__}: {e}", file=sys.stderr, flush=True)
            self.send_error_json(500, f"{type(e).__name__}: {e}")

    def do_POST(self):
        u = urlparse(self.path)
        try:
            with STORE.lock:
                return self.api_post(u.path[5:], self.body())
        except BrokenPipeError:
            pass
        except Exception as e:
            print(f"POST {self.path} -> {type(e).__name__}: {e}", file=sys.stderr, flush=True)
            self.send_error_json(400, f"{type(e).__name__}: {e}")

    # ---- static files
    def static(self, path):
        rel = "index.html" if path in ("/", "") else unquote(path.lstrip("/"))
        f = (WEB / rel).resolve()
        if not str(f).startswith(str(WEB.resolve())) or not f.is_file():
            return self.send_error_json(404, "not found")
        data = f.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", mimetypes.guess_type(f.name)[0] or "application/octet-stream")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    # ---- API
    def api_get(self, route, q):
        run = q.get("run", "")
        if route == "state":
            return self.send_json(dict(runs=list_runs(), status=RUNNER.status(int(q.get("since", 0))),
                                       namelists=[p.name for p in sorted(ROOT.glob("*.toml"))],
                                       species_groups=all_species_groups(),
                                       root=str(ROOT)))
        if route == "status":
            return self.send_json(RUNNER.status(int(q.get("since", 0))))
        if route == "namelist":
            name = q.get("name", "namelist.toml")
            p = (ROOT / name).resolve()
            if p.parent not in (ROOT, CONFIG_DIR.resolve()) or not p.is_file():
                return self.send_error_json(404, "no such namelist")
            return self.send_json({"name": name, "text": p.read_text()})

        if route in ("meta", "times", "series", "field", "agents", "range", "vectors", "figure", "validation", "valmap"):
            if not valid_run(run):
                return self.send_error_json(404, f"no run '{run}'")
        if route == "meta":
            return self.send_json(self.meta(run))
        if route == "times":
            return self.send_json(self.times(run))
        if route == "series":
            return self.send_json(self.series(run, q))
        if route == "field":
            return self.field(run, q)
        if route == "vectors":
            return self.vectors(run, q)
        if route == "agents":
            return self.agents(run, q)
        if route == "range":
            return self.send_json(self.range(run, q))
        if route == "validation":
            return self.validation(run, q)
        if route == "valmap":
            return self.valmap(run, q)
        if route == "figure":
            p = (run_dir(run) / "figures" / Path(q.get("name", "")).name)
            if not p.is_file():
                return self.send_error_json(404, "no such figure")
            return self.sendfile(p)
        return self.send_error_json(404, "unknown endpoint")

    def api_post(self, route, body):
        if route == "start":
            name = re.sub(r"[^A-Za-z0-9_.-]", "_", body.get("name") or "dashboard_run")
            src = body.get("text")
            if not src:
                p = (ROOT / (body.get("namelist") or "namelist.toml")).resolve()
                src = p.read_text()
            text = apply_overrides(src, body.get("overrides"))
            text = set_key(text, "[run]", "name", name)
            text = absolutise(text, ROOT)
            try:                                              # catch a bad namelist here, not in a traceback
                tomllib.loads(text)
            except tomllib.TOMLDecodeError as e:
                line = int(m[1]) if (m := re.search(r"at line (\d+)", str(e))) else 0
                near = text.splitlines()[line - 1].strip() if 0 < line <= len(text.splitlines()) else ""
                return self.send_error_json(400, f"{name}.toml is not valid TOML: {e}" + (f"\n  {near}" if near else ""))
            if body.get("fresh", True):
                shutil.rmtree(run_dir(name), ignore_errors=True)
            STORE.close_run(name)
            return self.send_json(RUNNER.start(text, name, int(body.get("ensemble", 0) or 0)))
        if route == "resume":                                 # continue a run from its newest restart file
            name = body.get("name", "")
            cfg_path = CONFIG_DIR / f"{name}.toml"
            if not valid_run(name) or not cfg_path.is_file():
                return self.send_error_json(404, f"no saved config for run '{name}' in {CONFIG_DIR}")
            if latest_restart_day(name) is None:
                return self.send_error_json(400, f"'{name}' has no restart files (runs/{name}/restart/); it was run "
                                                 f"without restart_every_days, so it can only be started again")
            text = cfg_path.read_text()
            if body.get("days"):
                text = set_key(text, "[run]", "days", int(body["days"]))
            STORE.close_run(name)
            return self.send_json(RUNNER.start(text, name, resume=True))
        if route == "stop":
            return self.send_json({"stopped": RUNNER.stop()})
        if route == "figures":
            run = body.get("run", "")
            if not valid_run(run):
                return self.send_error_json(404, "no such run")
            cmd = [sys.executable, "-u", str(ROOT / "plot.py"), str(run_dir(run))]
            if body.get("no_anim", True):
                cmd.append("--no-anim")
            subprocess.Popen(cmd, cwd=str(ROOT), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return self.send_json({"started": True})
        if route == "preview":
            src = (ROOT / (body.get("namelist") or "namelist.toml")).read_text()
            return self.send_json({"text": apply_overrides(src, body.get("overrides"))})
        return self.send_error_json(404, "unknown endpoint")

    def sendfile(self, p):
        data = p.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", mimetypes.guess_type(p.name)[0] or "application/octet-stream")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    # ---- data endpoints
    def meta(self, run):
        def grid(F):
            mask = np.asarray(F["mask"][:]) > 0.5
            bd = np.asarray(F["bottom_depth"][:], dtype="f8")
            return dict(lon=[float(x) for x in F["lon"][:]], lat=[float(x) for x in F["lat"][:]],
                        depth=[float(x) for x in F["depth"][:]],
                        nx=len(F.dimensions["lon"]), ny=len(F.dimensions["lat"]), nz=len(F.dimensions["depth"]),
                        mask=[int(x) for x in mask.ravel()],
                        fields=[v for v in F.variables if F[v].dimensions[:1] == ("time",)],
                        max_depth=float(np.nanmax(np.where(mask, bd, np.nan))) if mask.any() else 0.0,
                        **meta_attrs(F))
        m = STORE.read(run, "fields.nc", grid) or STORE.read(run, "series.nc", meta_attrs)
        if m is None:
            return {"error": "no output yet", "run": run}
        m = dict(run=run, **m)
        m.update(self.times(run))
        cfg = CONFIG_DIR / f"{run}.toml"
        m["config_path"] = str(cfg) if cfg.exists() else None
        return m

    def times(self, run):
        ft = STORE.read(run, "fields.nc",
                        lambda F: [float(x) for x in F["time"][:]] if len(F.dimensions["time"]) else [], []) or []
        sn = STORE.read(run, "series.nc",
                        lambda S: (len(S.dimensions["time"]), float(S["time"][-1]) if len(S.dimensions["time"]) else 0.0),
                        (0, 0.0))
        a = STORE.agent_index(run)
        return dict(field_times=ft, n_series=sn[0], last_day=sn[1] or (ft[-1] if ft else 0.0),
                    agent_times=a["times"], agent_records=a["n"],
                    running=RUNNER.running() and RUNNER.run_name == run)

    def series(self, run, q):
        def read(S):
            n = len(S.dimensions["time"])
            if n == 0:
                return {"n": 0, "time": [], "vars": {}}
            if int(q.get("have", -1)) == n:
                return {"n": n, "unchanged": True}
            want = [v for v in (q.get("vars") or "biomass,agents,eggs,plankton,budget_error,ibm_cells,catch").split(",")
                    if v in S.variables and v != "time"]
            maxpts = max(50, int(q.get("maxpts", 900)))
            idx = np.arange(n) if n <= maxpts else np.unique(np.r_[np.linspace(0, n - 1, maxpts).astype(int), n - 1])
            out = {"n": n, "time": [round(float(x), 4) for x in np.asarray(S["time"][:])[idx]], "vars": {}}
            for v in want:
                a = np.asarray(S[v][:])[idx]
                out["vars"][v] = np.round(np.where(np.isfinite(a), a, 0.0), 8).tolist()
            return out
        return STORE.read(run, "series.nc", read, {"error": "no series yet"})

    def field(self, run, q):
        def read(F):
            nt = len(F.dimensions["time"])
            if not nt:
                return None
            t = max(0, min(int(q.get("t", nt - 1)), nt - 1))
            a = field_slice(F, q.get("var", "fish_biomass"), t,
                            parse_species(q.get("species", -1)), int(q.get("stage", -1)),
                            int(q.get("depth", 0)), int(q.get("trait", 0)))
            finite = a[np.isfinite(a)]
            return (np.ascontiguousarray(a, dtype="<f4").tobytes(),
                    dict(t=t, day=float(F["time"][t]), shape=list(a.shape),
                         min=float(finite.min()) if finite.size else 0.0,
                         max=float(finite.max()) if finite.size else 0.0))
        got = STORE.read(run, "fields.nc", read)
        if got is None:
            return self.send_error_json(404, "no fields yet")
        return self.send_bin(*got)

    def vectors(self, run, q):
        def read(F):
            nt = len(F.dimensions["time"])
            if not nt:
                return None
            t = max(0, min(int(q.get("t", nt - 1)), nt - 1))
            k, s = int(q.get("depth", 0)), max(1, int(q.get("stride", 2)))
            u = np.asarray(F["u"][t, k], dtype="f4")[::s, ::s]
            v = np.asarray(F["v"][t, k], dtype="f4")[::s, ::s]
            return (np.ascontiguousarray(np.stack([u, v]), dtype="<f4").tobytes(),
                    dict(t=t, stride=s, shape=list(u.shape),
                         max=float(np.nanmax(np.hypot(u, v))) if np.isfinite(u).any() else 0.0))
        got = STORE.read(run, "fields.nc", read)
        if got is None:
            return self.send_error_json(404, "no fields yet")
        return self.send_bin(*got)

    def agents(self, run, q):
        idx = STORE.agent_index(run)
        if not idx["times"]:
            return self.send_bin(b"", dict(n=0, day=None, k=-1))

        def read(A):
            nrec = len(A.dimensions["rec"])
            k = max(0, min(int(q.get("k", len(idx["times"]) - 1)), len(idx["times"]) - 1))
            i0 = min(idx["start"][k], nrec)
            i1 = min(idx["start"][k + 1] if k + 1 < len(idx["start"]) else nrec, nrec)
            limit = max(500, int(q.get("max", 25000)))
            step = max(1, -(-(i1 - i0) // limit))
            sl = slice(i0, i1, step)
            cols = [np.asarray(A[c][sl], dtype="<f4") for c in ("lon", "lat", "n", "species", "stage", "depth", "L")]
            return (np.concatenate(cols).astype("<f4").tobytes(),
                    dict(n=int(cols[0].size), day=idx["times"][k], k=k, total=int(i1 - i0), subsample=step,
                         fields=["lon", "lat", "n", "species", "stage", "depth", "L"]))
        got = STORE.read(run, "agents.nc", read)
        return self.send_bin(*(got or (b"", dict(n=0, day=None, k=-1))))

    def validation(self, run, q):
        """Model-observation comparison for every species with observations (obsval.py)."""
        try:
            v = validation_data(run, float(q.get("window", 365)))
        except Exception as e:                                    # a missing or unreadable source
            return self.send_error_json(500, f"validation failed: {type(e).__name__}: {e}")
        if v is None:
            return self.send_error_json(404, "no gridded output yet")
        M, P = v["M"], v["P"]
        return self.send_json(finite(dict(
            run=run, window=float(q.get("window", 365)), t0=M["t0"], t1=M["t1"], frames=M["frames"],
            years=P["years"], built=P["built"], global_grid=bool(v["grid"].glob),
            fishing=bool(M["catch"] is not None and (M["catch"] > 0).any()), species=v["rows"])))

    def valmap(self, run, q):
        """Model juvenile+adult biomass and one observation layer for a species, NaN over land."""
        v = validation_data(run, float(q.get("window", 365)))
        if v is None:
            return self.send_error_json(404, "no gridded output yet")
        import obsval
        g, sp = v["grid"], int(q.get("species", 0))
        if not 0 <= sp < len(v["names"]):
            return self.send_error_json(400, "bad species index")
        src = q.get("source", "obis")
        mod = np.where(g.mask, v["M"]["juvad"][sp], np.nan)
        obs = np.where(g.mask, obsval.obs_field(g, v["P"], v["names"][sp], src), np.nan)
        return self.send_bin(np.ascontiguousarray(np.concatenate([mod.ravel(), obs.ravel()]), dtype="<f4").tobytes(),
                             dict(ny=g.ny, nx=g.nx, source=src, species=v["names"][sp]))

    def range(self, run, q):
        """Robust colour limits for a variable, sampled over the run so animations do not flicker."""
        spec = (run, q.get("var"), q.get("species"), q.get("stage"), q.get("depth"), q.get("trait"), q.get("log"))
        key = json.dumps(spec)

        def read(F):
            nt = len(F.dimensions["time"])
            if not nt:
                return None
            bucket = nt // 8
            cached = STORE.ranges.get(key)
            if cached and cached[0] == bucket:
                return cached[1]
            vals = []
            for t in sorted(set(int(x) for x in np.linspace(0, nt - 1, min(nt, 8)))):
                a = field_slice(F, q.get("var", "fish_biomass"), t, parse_species(q.get("species", -1)),
                                int(q.get("stage", -1)), int(q.get("depth", 0)), int(q.get("trait", 0)))
                v = a[np.isfinite(a)]
                vals.append(v[:: max(1, v.size // 4000)])
            lo, hi = percentiles(np.concatenate(vals) if vals else np.zeros(1), log=q.get("log") == "1")
            out = {"lo": lo, "hi": hi, "sampled": len(vals)}
            STORE.ranges[key] = (bucket, out)
            return out
        return STORE.read(run, "fields.nc", read) or {"lo": 0.0, "hi": 1.0}


VAL_CACHE, VAL_LOCK = {}, threading.Lock()


def validation_data(run, window):
    """Time-mean model fields and the comparison rows, recomputed only when fields.nc grows."""
    import obsval
    with VAL_LOCK:
        P = obsval.products()

        def read(F):
            nt = len(F.dimensions["time"])
            hit = VAL_CACHE.get((run, window))
            if hit and hit["nt"] == nt:
                return hit
            if not nt or "fish_biomass" not in F.variables:
                return None
            grid = obsval.Grid(F["lon"][:], F["lat"][:], np.asarray(F["mask"][:]) > 0.5)
            names = [n for n in F.getncattr("species").split(",") if n]
            M = obsval.model_means(F, STORE.open(run, "series.nc"), window)
            rows = obsval.compare(grid, names, M, P)
            yr = obsval.yearly(F, grid, names, P)
            for r in rows:
                if r["name"] in yr and "ram" in r:
                    r["ram"]["yearly"] = yr[r["name"]]
            out = dict(nt=nt, grid=grid, M=M, P=P, names=names, rows=rows)
            VAL_CACHE[(run, window)] = out
            return out
        return STORE.read(run, "fields.nc", read)


def finite(o):
    """NaN and inf are not JSON: send them as null."""
    if isinstance(o, dict):
        return {k: finite(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [finite(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    return o


def jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def main():
    ap = argparse.ArgumentParser(description="fishNET browser dashboard")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--no-browser", action="store_true")
    args = ap.parse_args()
    if not WEB.is_dir():
        sys.exit(f"missing web assets: {WEB}")
    srv = ThreadingHTTPServer((args.host, args.port), Handler)
    srv.daemon_threads = True
    url = f"http://{args.host}:{args.port}/"
    print(f"fishNET dashboard on {url}   (model: {ROOT / 'fishnet.py'}, runs: {RUNS_DEFAULT})")
    print("Ctrl-C to stop")
    if not args.no_browser:
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nstopping...")
        RUNNER.stop()
        srv.shutdown()


if __name__ == "__main__":
    main()
