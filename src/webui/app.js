/* fishNET dashboard — run control, live polling, playback and charts. */
(() => {
const $ = id => document.getElementById(id);
const api = (p, q) => fetch("/api/" + p + (q ? "?" + new URLSearchParams(q) : "")).then(r => r.json());
const post = (p, b) => fetch("/api/" + p, { method: "POST", body: JSON.stringify(b) }).then(r => r.json());

const SERIES_COLORS = () => {
  const s = getComputedStyle(document.body);
  return [1, 2, 3, 4, 5, 6, 7, 8].map(i => s.getPropertyValue(`--series-${i}`).trim());
};
const STAGE_COLORS = () => {                     // one hue, light -> dark: an ordinal ramp for egg->adult
  const dark = document.documentElement.getAttribute("data-theme") === "dark";
  return dark ? ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab"] : ["#9ec5f4", "#5598e7", "#2a78d6", "#184f95"];
};
const CAUSE_COLORS = () => SERIES_COLORS();
/* Eight validated hues, then composite encoding: the 9th-16th series reuse the hues with a dashed line, the
   17th-24th dotted, and so on, so a hue is never cycled bare. A style follows its entity (species index or
   group order), never its rank on screen. */
const DASHES = ["", "6 3", "1.5 3", "8 3 2 3", "3 2 1 2"];
const styleOf = k => ({ color: SERIES_COLORS()[k % 8], dash: DASHES[Math.floor(k / 8) % DASHES.length] });
/* Legend/list glyph: a dot for solid series, a short line sample for dashed ones (the dash is the identity). */
const glyph = (st, cls = "dot") => st.dash
  ? `<svg class="glyph" width="16" height="8" aria-hidden="true"><line x1="1" y1="4" x2="15" y2="4" stroke="${st.color}" stroke-width="2.4" stroke-dasharray="${st.dash}" stroke-linecap="butt"/></svg>`
  : `<i class="${cls}" style="background:${st.color}"></i>`;

// field variables offered on the map
const FIELDS = [
  { key: "fish_biomass", label: "fish biomass (g m⁻²)", species: 1, stage: 1, log: 1, cmap: "blues" },
  { key: "fish_numbers", label: "fish numbers (m⁻²)", species: 1, stage: 1, log: 1, cmap: "blues" },
  { key: "agent_biomass", label: "agent biomass (g m⁻²)", species: 1, log: 1, cmap: "blues" },
  { key: "temp", label: "temperature (°C)", depth: 1, log: 0, cmap: "thermal" },
  { key: "P", label: "phytoplankton (mmol N m⁻³)", depth: 1, log: 1, cmap: "viridis" },
  { key: "Z", label: "meso-zooplankton (mmol N m⁻³)", depth: 1, log: 1, cmap: "viridis" },
  { key: "krill", label: "krill (mmol N m⁻³)", depth: 1, log: 1, cmap: "viridis" },
  { key: "N", label: "nutrients (mmol N m⁻³)", depth: 1, log: 0, cmap: "blues" },
  { key: "D", label: "detritus (mmol N m⁻³)", depth: 1, log: 1, cmap: "inferno" },
  { key: "O2", label: "oxygen (mmol O₂ m⁻³)", depth: 1, log: 0, cmap: "blues" },
  { key: "par", label: "PAR (W m⁻²)", depth: 1, log: 0, cmap: "inferno" },
  { key: "speed", label: "current speed (m s⁻¹)", depth: 1, log: 0, cmap: "viridis" },
  { key: "w", label: "vertical velocity (m d⁻¹, up +)", depth: 1, log: 0, cmap: "thermal" },
  { key: "benthos", label: "benthic fauna (g m⁻²)", log: 1, cmap: "inferno" },
  { key: "trait_mean", label: "mean breeding value", species: 1, trait: 1, log: 0, cmap: "thermal" },
  { key: "ibm", label: "IBM cells", log: 0, cmap: "blues" },
  { key: "bottom_depth", label: "bottom depth (m)", log: 0, cmap: "blues" }
];

const S = {
  run: null, meta: null, times: null, idx: 0, playing: false, follow: true, fps: 10,
  tab: "map", logSince: 0, series: {}, seriesN: -1, seriesVars: "", status: null,
  frames: new Map(), agents: new Map(), vectors: new Map(), ibm: new Map(), ranges: new Map(),
  lastRange: 0, pendingMeta: false,
  hidden: new Set(), groupMode: "none"        // species filter + aggregation for every chart
};
const TAB_VARS = {                                   // biomass/agents are always fetched: they feed the header tiles
  map: "",
  pops: "numbers,eggs,catch,losses",
  eco: "plankton,benthos,oxygen,o2_min,export_sinking,export_fish,export_respired,metabolic_index,ibm_cells,budget_error,behavior",
  traits: "trait_mean,trait_sd,max_gen,fst",
  figs: "",
  val: ""
};
const varsFor = tab => ["biomass", "agents"].concat(TAB_VARS[tab] ? TAB_VARS[tab].split(",") : []).join(",");

// ----------------------------------------------------------------- helpers
/* "all shown species" sends the visible indices so the map sums the same set the charts do.
   A single named species is sent as its own index, unaffected by the filter. */
const speciesParam = () => {
  const v = $("spSelect").value;
  if (v !== "-1") return v;
  const vis = shown();
  return (!vis.length || vis.length === spCount()) ? "-1" : vis.join(",");
};
const spec = () => ({
  var: $("varSelect").value, species: speciesParam(), stage: $("stageSelect").value,
  depth: $("depthSelect").value, trait: $("traitSelect").value
});
const fieldDef = () => FIELDS.find(f => f.key === $("varSelect").value) || FIELDS[0];

/* ---- species filter and functional-group aggregation -------------------------------
   S.hidden holds the species indices the user switched off; S.groupMode is how series are
   combined for display. Everything that draws a per-species series goes through
   seriesFor(), so the filter and the grouping apply everywhere at once. */
const DEFAULT_SHOWN = 5;                    // species visible before the user touches the filter
const spCount = () => (S.meta && S.meta.species ? S.meta.species.length : 0);
const shown = () => Array.from({ length: spCount() }, (_, q) => q).filter(q => !S.hidden.has(q));
const groupOf = q => {
  const g = ((S.meta && S.meta.groups) || [])[q];
  if (g && g !== "other" && g !== "pelagic" && g !== "mammal") return g;
  const sp = ((S.meta && S.meta.species) || [])[q];
  return (S.speciesGroups && sp && S.speciesGroups[sp]) || g || "other";
};

const FEISTY_NAMES = { forage: "forage fish", mesopelagic: "mesopelagic fish", large_pelagic: "large pelagics",
                       midwater_predator: "midwater predators", demersal: "demersal fish", none: "not in FEISTY" };
const feistyOf = q => ((S.meta && S.meta.feisty) || [])[q] || "none";
const bucketOf = q => S.groupMode === "feisty" ? feistyOf(q) : groupOf(q);

function groupList() {                      // groups in the order their first species appears (all species, so a
  const seen = [];                          // group keeps its style while the filter changes)
  for (let q = 0; q < spCount(); q++) { const g = bucketOf(q); if (!seen.includes(g)) seen.push(g); }
  if (S.groupMode === "feisty")             // FEISTY's own order, with what it leaves out last
    seen.sort((a, b) => Object.keys(FEISTY_NAMES).indexOf(a) - Object.keys(FEISTY_NAMES).indexOf(b));
  return seen;
}

/* Buckets to draw: [{name, color, members:[speciesIndex,...]}, ...] honouring filter+mode. */
function buckets() {
  const cols = SERIES_COLORS(), vis = shown();
  if (S.groupMode === "total")
    return vis.length ? [{ name: "all species", color: cols[0], dash: "", members: vis }] : [];
  if (S.groupMode === "group" || S.groupMode === "feisty")
    return groupList().map((g, i) => ({ name: S.groupMode === "feisty" ? FEISTY_NAMES[g] || g : g.replace(/_/g, " "),
                                        ...styleOf(i), members: vis.filter(q => bucketOf(q) === g) }))
                      .filter(b => b.members.length);
  return vis.map(q => ({ name: S.meta.species[q], ...styleOf(q), members: [q] }));
}

/* Per-species dot style matching the map legend (same buckets); hidden species get none. Dashed buckets
   (the 9th and later) are drawn as rings so they stay distinct from the solid dot of the same hue. */
function agentColors() {
  const c = new Array(spCount()).fill(null);
  buckets().forEach(b => b.members.forEach(q => { c[q] = { color: b.color, ring: !!b.dash }; }));
  return c;
}

/* Build chart series from a (time x species) variable, summing within each bucket.
   `pick` pulls one species' value out of a row; `transform` is applied after summing. */
function seriesFor(varName, pick, transform) {
  const v = S.series.vars && S.series.vars[varName];
  if (!v || !S.meta) return [];
  const get = pick || ((row, q) => row[q]);
  return buckets().map(b => ({
    name: b.name, color: b.color, dash: b.dash,
    y: v.map(row => {
      let t = 0;
      for (const q of b.members) { const x = get(row, q); if (isFinite(x)) t += x; }
      return transform ? transform(t) : t;
    })
  }));
}

/* Means, not sums, for intensive quantities (a rate or an index cannot be added up). */
function meanSeriesFor(varName, pick, skip) {
  const v = S.series.vars && S.series.vars[varName];
  if (!v || !S.meta) return [];
  const get = pick || ((row, q) => row[q]);
  return buckets().map(b => ({ ...b, members: b.members.filter(q => !(skip && skip(q))) }))
    .filter(b => b.members.length)
    .map(b => ({
      name: b.name, color: b.color, dash: b.dash,
      y: v.map(row => {
        let t = 0, n = 0;
        for (const q of b.members) { const x = get(row, q); if (isFinite(x)) { t += x; n++; } }
        return n ? t / n : NaN;
      })
    }));
}
const keyOf = (sp, t) => `${sp.var}|${sp.species}|${sp.stage}|${sp.depth}|${sp.trait}|${t}`;

function dateOf(day) {
  if (!S.meta || !S.meta.start) return "";
  const t0 = Date.parse(S.meta.start.replace(" ", "T") + "Z");
  if (!isFinite(t0)) return "";
  return new Date(t0 + day * 864e5).toISOString().slice(0, 16).replace("T", " ");
}
/* Two output streams with different cadences: gridded fields (output_every_hours) and agent
   snapshots (agent_output_every_hours). Playback steps through one of them; the other is matched
   to the nearest available time, so agents hold still between snapshots unless you drive by them. */
const agentTimeline = () => $("timelineSelect").value === "agents" && S.times && S.times.agent_times.length > 1;
const timeline = () => (!S.times ? [] : agentTimeline() ? S.times.agent_times : S.times.field_times);
const nearestIdx = (arr, v) => {
  let best = 0, bd = Infinity;
  for (let i = 0; i < arr.length; i++) { const d = Math.abs(arr[i] - v); if (d < bd) { bd = d; best = i; } }
  return best;
};
const spacing = arr => (arr.length > 1 ? arr[arr.length - 1] - arr[arr.length - 2] : 0);

const trim = (n, d = 1) => (Math.abs(n) >= 1e5 || (n && Math.abs(n) < 1e-2) ? n.toExponential(1) : n.toFixed(d));
const mass = t => t >= 1e9 ? (t / 1e9).toFixed(2) + " Gt" : t >= 1e6 ? (t / 1e6).toFixed(2) + " Mt"
                : t >= 1e3 ? (t / 1e3).toFixed(1) + " kt" : t.toFixed(0) + " t";

async function bin(path, q) {
  const r = await fetch("/api/" + path + "?" + new URLSearchParams(q));
  if (!r.ok) throw new Error((await r.json()).error || r.status);
  const meta = JSON.parse(r.headers.get("X-Meta") || "{}");
  return { meta, data: new Float32Array(await r.arrayBuffer()) };
}

// ----------------------------------------------------------------- run control
async function refreshState() {
  const st = await api("state", { since: S.logSince });
  S.status = st.status;
  if (st.species_groups) S.speciesGroups = st.species_groups;
  renderRunList(st.runs);
  renderLog(st.status);
  renderTiles();
  if (!S.run && st.runs.length) selectRun(st.status.running ? st.status.run : st.runs[0].name);
  if (st.status.running && st.status.run && st.status.run !== S.run && S.autoFollowRun) selectRun(st.status.run);
  if (!$("namelistSelect").options.length) {
    $("namelistSelect").innerHTML = st.namelists.map(n => `<option>${n}</option>`).join("");
    if (st.namelists.includes("namelist.toml")) $("namelistSelect").value = "namelist.toml";
    loadNamelist();
  }
}

function renderRunList(runs) {
  const sel = $("runSelect"), cur = S.run;
  sel.innerHTML = runs.map(r => {
    const tag = r.running ? " ● running" : "";
    return `<option value="${r.name}"${r.name === cur ? " selected" : ""}>${r.name}${tag} — day ${(r.day || 0).toFixed(0)}</option>`;
  }).join("") || `<option value="">no runs yet</option>`;
  S.runsInfo = runs;
  const info = runs.find(r => r.name === cur);
  $("liveBadge").classList.toggle("hidden", !(info && info.running));
  if (S.tab === "figs") renderFigures();
}

function renderLog(status) {
  const box = $("log");
  if (status.first > S.logSince && S.logSince) box.textContent = "";      // log rolled over
  if (status.lines && status.lines.length) {
    const atBottom = box.scrollTop + box.clientHeight >= box.scrollHeight - 30;
    box.textContent += (box.textContent ? "\n" : "") + status.lines.join("\n");
    const lines = box.textContent.split("\n");
    if (lines.length > 1200) box.textContent = lines.slice(-1000).join("\n");
    if (atBottom) box.scrollTop = box.scrollHeight;
  }
  S.logSince = status.next;
  const p = status.progress || {};
  const frac = p.nsteps ? p.step / p.nsteps : (status.running ? 0 : 0);
  $("progressBar").style.width = (frac * 100).toFixed(1) + "%";
  $("startBtn").disabled = status.running;
  const sel = (S.runsInfo || []).find(r => r.name === S.run);
  $("resumeBtn").disabled = status.running || !(sel && sel.restart_day != null);
  $("stopBtn").disabled = !status.running;
  // The live table follows the same filter as the charts, so a 23-species run stays readable.
  const all = Object.keys(status.species || {});
  const visNames = (S.meta && S.meta.species) ? shown().map(q => S.meta.species[q]) : null;
  const names = visNames ? all.filter(n => visNames.includes(n)) : all;
  $("speciesNow").innerHTML = status.running && names.length
    ? `<div class="sn-row sn-head"><span>latest step</span><span>biomass</span><span>adult</span><span>agents</span></div>` +
      names.map(n => {
        const i = all.indexOf(n);
        const s = status.species[n], adult = s.stages && s.stages.length === 4 ? (s.stages[3] * 100).toFixed(0) + "%" : "—";
        return `<div class="sn-row"><span>${glyph(styleOf(i))}${n}</span>` +
               `<span>${mass(s.biomass_t)}</span><span>${adult}</span><span>${s.agents.toLocaleString()}</span></div>`;
      }).join("")
    : "";
  $("progressText").textContent = status.running
    ? `${p.step || 0}/${p.nsteps || "?"} steps · day ${trim(p.day || 0)} · ${trim(p.s_per_step || 0, 2)} s/step · ETA ${trim(p.eta_min || 0)} min`
    : (status.run ? `${status.run}: ${status.returncode === 0 ? "finished" : status.returncode === null ? "idle" : "exited (" + status.returncode + ")"}` : "idle");
}

function renderTiles() {
  const st = S.status || {}, p = st.progress || {};
  const B = S.series.vars && S.series.vars.biomass;
  let total = null, agents = p.agents;
  if (B && B.length) total = B[B.length - 1].flat().reduce((a, b) => a + b, 0) / 1e6;
  if (!st.running && S.series.vars && S.series.vars.agents) {
    const a = S.series.vars.agents; agents = a[a.length - 1].reduce((x, y) => x + y, 0);
  }
  const tl = timeline();
  const day = tl.length ? tl[Math.min(S.idx, tl.length - 1)] : p.day;
  const tiles = [
    ["model day", day === undefined ? "—" : trim(day)],
    ["date", day === undefined ? "—" : dateOf(day).slice(0, 10)],
    ["total biomass", total === null ? "—" : mass(total)],
    ["agents", agents === undefined ? "—" : Math.round(agents).toLocaleString()],
    ["IBM cells", p.ibm_cells === undefined ? "—" : p.ibm_cells],
    ["speed", p.s_per_step ? trim(p.s_per_step, 2) + " s/step" : "—"]
  ];
  $("tiles").innerHTML = tiles.map(([k, v]) => `<div class="tile"><div class="k">${k}</div><div class="v">${v}</div></div>`).join("");
}

async function loadNamelist() {
  const name = $("namelistSelect").value;
  const r = await api("namelist", { name });
  if (r.error) return;
  $("namelistText").value = r.text;
  const cfg = parseNamelist(r.text);
  $("optDays").value = cfg.run.days ?? 365;
  $("optDt").value = cfg.run.dt_hours ?? 12;
  $("optSeed").value = cfg.run.seed ?? 1;
  $("optOut").value = cfg.run.output_every_hours ?? 48;
  $("optAgents").value = cfg.run.agent_output_every_hours ?? 240;
  $("optLog").value = cfg.run.log_every ?? 1;
  $("optNx").value = cfg.grid.nx ?? 48;
  $("optNy").value = cfg.grid.ny ?? 32;
  $("optWarm").value = cfg.ocean.warming_per_year ?? 0;
  $("optSpin").value = cfg.run.spinup_days ?? 30;
  $("runName").value = (cfg.run.name || "run") + "_gui";
  const cols = SERIES_COLORS();
  const byGroup = new Map();
  cfg.species.forEach((s, i) => {
    const grp = (S.speciesGroups && S.speciesGroups[s.name]) || s.group || "other";
    if (!byGroup.has(grp)) byGroup.set(grp, []);
    byGroup.get(grp).push({ s, i });
  });
  let listHtml = `<div class="species-row"><span class="hd">species</span><span class="hd">B₀ g/m²</span><span class="hd">F /yr</span></div>`;
  for (const [grp, items] of byGroup.entries()) {
    if (byGroup.size > 1) {
      listHtml += `<div class="grp" style="font-size:10px;text-transform:uppercase;letter-spacing:.05em;color:var(--accent);margin:6px 0 2px;font-weight:600;">${grp.replace(/_/g, " ")}</div>`;
    }
    for (const { s, i } of items) {
      listHtml += `<div class="species-row">
         <span>${glyph(styleOf(i))}${s.name}</span>
         <input type="number" step="0.1" min="0" data-sp="${i}" data-key="biomass" value="${s.biomass ?? 1}">
         <input type="number" step="0.05" min="0" data-sp="${i}" data-key="fishing_F" value="${s.fishing_F ?? 0}">
       </div>`;
    }
  }
  $("speciesList").innerHTML = listHtml;
}

/* Just enough TOML for the launch form: scalars in [run]/[grid]/[ocean] and the [[species]] list. */
function parseNamelist(text) {
  const cfg = { run: {}, grid: {}, ocean: {}, species: [] };
  let sec = null, cur = null;
  for (const raw of text.split("\n")) {
    const line = raw.split("#")[0].trim();
    if (!line) continue;
    if (line === "[[species]]") { cur = {}; cfg.species.push(cur); sec = "species"; continue; }
    if (line.startsWith("[")) { sec = line.replace(/[\[\]]/g, ""); cur = cfg[sec] || null; continue; }
    const m = line.match(/^([A-Za-z_][\w]*)\s*=\s*(.+)$/);
    if (!m || !cur) continue;
    let v = m[2].trim();
    if (/^-?[\d.]+(e-?\d+)?$/i.test(v)) v = parseFloat(v);
    else if (v === "true" || v === "false") v = v === "true";
    else v = v.replace(/^"|"$/g, "");
    cur[m[1]] = v;
  }
  cfg.species.forEach(s => { s.name = String(s.file || "species").split("/").pop().replace(".toml", ""); });
  return cfg;
}

async function startRun() {
  const overrides = {
    "run.days": +$("optDays").value, "run.dt_hours": +$("optDt").value, "run.seed": +$("optSeed").value,
    "run.output_every_hours": +$("optOut").value, "run.agent_output_every_hours": +$("optAgents").value,
    "run.log_every": +$("optLog").value, "grid.nx": +$("optNx").value, "grid.ny": +$("optNy").value,
    "ocean.warming_per_year": +$("optWarm").value, "run.spinup_days": +$("optSpin").value
  };
  $("speciesList").querySelectorAll("input").forEach(inp => {
    overrides[`species.${inp.dataset.sp}.${inp.dataset.key}`] = +inp.value;
  });
  const body = { name: $("runName").value || "dashboard_run", text: $("namelistText").value, overrides, fresh: true };
  const r = await post("start", body);
  if (r.error) return alert("could not start: " + r.error);
  S.logSince = 0; $("log").textContent = "";
  S.autoFollowRun = true;
  resetRun();
  S.run = body.name;
  setTimeout(() => selectRun(body.name, true), 1200);
}

/* Continue the run chosen at the top from its newest restart file; asks how far to take it. */
async function resumeRun() {
  const r = (S.runsInfo || []).find(x => x.name === S.run);
  if (!r || r.restart_day == null) return alert("the selected run has no restart files to resume from");
  const ans = prompt(`Resume '${r.name}' from its restart at day ${r.restart_day}.\nRun until day:`, r.cfg_days ?? "");
  if (ans === null) return;
  const days = parseInt(ans, 10);
  if (!(days > r.restart_day)) return alert(`the end day must be past the restart day (${r.restart_day})`);
  const res = await post("resume", { name: r.name, days });
  if (res.error) return alert("could not resume: " + res.error);
  S.logSince = 0; $("log").textContent = "";
  S.autoFollowRun = true;
  resetRun();
  S.run = r.name;
  setTimeout(() => selectRun(r.name, true), 1200);
}

// ----------------------------------------------------------------- run selection & metadata
function resetRun() {
  S.meta = null; S.times = null; S.idx = 0; S.series = {}; S.seriesN = -1;
  S.frames.clear(); S.agents.clear(); S.vectors.clear(); S.ibm.clear(); S.ranges.clear();
}

async function selectRun(name, keepLog) {
  if (!name) return;
  if (name !== S.run) { resetRun(); S.fieldSig = null; S.run = name; if (!keepLog) S.autoFollowRun = false; Validation.reset(); }
  $("runSelect").value = name;
  await loadMeta();
  if (S.tab === "val") Validation.show(S.run);
}

async function refreshFields() {                // a run's field list, re-read without resetting the view
  const m = await api("meta", { run: S.run });
  if (m.error || !S.meta) return;
  S.meta.fields = m.fields;
  buildVarSelect(S.meta);
}

async function loadMeta() {
  const m = await api("meta", { run: S.run });
  if (m.error) { S.meta = null; Charts.empty($("miniBiomass"), "waiting for output…"); return; }
  const first = !S.meta;
  S.meta = m;
  MapView.setMeta(m);
  if (first) buildSelectors(m); else buildVarSelect(m);
  S.times = { field_times: m.field_times, agent_times: m.agent_times, n_series: m.n_series, last_day: m.last_day, running: m.running };
  S.idx = Math.max(0, timeline().length - 1);
  $("timeSlider").max = Math.max(0, timeline().length - 1);
  $("timeSlider").value = S.idx;
  await Promise.all([drawMap(), loadSeries(true)]);
}

function renderFilter() {
  const m = S.meta, sp = (m && m.species) || [], cols = SERIES_COLORS();
  if (!sp.length) { $("spFilter").innerHTML = ""; $("spCount").textContent = ""; return; }
  const byGroup = new Map();
  sp.forEach((n, q) => {
    const g = S.groupMode === "feisty" ? FEISTY_NAMES[feistyOf(q)] || feistyOf(q) : groupOf(q);
    if (!byGroup.has(g)) byGroup.set(g, []);
    byGroup.get(g).push(q);
  });
  const single = byGroup.size <= 1;          // no point labelling groups when there is only one
  $("spFilter").innerHTML = [...byGroup.entries()].map(([g, qs]) =>
    (single ? "" : `<div class="grp">${g.replace(/_/g, " ")}</div>`) +
    qs.map(q => `<label title="${sp[q]}">
        <input type="checkbox" data-sp="${q}" ${S.hidden.has(q) ? "" : "checked"}>
        ${glyph(styleOf(q))}${sp[q]}${(m.individual || [])[q]
          ? ` <span class="tag" title="always individuals: agents everywhere, never biomass">ind.</span>` : ""}</label>`).join("")
  ).join("");
  $("spCount").textContent = `${shown().length} of ${sp.length}`;
  $("spFilter").querySelectorAll("input").forEach(inp => inp.onchange = () => {
    const q = +inp.dataset.sp;
    inp.checked ? S.hidden.delete(q) : S.hidden.add(q);
    afterFilterChange();
  });
}

/* A filter or grouping change only affects what is drawn, so redraw rather than refetch. */
function afterFilterChange() {
  $("spCount").textContent = `${shown().length} of ${spCount()}`;
  refreshSpeciesSelects();
  renderCharts();
  if (S.tab === "map") drawMap();
}

/* The map and per-species panels list only the species still switched on. */
function refreshSpeciesSelects() {
  const sp = (S.meta && S.meta.species) || [], vis = shown();
  const keepSp = $("spSelect").value, keepPop = $("popSpecies").value;
  $("spSelect").innerHTML = `<option value="-1">all shown species</option>` +
    vis.map(q => `<option value="${q}">${sp[q]}</option>`).join("");
  $("popSpecies").innerHTML = vis.map(q => `<option value="${q}">${sp[q]}</option>`).join("");
  // keep the previous choice when it is still visible, else fall back
  // ("" is the value before any option existed, and +"" is 0, so it must not count as species 0)
  $("spSelect").value = (keepSp === "-1" || (keepSp !== "" && vis.includes(+keepSp))) ? keepSp : "-1";
  $("popSpecies").value = (keepPop !== "" && vis.includes(+keepPop)) ? keepPop : (vis.length ? String(vis[0]) : "");
  $("spLegend").innerHTML = buckets().map(b =>
    `<div class="row">${b.dash ? glyph(b) : `<i class="dot" style="background:${b.color}"></i>`}${b.name}</div>`).join("");
}

/* The map's field list follows what fields.nc holds. It is rebuilt whenever the run's metadata is reloaded,
   keeping the current choice: a file opened during a spin-up may not list its fields yet. */
function buildVarSelect(m) {
  const keep = $("varSelect").value;
  const avail = FIELDS.filter(f => !m.fields || f.key === "speed" || f.key === "bottom_depth" || m.fields.includes(f.key));
  const sig = avail.map(f => f.key).join(",");
  if (sig === S.fieldSig) return;
  S.fieldSig = sig;
  $("varSelect").innerHTML = avail
    .map(f => `<option value="${f.key}">${f.label}</option>`).join("");
  if (avail.some(f => f.key === keep)) $("varSelect").value = keep;
  else { const f = fieldDef(); $("logScale").checked = !!f.log; $("cmapSelect").value = f.cmap; }
  syncVarControls();
}

function buildSelectors(m) {
  const sp = m.species || [], colors = SERIES_COLORS();
  // A new run starts with only the first few species on: 23 lines is an unreadable legend.
  // The filter panel is right there to add the rest, or "all" to switch every species on.
  S.hidden = new Set(sp.map((n, q) => q).slice(DEFAULT_SHOWN));
  buildVarSelect(m);
  $("stageSelect").innerHTML = `<option value="-1">all stages</option>` + (m.stages || []).map((n, i) => `<option value="${i}">${n}</option>`).join("");
  $("depthSelect").innerHTML = (m.depth || []).map((z, i) => `<option value="${i}">${z.toFixed(0)} m</option>`).join("");
  $("traitSelect").innerHTML = (m.traits || []).map((n, i) => `<option value="${i}">${n}</option>`).join("");
  $("traitChartSelect").innerHTML = (m.traits || []).map((n, i) => `<option value="${i}">${n}</option>`).join("");
  renderFilter();
  refreshSpeciesSelects();
  syncVarControls();
}

function syncVarControls() {
  const f = fieldDef();
  $("speciesWrap").classList.toggle("hidden", !f.species);
  $("stageWrap").classList.toggle("hidden", !f.stage);
  $("depthWrap").classList.toggle("hidden", !f.depth);
  $("traitWrap").classList.toggle("hidden", !f.trait);
  $("spLegend").classList.toggle("hidden", !$("showAgents").checked);
}

// ----------------------------------------------------------------- map frames
async function frameAt(t) {
  const sp = spec(), k = keyOf(sp, t);
  if (S.frames.has(k)) return S.frames.get(k);
  const p = bin("field", { run: S.run, t, ...sp }).then(r => r.data).catch(() => null);
  S.frames.set(k, p);
  if (S.frames.size > 240) S.frames.delete(S.frames.keys().next().value);
  return p;
}

async function colorRange() {
  const sp = spec(), log = $("logScale").checked ? 1 : 0;
  if (!$("fixedScale").checked) return null;
  const k = keyOf(sp, "range") + "|" + log;
  const fresh = S.times && S.times.running && Date.now() - S.lastRange > 20000;
  if (!S.ranges.has(k) || fresh) {
    S.lastRange = Date.now();
    S.ranges.set(k, api("range", { run: S.run, log, ...sp }));
  }
  return S.ranges.get(k);
}

async function agentsAt(day) {
  if (!$("showAgents").checked || !S.times || !S.times.agent_times.length) return null;
  const k = nearestIdx(S.times.agent_times, day);
  if (!S.agents.has(k)) {
    S.agents.set(k, bin("agents", { run: S.run, k, max: 15000 }).then(r => {
      const n = r.meta.n, d = r.data;
      return { n, day: r.meta.day, lon: d.subarray(0, n), lat: d.subarray(n, 2 * n), num: d.subarray(2 * n, 3 * n),
               species: d.subarray(3 * n, 4 * n), stage: d.subarray(4 * n, 5 * n), total: r.meta.total };
    }).catch(() => null));
    if (S.agents.size > 120) S.agents.delete(S.agents.keys().next().value);
  }
  return S.agents.get(k);
}

async function vectorsAt(t) {
  if (!$("showCurrents").checked) return null;
  const k = t + "|" + $("depthSelect").value;
  if (!S.vectors.has(k)) {
    S.vectors.set(k, bin("vectors", { run: S.run, t, depth: $("depthSelect").value, stride: 2 }).then(r => {
      const [ny, nx] = r.meta.shape, half = ny * nx;
      return { ny, nx, stride: r.meta.stride, scale: r.meta.max || 1, u: r.data.subarray(0, half), w: r.data.subarray(half) };
    }).catch(() => null));
    if (S.vectors.size > 120) S.vectors.delete(S.vectors.keys().next().value);
  }
  return S.vectors.get(k);
}

async function ibmAt(t) {
  if (!$("showIbm").checked) return null;
  if (!S.ibm.has(t)) {
    S.ibm.set(t, bin("field", { run: S.run, t, var: "ibm" }).then(r => r.data).catch(() => null));
    if (S.ibm.size > 160) S.ibm.delete(S.ibm.keys().next().value);
  }
  return S.ibm.get(t);
}

let drawing = false, redrawQueued = false;
async function drawMap() {
  if (!S.meta || !S.times || !S.times.field_times.length) return;
  if (drawing) { redrawQueued = true; return; }
  drawing = true;
  try {
    const times = timeline();
    if (!times.length) return;
    const i = Math.min(S.idx, times.length - 1), day = times[i];
    const t = agentTimeline() ? nearestIdx(S.times.field_times, day) : i;   // the field frame to draw under it
    const [a, rng, ag, vec, ibm] = await Promise.all([frameAt(t), colorRange(), agentsAt(day), vectorsAt(t), ibmAt(t)]);
    if (!a) return;
    const log = $("logScale").checked;
    let lo, hi;
    if (rng) { lo = rng.lo; hi = rng.hi; }
    else {
      let mn = Infinity, mx = -Infinity;
      for (const v of a) if (isFinite(v) && (!log || v > 0)) { if (v < mn) mn = v; if (v > mx) mx = v; }
      lo = isFinite(mn) ? mn : 0; hi = isFinite(mx) ? mx : 1;
    }
    if (log) { lo = Math.max(lo, hi / 1e6, 1e-12); }
    if (!(hi > lo)) hi = lo + Math.abs(lo || 1) * 1e-3;
    MapView.render(a, {
      lo, hi, log, cmap: $("cmapSelect").value, label: fieldDef().label,
      agents: ag && ag.n ? ag : null, vectors: vec, ibm, colors: agentColors(), smooth: $("varSelect").value !== "ibm"
    });
    S.lastDay = day;
    const stale = ag && Math.abs(ag.day - day) > 1e-6;
    $("timeLabel").textContent = `day ${day.toFixed(2)}  ·  ${dateOf(day)}` +
      (ag ? `  ·  ${(ag.total || 0).toLocaleString()} agents${stale ? ` @ day ${ag.day.toFixed(0)}` : ""}` : "");
    $("timeSlider").value = i;
    playerNote(ag);
    renderTiles();
    if (S.tab === "map") renderMapCharts();
    prefetch(i);
  } finally {
    drawing = false;
    if (redrawQueued) { redrawQueued = false; drawMap(); }
  }
}

/* Agent snapshots are usually written far less often than fields, so playing along the field
   timeline redraws the same dots for several frames. Say so, and offer the other timeline. */
function playerNote(ag) {
  const df = spacing(S.times.field_times), da = spacing(S.times.agent_times);
  const note = $("playerNote");
  if (!ag || !da || da <= df * 1.5) { note.textContent = ""; return; }
  note.textContent = agentTimeline()
    ? `stepping agent snapshots (every ${trim(da)} d); the field under them is the nearest frame (every ${trim(df)} d).`
    : `agent snapshots are every ${trim(da)} d but fields every ${trim(df)} d, so the dots hold still for `
      + `${Math.round(da / df)} frames at a time — switch timeline to “agent snapshots”, or set `
      + `agents out = fields out when launching for smooth motion.`;
}

function prefetch(i) {
  if (agentTimeline()) return;
  const n = S.times.field_times.length;
  for (let k = 1; k <= 6; k++) if (i + k < n) frameAt(i + k);
}

// ----------------------------------------------------------------- playback
let lastTick = 0, acc = 0;
function tick(ts) {
  requestAnimationFrame(tick);
  const dt = (ts - lastTick) / 1000; lastTick = ts;
  if (!S.playing || !timeline().length) return;
  acc += Math.min(dt, 0.25);
  const period = 1 / S.fps;
  if (acc < period) return;
  acc = 0;
  const n = timeline().length;
  if (S.idx >= n - 1) {
    if (S.times.running && S.follow) return;      // sit at the leading edge, waiting for new output
    S.idx = 0;
  } else S.idx++;
  drawMap();
}
requestAnimationFrame(ts => { lastTick = ts; requestAnimationFrame(tick); });

function setPlaying(p) {
  S.playing = p;
  $("playBtn").textContent = p ? "❚❚" : "▶";
}

// ----------------------------------------------------------------- series + charts
async function loadSeries(force) {
  if (!S.run) return;
  const vars = varsFor(S.tab);
  const r = await api("series", { run: S.run, vars, maxpts: 900, have: (force || vars !== S.seriesVars) ? -1 : S.seriesN });
  if (r.error) return;
  if (!r.unchanged) { S.series = r; S.seriesN = r.n; S.seriesVars = vars; }
  renderCharts();
  renderTiles();
}

/* Extensive quantities (eggs, catch, agents) add up across a group; intensive ones
   (an index, a rate, a fraction) are averaged instead - summing them would be meaningless. */
const INTENSIVE = new Set(["metabolic_index", "fst"]);

function speciesSeries(name, transform) {
  if (name === "max_gen") return maxSeriesFor(name);     // a maximum stays a maximum across a group
  return INTENSIVE.has(name) ? meanSeriesFor(name) : seriesFor(name, null, transform);
}

function maxSeriesFor(varName) {
  const v = S.series.vars && S.series.vars[varName];
  if (!v || !S.meta) return [];
  return buckets().map(b => ({
    name: b.name, color: b.color, dash: b.dash,
    y: v.map(row => Math.max(...b.members.map(q => row[q]).filter(isFinite)))
  }));
}

function renderCharts() {
  if (!S.meta || !S.series.time) return;
  const now = timeline().length ? timeline()[Math.min(S.idx, timeline().length - 1)] : null;
  if (S.tab === "map") return renderMapCharts();
  const x = S.series.time;
  if (S.tab === "pops") {
    if (!shown().length) {                    // everything filtered out: draw nothing rather than a stale species
      ["chBiomass", "chStages", "chEggs", "chCatch", "chLosses"].forEach(id => Charts.draw($(id), { x, now, series: [] }));
      $("stageSpName").textContent = $("lossSpName").textContent = "(no species shown)";
      return;
    }
    const q = +$("popSpecies").value || shown()[0], sp = S.meta.species[q] || "";
    $("stageSpName").textContent = sp; $("lossSpName").textContent = sp;
    Charts.draw($("chBiomass"), { x, now, log: true, unit: "t",
      series: seriesFor("biomass", (row, q) => row[q].reduce((a, b) => a + b, 0), v => v / 1e6) });
    Charts.draw($("chStages"), { x, now, stacked: true, unit: "t",
      series: (S.meta.stages || []).map((st, k) => ({ name: st, color: STAGE_COLORS()[k], y: S.series.vars.biomass.map(r => r[q][k] / 1e6) })) });
    Charts.draw($("chEggs"), { x, now, log: true, series: speciesSeries("eggs") });
    Charts.draw($("chCatch"), { x, now, series: speciesSeries("catch", v => v / 1e6) });
    const causes = S.meta.causes || [];
    Charts.draw($("chLosses"), { x, now, stacked: true, unit: "t d⁻¹",
      series: causes.map((c, k) => ({ name: c, color: CAUSE_COLORS()[k % 8], y: S.series.vars.losses.map(r => r[q][k] / 1e6) }))
                    .filter((s, k) => causes[k] !== "alive" && s.y.some(v => v > 0)) });
  } else if (S.tab === "eco") {
    const cols = SERIES_COLORS();
    const npzd = ["N", "P", "Z", "K", "D"].map((n, k) => ({ name: n, color: cols[k], y: S.series.vars.plankton.map(r => r[k]) }));
    if (S.series.vars.benthos && S.series.vars.benthos.some(v => v > 0))      // benthic fauna pool, g -> mmol N
      npzd.push({ name: "benthos", color: cols[5], y: S.series.vars.benthos.map(v => v * 1.8) });
    Charts.draw($("chNPZD"), { x, now, log: true, series: npzd });
    Charts.draw($("chO2"), { x, now, series: [
      { name: "domain minimum", color: cols[0], y: S.series.vars.o2_min },
    ] });
    Charts.draw($("chExport"), { x, now, series: [
      { name: "sinking", color: cols[0], y: S.series.vars.export_sinking },
      { name: "fish faeces/carcasses", color: cols[1], y: S.series.vars.export_fish },
      { name: "fish respiration", color: cols[2], y: S.series.vars.export_respired }] });
    // air breathers (whales, dolphins) have no metabolic index: leave them out rather than plot a constant
    Charts.draw($("chPhi"), { x, now, series: meanSeriesFor("metabolic_index", null, q => (S.meta.air_breathing || [])[q]),
                              emptyMsg: "only air-breathing species shown" });
    Charts.draw($("chIbm"), { x, now, series: [{ name: "IBM cells", color: cols[0], y: S.series.vars.ibm_cells }] });
    Charts.draw($("chBudget"), { x, now, zeroLine: true, series: [{ name: "N budget error", color: cols[7], y: S.series.vars.budget_error }] });
    const q = +$("popSpecies").value || 0;
    $("behSpName").textContent = S.meta.species[q] || "";
    Charts.draw($("chBeh"), { x, now, stacked: true,
      series: (S.meta.behaviors || []).map((b, k) => ({ name: b, color: SERIES_COLORS()[k % 8], y: S.series.vars.behavior.map(r => r[q][k]) })) });
  } else if (S.tab === "traits") {
    const k = +$("traitChartSelect").value || 0, cols = SERIES_COLORS();
    $("traitName").textContent = (S.meta.traits || [])[k] || "";
    const tm = S.series.vars.trait_mean, ts = S.series.vars.trait_sd;
    Charts.draw($("chTrait"), { x, now, zeroLine: true,
      series: buckets().map(b => {
        const avg = (rows, f) => rows.map(r => {
          let t = 0, n = 0;
          for (const q of b.members) { const x = f(r, q); if (isFinite(x)) { t += x; n++; } }
          return n ? t / n : NaN;
        });
        const y = avg(tm, (r, q) => r[q][k] - tm[0][q][k]);
        const sd = ts ? avg(ts, (r, q) => r[q][k]) : null;
        return { name: b.name, color: b.color, dash: b.dash, y,
                 band: sd ? [y.map((v, i) => v - sd[i]), y.map((v, i) => v + sd[i])] : null };
      }) });
    Charts.draw($("chGen"), { x, now, series: speciesSeries("max_gen") });
    Charts.draw($("chFst"), { x, now, series: speciesSeries("fst") });
  }
}

function renderMapCharts() {
  if (!S.series.time || !S.meta) return;
  const x = S.series.time;
  const now = timeline().length ? timeline()[Math.min(S.idx, timeline().length - 1)] : null;
  // the small charts under the map drop their legend past 10 series: the map's legend names them, and the
  // tooltip does too, while a 39-row legend would leave no room for the plot
  const legend = buckets().length <= 10;
  if (S.series.vars.biomass) {
    Charts.draw($("miniBiomass"), { x, now, log: true, unit: "t", legend, emptyMsg: shown().length ? "no data yet" : "no species shown",
      series: seriesFor("biomass", (row, q) => row[q].reduce((a, b) => a + b, 0), v => v / 1e6) });
  }
  if (S.series.vars.agents) {
    Charts.draw($("miniAgents"), { x, now, legend, series: seriesFor("agents"), emptyMsg: shown().length ? "no data yet" : "no species shown" });
  }
}

// ----------------------------------------------------------------- figures tab
function renderFigures() {
  const info = (S.runsInfo || []).find(r => r.name === S.run);
  const figs = (info && info.figures) || [];
  $("figGrid").innerHTML = figs.length ? figs.map(f => {
    const url = `/api/figure?run=${encodeURIComponent(S.run)}&name=${encodeURIComponent(f)}`;
    return `<figure><figcaption>${f}</figcaption>` +
      (f.endsWith(".mp4") ? `<video src="${url}" controls loop muted></video>` : `<img src="${url}" loading="lazy">`) +
      `</figure>`;
  }).join("") : `<div class="empty">no figures yet — press “Generate figures”</div>`;
}

// ----------------------------------------------------------------- live polling
async function poll() {
  try {
    await refreshState();
    if (S.run) {
      const t = await api("times", { run: S.run });
      if (!t.error) {
        const grew = !S.times || t.field_times.length !== S.times.field_times.length || t.n_series !== S.times.n_series;
        S.times = t;
        $("timeSlider").max = Math.max(0, timeline().length - 1);
        if (grew) {
          if (!S.meta || !S.meta.species) await loadMeta();
          else if (!S.meta.fields || !S.meta.fields.includes("fish_biomass")) await refreshFields();
          if (S.follow && t.running) { S.idx = Math.max(0, timeline().length - 1); }
          await drawMap();
          await loadSeries(false);
        }
        if (S.tab === "val") Validation.tick(S.run, t);
      }
    }
  } catch (e) { /* server restarting or run switching */ }
}

// ----------------------------------------------------------------- wiring
function on(id, ev, fn) { const e = $(id); if (e) e.addEventListener(ev, fn); }

on("runSelect", "change", e => selectRun(e.target.value));
on("namelistSelect", "change", loadNamelist);
on("startBtn", "click", startRun);
on("stopBtn", "click", () => post("stop", {}));
on("resumeBtn", "click", resumeRun);
on("playBtn", "click", () => setPlaying(!S.playing));
on("timeSlider", "input", e => { S.idx = +e.target.value; $("followLive").checked = false; S.follow = false; drawMap(); });
on("fpsSelect", "change", e => { S.fps = +e.target.value; });
on("followLive", "change", e => { S.follow = e.target.checked; if (S.follow && S.times) { S.idx = timeline().length - 1; drawMap(); } });
["varSelect", "spSelect", "stageSelect", "depthSelect", "traitSelect", "cmapSelect"].forEach(id =>
  on(id, "change", () => { syncVarControls(); if (id === "varSelect") { const f = fieldDef(); $("logScale").checked = !!f.log; $("cmapSelect").value = f.cmap; } drawMap(); }));
["logScale", "fixedScale", "showAgents", "showCurrents", "showIbm"].forEach(id =>
  on(id, "change", () => { syncVarControls(); drawMap(); }));

// ---- species filter + aggregation controls ----
on("groupMode", "change", () => { S.groupMode = $("groupMode").value; renderFilter(); afterFilterChange(); });
const setAll = hide => { S.hidden = hide ? new Set(Array.from({ length: spCount() }, (_, q) => q)) : new Set();
                         renderFilter(); afterFilterChange(); };
on("spAll", "click", () => setAll(false));
on("spNone", "click", () => setAll(true));
on("spInvert", "click", () => {
  const inv = new Set();
  for (let q = 0; q < spCount(); q++) if (!S.hidden.has(q)) inv.add(q);
  S.hidden = inv; renderFilter(); afterFilterChange();
});
on("timelineSelect", "change", () => {                    // keep the same model day when switching streams
  const day = S.lastDay;
  const tl = timeline();
  S.idx = day === undefined ? Math.max(0, tl.length - 1) : nearestIdx(tl, day);
  $("timeSlider").max = Math.max(0, tl.length - 1);
  drawMap();
});
on("popSpecies", "change", renderCharts);
on("traitChartSelect", "change", renderCharts);
on("themeBtn", "click", () => {
  const root = document.documentElement;
  root.setAttribute("data-theme", root.getAttribute("data-theme") === "dark" ? "light" : "dark");
  localStorage.setItem("fishnet-theme", root.getAttribute("data-theme"));
  drawMap(); renderCharts(); if (S.tab === "val") Validation.draw();
});
on("figBtn", "click", async () => {
  $("figStatus").textContent = "running plot.py…";
  await post("figures", { run: S.run, no_anim: !$("figAnim").checked, feisty: $("figFeisty").checked });
  setTimeout(() => { $("figStatus").textContent = "figures appear as plot.py finishes each one"; }, 1500);
});
$("tabs").addEventListener("click", e => {
  const b = e.target.closest(".tab"); if (!b) return;
  S.tab = b.dataset.tab;
  document.querySelectorAll(".tab").forEach(t => t.classList.toggle("active", t === b));
  document.querySelectorAll(".view").forEach(v => v.classList.toggle("active", v.id === "view-" + S.tab));
  if (S.tab === "figs") renderFigures(); else loadSeries(true);
  if (S.tab === "map") drawMap();
  if (S.tab === "val") Validation.show(S.run);
});
window.addEventListener("resize", () => { clearTimeout(S.rt); S.rt = setTimeout(() => {
  drawMap(); renderCharts(); if (S.tab === "val") Validation.draw(); }, 150); });
document.addEventListener("keydown", e => {
  if (e.target.matches("input, textarea, select")) return;
  if (e.code === "Space") { e.preventDefault(); setPlaying(!S.playing); }
  if (e.key === "ArrowRight" && S.times) { S.idx = Math.min(S.idx + 1, timeline().length - 1); S.follow = false; $("followLive").checked = false; drawMap(); }
  if (e.key === "ArrowLeft" && S.times) { S.idx = Math.max(S.idx - 1, 0); S.follow = false; $("followLive").checked = false; drawMap(); }
});

window.addEventListener("unhandledrejection", e => console.error("unhandled:", e.reason && (e.reason.stack || e.reason.message || e.reason)));
window.FISHNET = S;                                 // handy for debugging in the browser console

const saved = localStorage.getItem("fishnet-theme");
if (saved) document.documentElement.setAttribute("data-theme", saved);
MapView.init($("mapCanvas"), $("mapTip"), $("mapLegend"));
Validation.wire();
setPlaying(false);
refreshState().then(() => setInterval(poll, 1800));
})();
