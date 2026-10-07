/* Validation tab: the run against observations/ (obsval.py does the matching on the server).
   Cross-species scatters, a scorecard, and model / observed / agreement maps for one species. */
const Validation = (() => {
  const $ = id => document.getElementById(id);
  const css = v => getComputedStyle(document.body).getPropertyValue(v).trim();
  const fmt = Charts.fmt;
  const nice = s => s.replace(/_/g, " ");
  let D = null, run = null, sel = null, busy = false, mapData = null, loadedKey = "";
  let loadedN = -1, lastFetch = 0, live = false;
  const LIVE_EVERY = 20000;                 // ms between recomputes while a run writes output

  const params = () => ({ run, window: $("valWindow").value });

  async function show(r, force) {
    run = r;
    if (!run) return;
    const key = JSON.stringify(params());
    if (!force && D && key === loadedKey) return draw();
    if (busy) return;
    busy = true;
    const times = window.FISHNET && window.FISHNET.times;
    const n = times ? times.field_times.length : -1;
    lastFetch = Date.now();
    live = !!(times && times.running);
    if (!D) $("valStatus").textContent = "matching model output to observations…";
    try {
      const res = await fetch("/api/validation?" + new URLSearchParams(params()));
      const d = await res.json();
      if (d.error) throw new Error(d.error);
      D = d; loadedKey = key; loadedN = n; mapData = null;
      const names = D.species.map(s => s.name);
      if (sel === null || !names.includes(sel)) sel = (D.species.find(s => s.obis) || D.species[0] || {}).name || null;
      $("valSpecies").innerHTML = D.species.map(s => `<option value="${s.name}">${nice(s.name)}</option>`).join("");
      if (sel) $("valSpecies").value = sel;
      const short = D.window > 0 && D.t1 - D.t0 + 1e-6 < D.window - 10 && D.t0 === 0;
      $("valStatus").innerHTML = (live ? `<span class="live-dot">● live</span> updates as the run writes output · ` : "") +
        `model days ${D.t0.toFixed(0)}–${D.t1.toFixed(0)} (${D.frames} snapshots) · observations ${D.years[0]}–${D.years[1]}` +
        (short ? ` · <span class="warn">run is shorter than the averaging window, so this includes the spin-up</span>` : "");
      draw();
    } catch (e) {
      if (D && loadedKey === key) {             // a live refresh hit a half-written file: keep what is shown, retry
        loadedN = -1;
        return;
      }
      D = null;
      const msg = /unknown endpoint/.test(e.message || "")
        ? "the dashboard server was started before the Validation tab existed: restart dashboard.py"
        : String(e.message || e);
      $("valStatus").textContent = msg;
      ["valTable", "valRam", "valCatch", "valLat", "valCells", "valProfile", "valIccat", "valStocks", "valYear"].forEach(id => Charts.empty($(id), msg));
    } finally { busy = false; }
  }

  const row = () => D && D.species.find(s => s.name === sel);
  const select = name => { sel = name; $("valSpecies").value = name; mapData = null; draw(); };

  function draw() {
    if (!D) return;
    if (!D.species.length) {
      $("valNote").textContent = "None of this run's species have observations in observations/ (matched by name: " +
        "anchoveta, herring, cod_atlantic, …). Nothing to compare.";
      ["valTable", "valRam", "valCatch", "valLat", "valCells", "valProfile", "valStocks", "valYear"].forEach(id => Charts.empty($(id), "no matching observations"));
      return;
    }
    note(); table(); crossScatters(); speciesPanels();
  }

  // ---------------------------------------------------------------- explanation
  function note() {
    const w = D.window > 0 ? `the last ${D.window >= 730 ? (D.window / 365).toFixed(0) + " years" : D.window + " days"}` : "the whole run";
    $("valNote").innerHTML =
      `Model: time mean of ${w} of gridded output (juveniles + adults unless stated). ` +
      `<b>OBIS</b> records are binned to model cells and compared as presence (records reflect survey effort, not abundance). ` +
      `<b>RAM Legacy</b> stock assessments are averaged over ${D.years[0]}–${D.years[1]}, summed over each species' current stocks, ` +
      `and compared with model biomass summed over the FAO areas those stocks are in. Assessments cover only part of those areas, ` +
      `so an unbiased model sits somewhat <i>above</i> the 1:1 line. <b>ICCAT</b> tuna catch on 5° squares is compared by rank with model biomass in the same squares.` +
      (D.fishing ? "" : ` <span class="warn">This run has no fishing, so the model should overshoot assessed (fished) biomass.</span>`) +
      (D.global_grid ? "" : ` <span class="warn">Regional domain: totals cover only the part of each area inside the grid.</span>`);
  }

  // ---------------------------------------------------------------- scorecard
  const grade = (v, good, fair, higherBetter = true) => v === null || v === undefined ? "" :
    (higherBetter ? (v >= good ? "g-good" : v >= fair ? "g-fair" : "g-poor") : (v <= good ? "g-good" : v <= fair ? "g-fair" : "g-poor"));
  const ratioOf = s => s && s.obs_t && s.model_t !== null ? s.model_t / s.obs_t : null;

  function table() {
    const f = (v, d = 2) => v === null || v === undefined ? "—" : v.toFixed(d);
    const rows = D.species.map(s => {
      const o = s.obis || {}, r = s.ram || {}, c = s.iccat || {}, ratio = ratioOf(r);
      const lr = ratio > 0 ? Math.abs(Math.log10(ratio)) : null;
      const dl = o.lat_model && o.lat_obs ? Math.abs(o.lat_model.mean - o.lat_obs.mean) : null;
      return `<tr data-sp="${s.name}" class="${s.name === sel ? "sel" : ""}">
        <td>${nice(s.name)}</td>
        <td class="num">${o.cells ?? "—"}</td>
        <td class="num ${grade(o.auc, 0.8, 0.7)}">${f(o.auc)}</td>
        <td class="num ${grade(o.rho, 0.4, 0.2)}">${f(o.rho)}</td>
        <td class="num ${grade(o.sens, 0.7, 0.4)}">${o.sens === undefined || o.sens === null ? "—" : (100 * o.sens).toFixed(0) + "%"}</td>
        <td class="num ${grade(dl, 5, 15, false)}">${o.lat_model ? f(o.lat_model.mean, 0) + "° / " + f(o.lat_obs && o.lat_obs.mean, 0) + "°" : "—"}</td>
        <td class="num ${ratio === 0 ? "g-poor" : grade(lr, Math.log10(3), 1, false)}">${ratio === null ? "—" : ratio === 0 ? "0" : "×" + fmt(ratio, 2)}</td>
        <td class="num ${grade(c.rho, 0.4, 0.2)}">${f(c.rho)}</td></tr>`;
    }).join("");
    $("valTable").innerHTML = `<table class="valtable"><thead><tr>
      <th>species</th><th title="model cells with at least one OBIS record">OBIS cells</th>
      <th title="chance that a cell with records has more model biomass than one without (0.5 = no skill, 1 = perfect)">AUC</th>
      <th title="Spearman rank correlation of model biomass with record count, over cells that have both">ρ cells</th>
      <th title="share of OBIS cells inside the model's core range (cells holding 95% of its biomass)">in range</th>
      <th title="mean latitude of the occupied range, model / OBIS">range centre</th>
      <th title="model biomass in the stocks' FAO areas ÷ assessed biomass">biomass model÷obs</th>
      <th title="Spearman rank correlation of model biomass with ICCAT catch over 5° squares">ICCAT ρ</th>
      </tr></thead><tbody>${rows}</tbody></table>`;
    $("valTable").querySelectorAll("tr[data-sp]").forEach(tr => tr.onclick = () => select(tr.dataset.sp));
  }

  // ---------------------------------------------------------------- cross-species scatters
  function stats(pts, log) {
    const p = pts.filter(q => isFinite(q.x) && isFinite(q.y) && (!log || (q.x > 0 && q.y > 0)));
    if (p.length < 3) return { n: p.length };
    const xs = p.map(q => log ? Math.log10(q.x) : q.x), ys = p.map(q => log ? Math.log10(q.y) : q.y);
    const mx = xs.reduce((a, b) => a + b) / xs.length, my = ys.reduce((a, b) => a + b) / ys.length;
    let sxy = 0, sxx = 0, syy = 0;
    xs.forEach((x, i) => { sxy += (x - mx) * (ys[i] - my); sxx += (x - mx) ** 2; syy += (ys[i] - my) ** 2; });
    const d = xs.map((x, i) => ys[i] - x).sort((a, b) => a - b);
    return { n: p.length, r: sxy / Math.sqrt(sxx * syy), med: d[Math.floor(d.length / 2)], d };
  }
  const pointColor = name => name === sel ? css("--series-2") : css("--series-1");

  function crossScatters() {
    const onClick = p => select(p.name);
    // biomass vs assessments
    const ram = D.species.filter(s => s.ram && s.ram.obs_t);
    const pts = ram.map(s => ({ x: s.ram.obs_t, y: s.ram.model_t, name: s.name, color: pointColor(s.name), r: s.name === sel ? 6 : undefined,
                               extra: `<div class="r"><span>${s.ram.metric}, ${s.ram.n_used}/${s.ram.n_stocks} stocks</span></div>` }));
    const zero = ram.filter(s => !(s.ram.model_t > 0)).map(s => nice(s.name));
    const st = stats(pts, true);
    $("valRamS").textContent = st.n >= 3 ? `r(log) ${st.r.toFixed(2)} · median model÷obs ×${fmt(Math.pow(10, st.med), 2)} · ` +
      `${st.d.filter(v => Math.abs(v) <= 1).length}/${st.n} within ×10` : "";
    Charts.scatter($("valRam"), { points: pts, logX: true, logY: true, diag: true, factor: 10, onClick,
      xLabel: "assessed (t)", yLabel: "model (t)", emptyMsg: "no species with RAM Legacy assessments" });
    $("valRamZero").textContent = zero.length ? `not plotted, model has none in the assessed areas: ${zero.join(", ")}` : "";

    // catch
    const cs = D.species.filter(s => s.ram && s.ram.catch_obs);
    const cp = cs.map(s => ({ x: s.ram.catch_obs, y: s.ram.catch_model, name: s.name, color: pointColor(s.name) }));
    const cst = stats(cp, true);
    $("valCatchS").textContent = D.fishing && cst.n >= 3 ? `r(log) ${cst.r.toFixed(2)} · median model÷obs ×${fmt(Math.pow(10, cst.med), 2)}` : "";
    Charts.scatter($("valCatch"), { points: D.fishing ? cp : [], logX: true, logY: true, diag: true, factor: 10, onClick,
      xLabel: "reported (t yr⁻¹)", yLabel: "model (t yr⁻¹)",
      emptyMsg: D.fishing ? "no species with reported catch" : "no fishing in this run (fishing_F = 0): nothing to compare with reported catch" });

    // range centre latitude
    const ls = D.species.filter(s => s.obis && s.obis.lat_model && s.obis.lat_obs);
    const lp = ls.map(s => ({ x: s.obis.lat_obs.mean, y: s.obis.lat_model.mean, name: s.name, color: pointColor(s.name),
      extra: `<div class="r"><span>OBIS 5–95%</span><span>${s.obis.lat_obs.p5.toFixed(0)}…${s.obis.lat_obs.p95.toFixed(0)}°</span></div>` +
             `<div class="r"><span>model 5–95%</span><span>${s.obis.lat_model.p5.toFixed(0)}…${s.obis.lat_model.p95.toFixed(0)}°</span></div>` }));
    const lst = stats(lp, false);
    $("valLatS").textContent = lst.n >= 3 ? `r ${lst.r.toFixed(2)} · mean |Δ| ${(lst.d.reduce((a, b) => a + Math.abs(b), 0) / lst.n).toFixed(0)}°` : "";
    Charts.scatter($("valLat"), { points: lp, diag: true, onClick, xLabel: "OBIS (°N)", yLabel: "model (°N)", emptyMsg: "no OBIS data" });
  }

  // ---------------------------------------------------------------- one species
  async function speciesPanels() {
    const s = row();
    $("valSpName").textContent = s ? `${nice(s.name)}${s.sci ? " · " + s.sci : ""}` : "";
    if (!s) return;
    const o = s.obis;
    // per-cell scatter
    if (o && o.scatter.obs.length) {
      const pts = o.scatter.obs.map((x, i) => ({ x, y: o.scatter.model[i], color: css("--series-1") }));
      $("valCellsS").textContent = `${o.n_rho} cells · Spearman ρ ${o.rho === null ? "—" : o.rho.toFixed(2)}`;
      Charts.scatter($("valCells"), { points: pts, logX: true, logY: true, xLabel: "OBIS records", yLabel: "model g m⁻²" });
    } else { $("valCellsS").textContent = ""; Charts.empty($("valCells"), "no cells with both OBIS records and model biomass"); }
    // latitude occupancy
    if (o) {
      const nan = a => a.map(v => v === null ? NaN : v);
      Charts.draw($("valProfile"), { x: o.profile.lat, xLabel: "lat", series: [
        { name: "model (core range)", color: css("--series-1"), y: nan(o.profile.model) },
        { name: "OBIS (any record)", color: css("--series-2"), y: nan(o.profile.obs) }] });
    } else Charts.empty($("valProfile"), "no OBIS data for this species");
    // ICCAT
    const c = s.iccat;
    $("valIccatCard").classList.toggle("hidden", !c);
    if (c) {
      const pts = c.scatter.obs.map((x, i) => ({ x, y: c.scatter.model[i], color: css("--series-3") }));
      $("valIccatS").textContent = `${c.boxes} squares · Spearman ρ ${c.rho === null ? "—" : c.rho.toFixed(2)}` +
        (c.r_log === null ? "" : ` · r(log) ${c.r_log.toFixed(2)}`);
      Charts.scatter($("valIccat"), { points: pts, logX: true, logY: true, xLabel: "ICCAT catch (t yr⁻¹)", yLabel: "model g m⁻²",
        emptyMsg: "model has no biomass in the ICCAT squares" });
    }
    stocks(s);
    yearly(s);
    const src = $("valSource");
    src.querySelector('option[value="iccat"]').disabled = !c;
    if (!c && src.value === "iccat") src.value = "obis";
    await maps(s);
  }

  function yearly(s) {
    const y = s.ram && s.ram.yearly;
    if (!y || !y.years.length) { $("valYearS").textContent = ""; return Charts.empty($("valYear"), "no assessed yearly series for this species"); }
    const nan = a => a.map(v => v === null ? NaN : v);
    const what = y.metric === "tb" ? "total biomass (juveniles + adults)" : "spawning biomass (adults)";
    $("valYearS").textContent = `${what}, t, in the stocks' FAO areas · RAM = sum of the stocks reporting each year (of ${y.of})` +
      (y.r_change === null ? "" : ` · year-to-year change correlation ${y.r_change.toFixed(2)}`) +
      " · meaningful only when the run uses real (dated) forcing";
    Charts.draw($("valYear"), { x: y.years, xLabel: "year", xFmt: v => String(Math.round(v)), log: true, unit: "t", series: [
      { name: "model", color: css("--series-1"), y: nan(y.model) },
      { name: "RAM Legacy", color: css("--series-2"), y: nan(y.obs) }] });
  }

  function stocks(s) {
    const r = s.ram;
    if (!r) { Charts.empty($("valStocks"), "no RAM Legacy assessments for this species"); $("valStocksS").textContent = ""; return; }
    $("valStocksS").textContent = `${r.metric}; model over FAO areas ${r.areas.join(", ") || "—"} = ${fmt(r.model_t)} t (whole domain ${fmt(r.model_domain_t)} t)`;
    $("valStocks").innerHTML = `<table class="valtable"><thead><tr><th>stock</th><th>region</th><th>FAO</th>
      <th>total biomass (t)</th><th>SSB (t)</th><th>catch (t yr⁻¹)</th></tr></thead><tbody>` +
      r.stocks.map(k => `<tr><td title="${k.id}">${k.name}</td><td>${k.region}</td><td class="num">${k.areas.join(",")}</td>
        <td class="num">${k.tb === undefined ? "—" : fmt(k.tb)}</td><td class="num">${k.ssb === undefined ? "—" : fmt(k.ssb)}</td>
        <td class="num">${k.catch === undefined ? "—" : fmt(k.catch)}</td></tr>`).join("") +
      `<tr class="tot"><td>sum used (${r.n_used} of ${r.n_stocks})</td><td></td><td></td><td colspan="2" class="num">${fmt(r.obs_t)}</td>` +
      `<td class="num">${fmt(r.catch_obs)}</td></tr></tbody></table>`;
  }

  // ---------------------------------------------------------------- maps
  async function maps(s) {
    const meta = window.FISHNET && window.FISHNET.meta;
    if (!meta || !meta.nx) return;
    const src = $("valSource").value;
    const key = `${run}|${s.q}|${src}|${$("valWindow").value}|${D.t1}`;
    if (!mapData || mapData.key !== key) {
      const r = await fetch("/api/valmap?" + new URLSearchParams({ ...params(), species: s.q, source: src }));
      if (!r.ok) return;
      const a = new Float32Array(await r.arrayBuffer()), n = meta.nx * meta.ny;
      mapData = { key, model: a.subarray(0, n), obs: a.subarray(n, 2 * n) };
    }
    const { model, obs } = mapData, n = model.length;
    // model core range: smallest set of cells holding 95% of the biomass (same rule as the server)
    const idx = Array.from({ length: n }, (_, i) => i).filter(i => model[i] > 0).sort((a, b) => model[b] - model[a]);
    const tot = idx.reduce((a, i) => a + model[i], 0), core = new Uint8Array(n);
    let acc = 0;
    for (const i of idx) { core[i] = 1; acc += model[i]; if (acc >= 0.95 * tot) break; }
    const cat = new Float32Array(n);
    for (let i = 0; i < n; i++) cat[i] = !isFinite(model[i]) ? NaN : (core[i] ? 1 : 0) + (obs[i] > 0 ? 2 : 0);
    const obsLabel = src === "iccat" ? "ICCAT catch (t yr⁻¹ per cell)" : "OBIS records per cell";
    grid($("valMapModel"), meta, model, { cmap: "blues", log: true, label: "model g m⁻²" });
    grid($("valMapObs"), meta, obs, { cmap: "inferno", log: true, label: obsLabel });
    const colors = ["", css("--series-1"), css("--series-2"), css("--series-3")];
    grid($("valMapAgree"), meta, cat, { categories: colors });
    $("valMapModelL").innerHTML = legendBar("blues", model);
    $("valMapObsL").innerHTML = legendBar("inferno", obs);
    const cnt = [0, 0, 0, 0];
    for (let i = 0; i < n; i++) if (isFinite(cat[i])) cnt[cat[i]]++;
    $("valMapAgreeL").innerHTML = [[3, "both"], [1, "model only"], [2, "observed only"]].map(([k, t]) =>
      `<span class="row"><i class="swatch" style="background:${colors[k]}"></i>${t} ${cnt[k]}</span>`).join("");
    const tipOf = i => `model ${fmt(model[i])} g m⁻²<br>${src === "iccat" ? "ICCAT" : "OBIS"} ${fmt(obs[i])}`;
    [$("valMapModel"), $("valMapObs"), $("valMapAgree")].forEach(cv => cv._tip = tipOf);
  }

  function range(a) {
    const v = []; for (const x of a) if (x > 0 && isFinite(x)) v.push(x);
    if (!v.length) return [1e-3, 1];
    v.sort((p, q) => p - q);
    const hi = v[Math.floor(0.995 * (v.length - 1))];
    const lo = Math.max(v[Math.floor(0.02 * (v.length - 1))], hi * 1e-4);   // four decades: traces stay visible
    return [lo, hi > lo ? hi : lo * 10];
  }
  function legendBar(cmap, a) {
    const [lo, hi] = range(a);
    return `<div class="bar" style="background:linear-gradient(90deg,${MapView.stops[cmap].join(",")})"></div>` +
           `<div class="ends"><span>${fmt(lo)}</span><span>log</span><span>${fmt(hi)}</span></div>`;
  }

  /* Draw a (ny*nx, row 0 south) grid into a canvas: log colour ramp, or fixed colours per integer category. */
  function grid(cv, meta, a, opt) {
    const wrap = cv.parentElement, dpr = window.devicePixelRatio || 1;
    const w = Math.max(wrap.clientWidth, 120);
    const latm = (meta.lat[0] + meta.lat[meta.ny - 1]) / 2;
    const aspect = (meta.nx * Math.abs(meta.lon[1] - meta.lon[0]) * Math.cos(latm * Math.PI / 180)) /
                   (meta.ny * Math.abs(meta.lat[1] - meta.lat[0]));
    const h = Math.round(w / Math.max(aspect, 0.5));
    cv.style.width = w + "px"; cv.style.height = h + "px"; cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
    const off = document.createElement("canvas"); off.width = meta.nx; off.height = meta.ny;
    const oc = off.getContext("2d"), img = oc.createImageData(meta.nx, meta.ny);
    const hex = c => { c = c.startsWith("#") ? c : "#888888"; return [1, 3, 5].map(k => parseInt(c.slice(k, k + 2), 16)); };
    const land = hex(css("--land")), sea = hex(css("--surface-2"));
    const table = opt.categories ? null : MapView.lut(opt.cmap);
    const cats = opt.categories ? opt.categories.map(c => c ? hex(c) : sea) : null;
    const [lo, hi] = opt.categories ? [0, 1] : range(a);
    const l0 = Math.log10(lo), l1 = Math.log10(hi);
    for (let j = 0; j < meta.ny; j++) {
      const src = (meta.ny - 1 - j) * meta.nx;
      for (let i = 0; i < meta.nx; i++) {
        const v = a[src + i], p = (j * meta.nx + i) * 4;
        let c;
        if (!isFinite(v)) c = land;
        else if (cats) c = cats[v | 0];
        else if (!(v > 0)) c = sea;
        else { const t = Math.min(Math.max((Math.log10(v) - l0) / (l1 - l0 || 1), 0), 1), k = (t * 255 | 0) * 3; c = [table[k], table[k + 1], table[k + 2]]; }
        img.data[p] = c[0]; img.data[p + 1] = c[1]; img.data[p + 2] = c[2]; img.data[p + 3] = 255;
      }
    }
    oc.putImageData(img, 0, 0);
    const ctx = cv.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(off, 0, 0, w, h);
    if (!cv._hover) {
      cv._hover = true;
      const tp = document.querySelector(".chart-tip") || (() => { const d = document.createElement("div"); d.className = "chart-tip"; document.body.appendChild(d); return d; })();
      cv.addEventListener("mousemove", ev => {
        const m = window.FISHNET.meta, r = cv.getBoundingClientRect();
        const i = Math.floor((ev.clientX - r.left) / r.width * m.nx), j = m.ny - 1 - Math.floor((ev.clientY - r.top) / r.height * m.ny);
        if (i < 0 || j < 0 || i >= m.nx || j >= m.ny || !cv._tip) return;
        const k = j * m.nx + i, lat = m.lat[j], lon = m.lon[i];
        tp.innerHTML = `<div class="t">${Math.abs(lat).toFixed(0)}°${lat >= 0 ? "N" : "S"} ${Math.abs(lon).toFixed(0)}°${lon >= 0 ? "E" : "W"}</div>` +
                       (isFinite(mapData.model[k]) ? cv._tip(k) : "land");
        tp.style.display = "block";
        tp.style.left = Math.min(ev.clientX + 14, window.innerWidth - tp.offsetWidth - 8) + "px";
        tp.style.top = (ev.clientY + 12) + "px";
      });
      cv.addEventListener("mouseleave", () => { tp.style.display = "none"; });
    }
  }

  /* Called on every dashboard poll: recompute when new gridded output has appeared, at most every
     LIVE_EVERY ms while the run is going, and once more when it stops so the final state is shown. */
  function tick(r, times) {
    if (!r || r !== run || !times || busy) return;
    const n = times.field_times.length;
    if (n === loadedN) return;
    if (times.running && Date.now() - lastFetch < LIVE_EVERY) return;
    show(r, true);
  }

  function wire() {
    $("valSpecies").addEventListener("change", e => select(e.target.value));
    $("valWindow").addEventListener("change", () => show(run, true));
    $("valSource").addEventListener("change", () => { mapData = null; const s = row(); if (s) maps(s); });
    $("valRefresh").addEventListener("click", () => show(run, true));
  }

  return { show, draw, wire, tick, reset: () => { D = null; loadedKey = ""; loadedN = -1; mapData = null; } };
})();
