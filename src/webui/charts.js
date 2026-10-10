/* Minimal SVG chart layer: multi-series lines, stacked areas, bands.
   Every chart gets a legend, selective direct labels, a crosshair and a tooltip. */
const Charts = (() => {
  const NS = "http://www.w3.org/2000/svg";
  const el = (n, a = {}) => { const e = document.createElementNS(NS, n); for (const k in a) e.setAttribute(k, a[k]); return e; };

  let tip = null;
  function tooltip() {
    if (!tip) { tip = document.createElement("div"); tip.className = "chart-tip"; tip.style.display = "none"; document.body.appendChild(tip); }
    return tip;
  }

  const fmt = (v, digits = 3) => {
    if (v === null || v === undefined || !isFinite(v)) return "—";
    const a = Math.abs(v);
    if (a === 0) return "0";
    if (a >= 1e5 || a < 1e-3) return v.toExponential(2).replace("e+", "e");
    return Number(v.toPrecision(digits)).toLocaleString();
  };

  function niceTicks(lo, hi, n = 5) {
    if (!(hi > lo)) return [lo];
    const raw = (hi - lo) / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => s >= raw) || 10 * mag;
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-6; v += step) out.push(+v.toFixed(12));
    return out;
  }
  function logTicks(lo, hi) {
    const out = [];
    for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) {
      const v = Math.pow(10, e);
      if (v >= lo * 0.999 && v <= hi * 1.001) out.push(v);
    }
    return out.length > 1 ? out : niceTicks(lo, hi, 4);
  }
  const logFmt = v => {
    const e = Math.log10(v);
    return Math.abs(e - Math.round(e)) < 1e-9 ? (Math.abs(e) <= 3 ? fmt(v) : "1e" + Math.round(e)) : fmt(v);
  };

  function empty(node, msg) {
    node.innerHTML = `<div class="empty">${msg}</div>`;
  }

  /* opts: {x, series:[{name,color,y,band?:[lo,hi]}], log, stacked, now, unit, yTicks, zeroLine} */
  function draw(node, opts) {
    const x = opts.x || [], series = (opts.series || []).filter(s => s.y && s.y.length);
    if (!x.length || !series.length) return empty(node, opts.emptyMsg || "no data yet");

    const W = Math.max(node.clientWidth, 220), H = Math.max(node.clientHeight, 110);
    const showLegend = series.length >= 2 && opts.legend !== false;
    // Lay the legend out first and measure it, so the plot shrinks to make room instead of the legend
    // spilling under the next card (a dozen functional groups wrap onto several rows).
    node.innerHTML = "";
    let lg = null;
    if (showLegend) {
      lg = document.createElement("div");
      lg.className = "chart-legend";
      lg.innerHTML = series.map(s => `<span class="row">${swatch(s, opts.stacked)}${s.name}</span>`).join("");
      node.appendChild(lg);
    }
    if (lg && lg.offsetHeight > 0.45 * H) {             // too many entries for this card: say so, keep the plot
      lg.innerHTML = `<span class="row muted">${series.length} series · hover for names (the species panel and map legend list them)</span>`;
    }
    const svgH = Math.max(80, H - (lg ? lg.offsetHeight + 2 : 0));
    // direct labels only pay for themselves with 2–4 lines; one series is named by the title
    const labelled = !opts.stacked && series.length >= 2 && series.length <= 4;
    const bands = series.length <= 4;                   // more than a few ±SD bands only blur into one another
    const labelW = labelled ? Math.min(10 + Math.max(...series.map(s => s.name.length)) * 6.1, 0.3 * W) : 0;
    const m = { l: 48, r: labelled ? labelW : 12, t: 8, b: 20 };
    const iw = Math.max(10, W - m.l - m.r), ih = Math.max(10, svgH - m.t - m.b);

    // ---- scales
    const x0 = x[0], x1 = x[x.length - 1] || x0 + 1;
    const sx = v => m.l + (x1 === x0 ? 0 : (v - x0) / (x1 - x0)) * iw;

    let stackTop = null;
    if (opts.stacked) {
      stackTop = new Array(x.length).fill(0);
      for (const s of series) for (let i = 0; i < x.length; i++) stackTop[i] += Math.max(s.y[i] || 0, 0);
    }
    let lo = Infinity, hi = -Infinity;
    const consider = v => { if (isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; } };
    if (opts.stacked) { stackTop.forEach(consider); lo = 0; }
    else for (const s of series) {
      for (let i = 0; i < x.length; i++) {
        consider(s.y[i]);
        if (s.band && series.length <= 4) { consider(s.band[0][i]); consider(s.band[1][i]); }
      }
    }
    if (!isFinite(lo) || !isFinite(hi)) return empty(node, "no finite values");

    const log = !!opts.log;
    if (log) {
      let pos = Infinity;
      for (const s of series) for (const v of s.y) if (v > 0 && v < pos) pos = v;
      if (!isFinite(pos)) return empty(node, "no positive values to plot on a log scale");
      lo = Math.max(pos, hi / 1e8); hi = Math.max(hi, lo * 10);
      lo = Math.pow(10, Math.floor(Math.log10(lo))); hi = Math.pow(10, Math.ceil(Math.log10(hi)));
    } else {
      if (hi === lo) { hi = lo + (Math.abs(lo) || 1) * 0.1; lo -= (Math.abs(lo) || 1) * 0.1; }
      const pad = (hi - lo) * 0.08; hi += pad; lo -= (lo >= 0 && lo - pad < 0 ? lo : pad);
    }
    const clampLog = v => Math.max(v, lo);
    const sy = v => {
      const t = log ? (Math.log10(clampLog(v)) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo))
                    : (v - lo) / (hi - lo);
      return m.t + ih - Math.min(Math.max(t, 0), 1) * ih;
    };

    // ---- svg
    const svg = el("svg", { viewBox: `0 0 ${W} ${svgH}`, preserveAspectRatio: "none" });
    svg.setAttribute("width", W); svg.setAttribute("height", svgH);
    const ticks = [];                                      // thin labels that would collide (many-decade log axes)
    for (const t of (log ? logTicks(lo, hi) : niceTicks(lo, hi, ih > 150 ? 5 : 4))) {
      if (!ticks.length || Math.abs(sy(t) - sy(ticks[ticks.length - 1])) >= 20) ticks.push(t);
    }
    for (const t of ticks) {
      const y = sy(t);
      svg.appendChild(el("line", { class: "gridline", x1: m.l, x2: m.l + iw, y1: y, y2: y }));
      const lbl = el("text", { class: "tick", x: m.l - 6, y: y + 3, "text-anchor": "end" });
      lbl.textContent = log ? logFmt(t) : fmt(t);
      svg.appendChild(lbl);
    }
    for (const t of niceTicks(x0, x1, Math.max(2, Math.floor(iw / 90)))) {
      const px = sx(t);
      const lbl = el("text", { class: "tick", x: px, y: svgH - 6, "text-anchor": "middle" });
      lbl.textContent = opts.xFmt ? opts.xFmt(t) : fmt(t, 4);
      svg.appendChild(lbl);
    }
    svg.appendChild(el("line", { class: "axis-line", x1: m.l, x2: m.l + iw, y1: m.t + ih, y2: m.t + ih }));
    if (opts.zeroLine && lo < 0 && hi > 0) svg.appendChild(el("line", { class: "axis-line", x1: m.l, x2: m.l + iw, y1: sy(0), y2: sy(0) }));

    // ---- marks
    const path = pts => pts.length ? "M" + pts.map(p => `${p[0].toFixed(1)},${p[1].toFixed(1)}`).join("L") : "";
    if (opts.stacked) {
      const base = new Array(x.length).fill(0);
      series.forEach(s => {
        const top = base.map((b, i) => b + Math.max(s.y[i] || 0, 0));
        const up = x.map((xv, i) => [sx(xv), sy(top[i])]);
        const dn = x.map((xv, i) => [sx(xv), sy(base[i])]).reverse();
        svg.appendChild(el("path", { d: path(up) + "L" + path(dn).slice(1) + "Z", fill: s.color, "fill-opacity": .9,
                                     stroke: "var(--surface-1)", "stroke-width": 2 }));
        for (let i = 0; i < x.length; i++) base[i] = top[i];
      });
    } else {
      series.forEach(s => {
        if (s.band && bands) {
          const up = x.map((xv, i) => [sx(xv), sy(s.band[1][i])]);
          const dn = x.map((xv, i) => [sx(xv), sy(s.band[0][i])]).reverse();
          svg.appendChild(el("path", { d: path(up) + "L" + path(dn).slice(1) + "Z", fill: s.color, "fill-opacity": .13, stroke: "none" }));
        }
        const pts = [];
        for (let i = 0; i < x.length; i++) if (isFinite(s.y[i]) && (!log || s.y[i] > 0)) pts.push([sx(x[i]), sy(s.y[i])]);
        const line = el("path", { class: "series-line", d: path(pts), stroke: s.color });
        if (s.dash) { line.setAttribute("stroke-dasharray", s.dash); line.style.strokeLinecap = "butt"; }
        svg.appendChild(line);
      });
    }

    if (opts.now !== undefined && opts.now !== null && opts.now >= x0 && opts.now <= x1) {
      svg.appendChild(el("line", { class: "nowline", x1: sx(opts.now), x2: sx(opts.now), y1: m.t, y2: m.t + ih }));
    }

    if (labelled) {
      const used = [];
      series.forEach(s => {
        let i = x.length - 1;
        while (i >= 0 && !(isFinite(s.y[i]) && (!log || s.y[i] > 0))) i--;
        if (i < 0) return;
        let y = sy(s.y[i]);
        while (used.some(u => Math.abs(u - y) < 11)) y += 11;
        used.push(y);
        const t = el("text", { class: "dlabel", x: m.l + iw + 5, y: y + 3, fill: s.color });
        t.textContent = s.name;
        svg.appendChild(t);
      });
    }

    // ---- hover
    const cross = el("line", { class: "crosshair", y1: m.t, y2: m.t + ih, x1: 0, x2: 0, opacity: 0 });
    svg.appendChild(cross);
    const dots = series.map(s => {
      const c = el("circle", { r: 3.5, fill: s.color, stroke: "var(--surface-1)", "stroke-width": 2, opacity: 0 });
      svg.appendChild(c); return c;
    });
    const hit = el("rect", { x: m.l, y: m.t, width: iw, height: ih, fill: "transparent" });
    svg.appendChild(hit);
    const tp = tooltip();
    hit.addEventListener("mousemove", ev => {
      const r = svg.getBoundingClientRect();
      const px = (ev.clientX - r.left) * (W / r.width);
      const frac = (px - m.l) / iw;
      let i = Math.round(frac * (x.length - 1));
      i = Math.min(Math.max(i, 0), x.length - 1);
      cross.setAttribute("x1", sx(x[i])); cross.setAttribute("x2", sx(x[i])); cross.setAttribute("opacity", 1);
      let rows = "";
      const base = opts.stacked ? new Array(series.length) : null;
      if (opts.stacked) { let acc = 0; series.forEach((s, k) => { acc += Math.max(s.y[i] || 0, 0); base[k] = acc; }); }
      series.forEach((s, k) => {
        const v = s.y[i];
        const yv = opts.stacked ? base[k] : v;
        if (isFinite(yv) && (!log || yv > 0)) {
          dots[k].setAttribute("cx", sx(x[i])); dots[k].setAttribute("cy", sy(yv)); dots[k].setAttribute("opacity", 1);
        } else dots[k].setAttribute("opacity", 0);
        rows += `<div class="r"><span>${swatch(s, opts.stacked, "sw")}${s.name}</span><span>${fmt(v)}</span></div>`;
      });
      tp.innerHTML = `<div class="t">${opts.xLabel || "day"} ${opts.xFmt ? opts.xFmt(x[i]) : fmt(x[i], 5)}${opts.unit ? " · " + opts.unit : ""}</div>${rows}`;
      tp.style.display = "block";
      const tw = tp.offsetWidth, th = tp.offsetHeight;
      tp.style.left = Math.min(ev.clientX + 14, window.innerWidth - tw - 8) + "px";
      tp.style.top = Math.max(8, Math.min(ev.clientY - th / 2, window.innerHeight - th - 8)) + "px";
    });
    hit.addEventListener("mouseleave", () => {
      cross.setAttribute("opacity", 0); dots.forEach(d => d.setAttribute("opacity", 0)); tp.style.display = "none";
    });

    svg.style.height = svgH + "px";
    node.insertBefore(svg, lg);
  }

  /* Legend and tooltip key: a square for areas and solid lines, a short line sample for dashed ones. */
  function swatch(s, area, cls = "swatch") {
    if (s.dash && !area)
      return `<svg class="glyph" width="16" height="8" aria-hidden="true"><line x1="1" y1="4" x2="15" y2="4" stroke="${s.color}" stroke-width="2.4" stroke-dasharray="${s.dash}"/></svg>`;
    return `<i class="${cls}" style="background:${s.color}"></i>`;
  }

  /* opts: {points:[{x,y,name?,color,r?}], logX, logY, diag (1:1 line), factor (dashed lines at ×/÷ factor, log axes),
             xLabel, yLabel, onClick(point), emptyMsg}. Non-finite points (and ≤0 on a log axis) are dropped. */
  function scatter(node, opts) {
    const logX = !!opts.logX, logY = !!opts.logY;
    const pts = (opts.points || []).filter(p => isFinite(p.x) && isFinite(p.y) && p.x !== null && p.y !== null &&
                                                (!logX || p.x > 0) && (!logY || p.y > 0));
    if (!pts.length) return empty(node, opts.emptyMsg || "no data");
    const W = Math.max(node.clientWidth, 220), H = Math.max(node.clientHeight, 140);
    const m = { l: 52, r: 12, t: 10, b: 32 };
    const iw = W - m.l - m.r, ih = H - m.t - m.b;
    const tr = (v, lg) => lg ? Math.log10(v) : v;
    const span = (key, lg) => {
      let lo = Math.min(...pts.map(p => tr(p[key], lg))), hi = Math.max(...pts.map(p => tr(p[key], lg)));
      if (hi === lo) { lo -= lg ? 0.5 : Math.abs(lo) * 0.1 || 1; hi += lg ? 0.5 : Math.abs(hi) * 0.1 || 1; }
      const pad = (hi - lo) * 0.06;
      return lg ? [Math.floor(lo - pad), Math.ceil(hi + pad)] : [lo - pad, hi + pad];
    };
    let [x0, x1] = span("x", logX), [y0, y1] = span("y", logY);
    if (opts.diag && logX === logY) { x0 = y0 = Math.min(x0, y0); x1 = y1 = Math.max(x1, y1); }   // square 1:1 frame
    const sx = v => m.l + (tr(v, logX) - x0) / (x1 - x0) * iw;
    const sy = v => m.t + ih - (tr(v, logY) - y0) / (y1 - y0) * ih;
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}` });
    svg.setAttribute("width", W); svg.setAttribute("height", H);
    const ticks = (lo, hi, lg, n) => lg ? Array.from({ length: hi - lo + 1 }, (_, k) => Math.pow(10, lo + k))
                                        .filter((_, k, a) => a.length <= 7 || k % Math.ceil(a.length / 6) === 0)
                                    : niceTicks(lo, hi, n);
    for (const t of ticks(y0, y1, logY, 4)) {
      const y = sy(t);
      if (y < m.t - 1 || y > m.t + ih + 1) continue;
      svg.appendChild(el("line", { class: "gridline", x1: m.l, x2: m.l + iw, y1: y, y2: y }));
      const l = el("text", { class: "tick", x: m.l - 6, y: y + 3, "text-anchor": "end" });
      l.textContent = logY ? logFmt(t) : fmt(t); svg.appendChild(l);
    }
    for (const t of ticks(x0, x1, logX, Math.max(2, Math.floor(iw / 80)))) {
      const x = sx(t);
      if (x < m.l - 1 || x > m.l + iw + 1) continue;
      svg.appendChild(el("line", { class: "gridline", x1: x, x2: x, y1: m.t, y2: m.t + ih }));
      const l = el("text", { class: "tick", x, y: m.t + ih + 13, "text-anchor": "middle" });
      l.textContent = logX ? logFmt(t) : fmt(t); svg.appendChild(l);
    }
    svg.appendChild(el("line", { class: "axis-line", x1: m.l, x2: m.l + iw, y1: m.t + ih, y2: m.t + ih }));
    const xl = el("text", { class: "tick", x: m.l + iw, y: H - 3, "text-anchor": "end" }); xl.textContent = opts.xLabel || "";
    const yl = el("text", { class: "tick", x: 4, y: m.t + 2, transform: `rotate(-90 4 ${m.t + 2})`, "text-anchor": "end" });
    yl.textContent = opts.yLabel || "";
    svg.appendChild(xl); svg.appendChild(yl);
    if (opts.diag) {                                           // reference lines, clipped to the frame
      const clip = el("clipPath", { id: "c" + Math.random().toString(36).slice(2) });
      clip.appendChild(el("rect", { x: m.l, y: m.t, width: iw, height: ih }));
      svg.appendChild(clip);
      const g = el("g", { "clip-path": `url(#${clip.id})` });
      const lo = Math.min(x0, y0), hi = Math.max(x1, y1);
      const line = (k, dash) => {                              // y = x (+k decades on log axes)
        const a = logX ? Math.pow(10, lo) : lo, b = logX ? Math.pow(10, hi) : hi, f = logY ? Math.pow(10, k) : 1;
        g.appendChild(el("line", { x1: sx(a), y1: sy(a * f), x2: sx(b), y2: sy(b * f), class: dash ? "ref-dash" : "ref-line" }));
      };
      line(0, false);
      if (opts.factor && logX && logY) { const k = Math.log10(opts.factor); line(k, true); line(-k, true); }
      svg.appendChild(g);
    }
    const tp = tooltip();
    const small = pts.length > 60;
    for (const p of pts) {
      const c = el("circle", { cx: sx(p.x), cy: sy(p.y), r: p.r || (small ? 2.6 : 4.5), fill: p.color || "var(--series-1)",
                               "fill-opacity": small ? 0.55 : 0.9, stroke: small ? "none" : "var(--surface-1)", "stroke-width": 1.2 });
      if (p.name || opts.onClick) {
        c.style.cursor = opts.onClick ? "pointer" : "default";
        c.addEventListener("mousemove", ev => {
          tp.innerHTML = `<div class="t">${p.name || ""}</div>` +
            `<div class="r"><span>${opts.xLabel || "x"}</span><span>${fmt(p.x)}</span></div>` +
            `<div class="r"><span>${opts.yLabel || "y"}</span><span>${fmt(p.y)}</span></div>` + (p.extra || "");
          tp.style.display = "block";
          tp.style.left = Math.min(ev.clientX + 14, window.innerWidth - tp.offsetWidth - 8) + "px";
          tp.style.top = Math.max(8, ev.clientY - tp.offsetHeight / 2) + "px";
        });
        c.addEventListener("mouseleave", () => { tp.style.display = "none"; });
        if (opts.onClick) c.addEventListener("click", () => opts.onClick(p));
      }
      svg.appendChild(c);
    }
    if (opts.labels) for (const p of pts) {
      const t = el("text", { class: "dlabel", x: sx(p.x) + 6, y: sy(p.y) + 3 }); t.textContent = p.name; svg.appendChild(t);
    }
    node.innerHTML = "";
    node.appendChild(svg);
  }

  return { draw, scatter, fmt, empty };
})();
