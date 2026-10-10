/* Canvas map: gridded field + agent dots + current vectors + IBM cells. */
const MapView = (() => {
  // sequential ramps (single hue where the palette provides one), sampled as anchor stops
  const CMAPS = {
    blues:   ["#f2f7fe", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"],
    viridis: ["#440154", "#472d7b", "#3b528b", "#2c728e", "#21918c", "#28ae80", "#5ec962", "#addc30", "#fde725"],
    inferno: ["#000004", "#1b0c41", "#4a0c6b", "#781c6d", "#a52c60", "#cf4446", "#ed6925", "#fb9b06", "#f7d13d", "#fcffa4"],
    thermal: ["#184f95", "#3987e5", "#9ec5f4", "#e6e6e2", "#f0b9a6", "#e34948", "#a11f1f"],  // diverging, gray midpoint
    // matplotlib's cyclic "banded" maps, all 256 of its steps (they repeat colours many times, so fewer stops
    // would blur the bands); good for picking out contours, not for reading values
    prism: [
      "#ff0000", "#ff0000", "#ff2100", "#ff5200", "#ff8200", "#ffb000", "#ffd800", "#fff700", "#e3ff00", "#b2ff00", "#81ff00", "#53fe00",
      "#2be200", "#0bbd39", "#00917d", "#0061b9", "#0030e9", "#0001ff", "#1a00ff", "#3e00ff", "#6a00fe", "#9a00d7", "#cb00a3", "#fa0063",
      "#ff001d", "#ff0000", "#ff0e00", "#ff3e00", "#ff6f00", "#ff9e00", "#ffc900", "#ffec00", "#f5ff00", "#c6ff00", "#95ff00", "#65ff00",
      "#3aef00", "#17cd1d", "#00a363", "#0074a2", "#0043d7", "#0013fe", "#0d00ff", "#2e00ff", "#5700ff", "#8600e9", "#b700b9", "#e7007e",
      "#ff003a", "#ff0000", "#ff0000", "#ff2a00", "#ff5b00", "#ff8c00", "#ffb900", "#ffdf00", "#fffc00", "#d9ff00", "#a9ff00", "#78ff00",
      "#4bf900", "#24dc00", "#06b547", "#00878a", "#0057c4", "#0026f0", "#0300ff", "#2000ff", "#4600ff", "#7300f8", "#a300ce", "#d40097",
      "#ff0056", "#ff000f", "#ff0000", "#ff1700", "#ff4800", "#ff7900", "#ffa700", "#ffd000", "#fff100", "#ecff00", "#bcff00", "#8bff00",
      "#5cff00", "#32e900", "#11c52b", "#009a70", "#006bae", "#0039e0", "#000aff", "#1300ff", "#3600ff", "#6000ff", "#8f00e1", "#c100af",
      "#f00071", "#ff002c", "#ff0000", "#ff0500", "#ff3400", "#ff6500", "#ff9500", "#ffc100", "#ffe500", "#ffff00", "#d0ff00", "#9fff00",
      "#6fff00", "#42f400", "#1dd50e", "#01ac55", "#007e96", "#004dcd", "#001df7", "#0800ff", "#2700ff", "#4e00ff", "#7c00f1", "#ad00c4",
      "#dd008b", "#ff0048", "#ff0001", "#ff0000", "#ff2100", "#ff5100", "#ff8200", "#ffb000", "#ffd700", "#fff700", "#e3ff00", "#b3ff00",
      "#82ff00", "#54fe00", "#2be300", "#0bbd39", "#00917d", "#0061b9", "#0030e8", "#0001ff", "#1900ff", "#3e00ff", "#6900fe", "#9900d8",
      "#ca00a3", "#f90064", "#ff001e", "#ff0000", "#ff0e00", "#ff3d00", "#ff6f00", "#ff9e00", "#ffc800", "#ffeb00", "#f6ff00", "#c6ff00",
      "#95ff00", "#66ff00", "#3bef00", "#17cd1c", "#00a462", "#0075a2", "#0044d7", "#0014fe", "#0d00ff", "#2e00ff", "#5700ff", "#8600e9",
      "#b700ba", "#e7007e", "#ff003a", "#ff0000", "#ff0000", "#ff2a00", "#ff5b00", "#ff8b00", "#ffb800", "#ffde00", "#fffb00", "#daff00",
      "#a9ff00", "#78ff00", "#4bfa00", "#24dc00", "#06b546", "#008889", "#0057c3", "#0027f0", "#0300ff", "#2000ff", "#4600ff", "#7200f8",
      "#a300cf", "#d40097", "#ff0056", "#ff0010", "#ff0000", "#ff1700", "#ff4700", "#ff7800", "#ffa700", "#ffd000", "#fff100", "#edff00",
      "#bdff00", "#8cff00", "#5dff00", "#33e900", "#11c62a", "#009b6f", "#006bad", "#003ae0", "#000bff", "#1300ff", "#3500ff", "#6000ff",
      "#8f00e1", "#c000af", "#f00071", "#ff002c", "#ff0000", "#ff0500", "#ff3300", "#ff6500", "#ff9400", "#ffc000", "#ffe500", "#ffff00",
      "#d0ff00", "#9fff00", "#6fff00", "#43f500", "#1dd50d", "#01ad54", "#007f95", "#004ecd", "#001df7", "#0800ff", "#2600ff", "#4e00ff",
      "#7c00f1", "#ac00c5", "#dd008b", "#ff0049", "#ff0001", "#ff0000", "#ff2000", "#ff5100", "#ff8200", "#ffaf00", "#ffd700", "#fff600",
      "#e4ff00", "#b3ff00", "#82ff00", "#54ff00"],
    gist_ncar: [
      "#000080", "#000777", "#000f6d", "#001664", "#001d5a", "#002451", "#002c48", "#00333e", "#003a35", "#00422b", "#004922", "#005019",
      "#00580f", "#005f06", "#005816", "#005127", "#004b37", "#004448", "#003d59", "#003669", "#002f7a", "#00298b", "#00229b", "#001bac",
      "#0014bc", "#000ecd", "#0007de", "#0000ee", "#000eff", "#001cff", "#002aff", "#0038ff", "#0047ff", "#0055ff", "#0063ff", "#0071ff",
      "#007fff", "#008dff", "#009bff", "#00a9ff", "#00b8ff", "#00c0ff", "#00c6ff", "#00caff", "#00ceff", "#00d3ff", "#00d7ff", "#00dcff",
      "#00e0ff", "#00e5ff", "#00e9ff", "#00edff", "#00f2ff", "#00f6f8", "#00fbf2", "#00ffeb", "#00ffe5", "#00fede", "#00fed8", "#00fdd1",
      "#00fdcb", "#00fcc4", "#00fcbd", "#00fbb7", "#00fbb0", "#00faaa", "#00faa3", "#00fa9d", "#00fa92", "#00fa88", "#00fa7d", "#00fb73",
      "#00fb68", "#00fc5e", "#00fc54", "#00fd49", "#00fd3f", "#00fd34", "#00fe2a", "#00fe1f", "#06ff15", "#0dff0b", "#13fb00", "#19f700",
      "#20f400", "#26f000", "#2dec00", "#33e800", "#39e500", "#40e100", "#46dd00", "#4cd900", "#53d600", "#59d200", "#60ce00", "#66d100",
      "#68d500", "#6ad800", "#6cdb00", "#6ede00", "#70e200", "#72e500", "#74e800", "#76eb00", "#78ef00", "#7af200", "#7cf500", "#7ef804",
      "#80fc08", "#84ff0c", "#89ff10", "#8dff14", "#92ff18", "#96ff1c", "#9bff20", "#9fff24", "#a4ff28", "#a9ff2c", "#adff30", "#b2ff34",
      "#b6ff38", "#bbff3c", "#bfff38", "#c4ff34", "#c8ff30", "#cdff2c", "#d2ff28", "#d6ff24", "#dbff20", "#dfff1c", "#e4ff18", "#e8ff14",
      "#edff10", "#f1ff0c", "#f6fd08", "#fafa04", "#fff800", "#fff500", "#fff300", "#fff000", "#ffee00", "#ffeb00", "#ffe900", "#ffe600",
      "#ffe400", "#ffe100", "#ffdf00", "#ffdc00", "#ffda00", "#ffd801", "#ffd502", "#ffd303", "#ffd004", "#ffce05", "#ffcb06", "#ffc908",
      "#ffc609", "#ffc40a", "#ffc10b", "#ffbf0c", "#ffbc0d", "#ffba0e", "#ffb20d", "#ffaa0c", "#ffa10b", "#ff990a", "#ff9109", "#ff8908",
      "#ff8107", "#ff7807", "#ff7006", "#ff6805", "#ff6004", "#ff5803", "#ff5002", "#ff4701", "#ff4300", "#ff3e00", "#ff3900", "#ff3400",
      "#ff3000", "#ff2b00", "#ff2600", "#ff2100", "#ff1d00", "#ff1800", "#ff1300", "#ff0e00", "#ff0a00", "#ff0512", "#ff0023", "#ff0035",
      "#ff0047", "#ff0058", "#ff006a", "#ff007c", "#ff008e", "#ff009f", "#ff00b1", "#ff00c3", "#ff00d5", "#ff00e6", "#ff00f8", "#f803fc",
      "#f107ff", "#ea0aff", "#e40eff", "#dd11ff", "#d615ff", "#cf18ff", "#c81cff", "#c11fff", "#ba22ff", "#b326ff", "#ac29ff", "#a62dff",
      "#9f33fe", "#a439fd", "#aa3ffb", "#b044fa", "#b64af9", "#bc50f8", "#c256f7", "#c85cf5", "#ce62f4", "#d468f3", "#da6ef2", "#e074f1",
      "#e67aef", "#eb80ee", "#ec84ef", "#ed89ef", "#ee8df0", "#ee92f1", "#ef97f1", "#f09bf2", "#f0a0f2", "#f1a5f3", "#f2a9f4", "#f3aef4",
      "#f3b3f5", "#f4b7f5", "#f5bcf6", "#f5c0f7", "#f6c5f7", "#f7caf8", "#f8cef9", "#f8d3f9", "#f9d8fa", "#fadcfa", "#fae1fb", "#fbe5fc",
      "#fceafc", "#fdeffd", "#fdf3fd", "#fef8fe"],
    flag: [
      "#ff0000", "#ff6035", "#ffb37e", "#ffeac6", "#ffffff", "#cdeeff", "#85b9ff", "#3c69ff", "#0009ff", "#0000d0", "#000088", "#00003f",
      "#000000", "#2c0000", "#730000", "#bc0000", "#fc0000", "#ff4f29", "#ffa570", "#ffe2b9", "#fffefa", "#d9f4ff", "#93c6ff", "#497aff",
      "#081cff", "#0000dd", "#000096", "#00004d", "#00000b", "#1f0000", "#650000", "#af0000", "#f10000", "#ff3d1c", "#ff9662", "#ffd9ab",
      "#fffbee", "#e6f9ff", "#a1d1ff", "#578aff", "#132fff", "#0000e9", "#0000a4", "#00005b", "#000016", "#130000", "#570000", "#a10000",
      "#e60000", "#ff2a11", "#ff8654", "#ffce9d", "#fff8e3", "#f1fcff", "#afdbff", "#659aff", "#1f41ff", "#0000f4", "#0000b2", "#000069",
      "#000022", "#080000", "#490000", "#930000", "#d90000", "#ff1805", "#ff7646", "#ffc38f", "#fff3d6", "#fcfeff", "#bce4ff", "#73a8ff",
      "#2c53ff", "#0000ff", "#0000c0", "#000077", "#00002f", "#000000", "#3c0000", "#850000", "#cd0000", "#ff0500", "#ff6539", "#ffb681",
      "#ffecca", "#ffffff", "#caecff", "#81b6ff", "#3965ff", "#0005ff", "#0000cd", "#000085", "#00003c", "#000000", "#2f0000", "#770000",
      "#c00000", "#ff0000", "#ff532c", "#ffa873", "#ffe4bc", "#fffefc", "#d6f3ff", "#8fc3ff", "#4676ff", "#0518ff", "#0000d9", "#000093",
      "#000049", "#000008", "#220000", "#690000", "#b20000", "#f40000", "#ff411f", "#ff9a65", "#ffdbaf", "#fffcf1", "#e3f8ff", "#9dceff",
      "#5486ff", "#112aff", "#0000e6", "#0000a1", "#000057", "#000013", "#160000", "#5b0000", "#a40000", "#e90000", "#ff2f13", "#ff8a57",
      "#ffd1a1", "#fff9e6", "#eefbff", "#abd9ff", "#6296ff", "#1c3dff", "#0000f1", "#0000af", "#000065", "#00001f", "#0b0000", "#4d0000",
      "#960000", "#dd0000", "#ff1c08", "#ff7a49", "#ffc693", "#fff4d9", "#fafeff", "#b9e2ff", "#70a5ff", "#294fff", "#0000fc", "#0000bc",
      "#000073", "#00002c", "#000000", "#3f0000", "#880000", "#d00000", "#ff0900", "#ff693c", "#ffb985", "#ffeecd", "#ffffff", "#c6eaff",
      "#7eb3ff", "#3560ff", "#0000ff", "#0000ca", "#000081", "#000039", "#000000", "#320000", "#7a0000", "#c30000", "#ff0000", "#ff582f",
      "#ffac77", "#ffe6c0", "#ffffff", "#d3f1ff", "#8cc0ff", "#4372ff", "#0313ff", "#0000d6", "#00008f", "#000046", "#000005", "#260000",
      "#6c0000", "#b60000", "#f70000", "#ff4622", "#ff9d69", "#ffdeb2", "#fffdf4", "#e0f7ff", "#9acbff", "#5082ff", "#0e26ff", "#0000e3",
      "#00009d", "#000054", "#000011", "#190000", "#5e0000", "#a80000", "#ec0000", "#ff3316", "#ff8e5b", "#ffd4a4", "#fffae9", "#ecfbff",
      "#a8d6ff", "#5e92ff", "#1938ff", "#0000ee", "#0000ab", "#000062", "#00001c", "#0e0000", "#500000", "#9a0000", "#e00000", "#ff210b",
      "#ff7e4d", "#ffc996", "#fff5dd", "#f7fdff", "#b6e0ff", "#6ca1ff", "#264aff", "#0000fa", "#0000b9", "#000070", "#000029", "#030000",
      "#430000", "#8c0000", "#d30000", "#ff0e00", "#ff6d3f", "#ffbc88", "#ffefd0", "#ffffff", "#c3e8ff", "#7aafff", "#325cff", "#0000ff",
      "#0000c6", "#00007e", "#000035", "#000000"]
  };
  const LUT = {};
  const hex2rgb = h => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
  function lut(name) {
    if (LUT[name]) return LUT[name];
    const stops = (CMAPS[name] || CMAPS.blues).map(hex2rgb), n = 256, out = new Uint8ClampedArray(n * 3);
    for (let i = 0; i < n; i++) {
      const t = i / (n - 1) * (stops.length - 1), k = Math.min(Math.floor(t), stops.length - 2), f = t - k;
      for (let c = 0; c < 3; c++) out[i * 3 + c] = stops[k][c] + f * (stops[k + 1][c] - stops[k][c]);
    }
    return (LUT[name] = out);
  }
  const css = v => getComputedStyle(document.body).getPropertyValue(v).trim();

  let cv, ctx, tipEl, legendEl, meta = null, off = null, offCtx = null;
  let state = {};                        // last render args, kept for hover readout
  let box = { x: 0, y: 0, w: 0, h: 0 };  // where the map sits inside the canvas

  function init(canvas, tip, legend) {
    cv = canvas; ctx = cv.getContext("2d"); tipEl = tip; legendEl = legend;
    cv.addEventListener("mousemove", hover);
    cv.addEventListener("mouseleave", () => tipEl.classList.add("hidden"));
  }

  function setMeta(m) {
    meta = m;
    off = document.createElement("canvas");
    off.width = m.nx; off.height = m.ny;
    offCtx = off.getContext("2d");
  }

  function fit() {
    const wrap = cv.parentElement, dpr = window.devicePixelRatio || 1;
    const latm = (meta.lat[0] + meta.lat[meta.lat.length - 1]) / 2;
    const dlon = meta.nx > 1 ? meta.lon[1] - meta.lon[0] : 1, dlat = meta.ny > 1 ? meta.lat[1] - meta.lat[0] : 1;
    const aspect = (meta.nx * Math.abs(dlon) * Math.cos(latm * Math.PI / 180)) / (meta.ny * Math.abs(dlat));
    let w = wrap.clientWidth, h = wrap.clientHeight;
    if (w / h > aspect) w = h * aspect; else h = w / aspect;
    cv.style.width = w + "px"; cv.style.height = h + "px";
    cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    box = { x: 0, y: 0, w, h };
    const pad = (cv.parentElement.clientWidth - w) / 2;             // keep the overlays on the map itself
    for (const [node, side] of [[legendEl, "left"], [document.getElementById("spLegend"), "right"]]) {
      if (node) node.style[side] = Math.max(pad + 8, 8) + "px";
    }
  }

  const toPx = (lon, lat) => {
    const dlon = meta.nx > 1 ? meta.lon[1] - meta.lon[0] : 1, dlat = meta.ny > 1 ? meta.lat[1] - meta.lat[0] : 1;
    const l0 = meta.lon[0] - dlon / 2, l1 = meta.lon[meta.nx - 1] + dlon / 2;
    const b0 = meta.lat[0] - dlat / 2, b1 = meta.lat[meta.ny - 1] + dlat / 2;
    return [box.x + (lon - l0) / (l1 - l0) * box.w, box.y + (1 - (lat - b0) / (b1 - b0)) * box.h];
  };
  const fromPx = (px, py) => {
    const dlon = meta.nx > 1 ? meta.lon[1] - meta.lon[0] : 1, dlat = meta.ny > 1 ? meta.lat[1] - meta.lat[0] : 1;
    const l0 = meta.lon[0] - dlon / 2, l1 = meta.lon[meta.nx - 1] + dlon / 2;
    const b0 = meta.lat[0] - dlat / 2, b1 = meta.lat[meta.ny - 1] + dlat / 2;
    return [l0 + (px - box.x) / box.w * (l1 - l0), b0 + (1 - (py - box.y) / box.h) * (b1 - b0)];
  };

  /* a: Float32Array(ny*nx) with NaN over land; row 0 = southernmost latitude */
  function render(a, opt) {
    if (!meta || !a) return;
    fit();
    state = Object.assign({}, opt, { field: a });
    const { lo, hi, log } = opt, table = lut(opt.cmap);
    const land = hex2rgb(rgbOf(css("--land")));
    const img = offCtx.createImageData(meta.nx, meta.ny);
    const llo = log ? Math.log10(Math.max(lo, 1e-30)) : lo, lhi = log ? Math.log10(Math.max(hi, 1e-29)) : hi;
    const span = (lhi - llo) || 1;
    for (let j = 0; j < meta.ny; j++) {
      const srcRow = (meta.ny - 1 - j) * meta.nx;                 // flip: north at top
      for (let i = 0; i < meta.nx; i++) {
        const v = a[srcRow + i], p = (j * meta.nx + i) * 4;
        if (!isFinite(v)) { img.data[p] = land[0]; img.data[p + 1] = land[1]; img.data[p + 2] = land[2]; img.data[p + 3] = 255; continue; }
        let t = ((log ? Math.log10(Math.max(v, 1e-30)) : v) - llo) / span;
        t = t < 0 ? 0 : t > 1 ? 1 : t;
        const k = (t * 255 | 0) * 3;
        img.data[p] = table[k]; img.data[p + 1] = table[k + 1]; img.data[p + 2] = table[k + 2]; img.data[p + 3] = 255;
      }
    }
    offCtx.putImageData(img, 0, 0);
    ctx.clearRect(0, 0, box.w, box.h);
    ctx.imageSmoothingEnabled = opt.smooth !== false;
    ctx.drawImage(off, box.x, box.y, box.w, box.h);

    if (opt.ibm) drawIbm(opt.ibm);
    if (opt.vectors) drawVectors(opt.vectors);
    if (opt.agents) drawAgents(opt.agents, opt.colors);
    drawLegend(opt);
  }

  function rgbOf(c) {                                  // accept #rrggbb or rgb(...)
    if (c.startsWith("#")) return c;
    const m = c.match(/(\d+)\D+(\d+)\D+(\d+)/);
    return m ? "#" + [m[1], m[2], m[3]].map(x => (+x).toString(16).padStart(2, "0")).join("") : "#888888";
  }

  function drawIbm(ibm) {
    const cw = box.w / meta.nx, ch = box.h / meta.ny;
    ctx.strokeStyle = "rgba(64,224,224,0.85)"; ctx.lineWidth = 1;
    for (let j = 0; j < meta.ny; j++) for (let i = 0; i < meta.nx; i++) {
      if (ibm[(meta.ny - 1 - j) * meta.nx + i] > 0.5) ctx.strokeRect(box.x + i * cw + .5, box.y + j * ch + .5, cw - 1, ch - 1);
    }
  }

  function drawVectors(v) {
    const { u, w, ny, nx, stride, scale } = v;
    const cw = box.w / meta.nx * stride, ch = box.h / meta.ny * stride;
    const L = Math.min(cw, ch) * 0.9;
    const arrows = new Path2D();
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const uu = u[j * nx + i], vv = w[j * nx + i];
      if (!isFinite(uu) || !isFinite(vv)) continue;
      const sp = Math.hypot(uu, vv);
      if (sp < 1e-6) continue;
      const [px, py] = toPx(meta.lon[Math.min(i * stride, meta.nx - 1)], meta.lat[Math.min(j * stride, meta.ny - 1)]);
      const f = Math.min(sp / scale, 1) * L, ang = Math.atan2(-vv, uu);
      const ex = px + Math.cos(ang) * f, ey = py + Math.sin(ang) * f;
      arrows.moveTo(px, py); arrows.lineTo(ex, ey);
      arrows.lineTo(ex - Math.cos(ang - 0.4) * f * 0.3, ey - Math.sin(ang - 0.4) * f * 0.3);
      arrows.moveTo(ex, ey);
      arrows.lineTo(ex - Math.cos(ang + 0.4) * f * 0.3, ey - Math.sin(ang + 0.4) * f * 0.3);
    }
    ctx.lineJoin = "round";                                  // dark backing keeps arrows legible on light cells
    ctx.strokeStyle = "rgba(0,0,0,0.45)"; ctx.lineWidth = 2.4; ctx.stroke(arrows);
    ctx.strokeStyle = "rgba(255,255,255,0.95)"; ctx.lineWidth = 1; ctx.stroke(arrows);
  }

  /* Dots scale with the cell size so a coarse grid does not turn into a blob;
     radius grows with the fish each super-individual represents. */
  function drawAgents(ag, colors) {
    const cell = Math.min(box.w / meta.nx, box.h / meta.ny);
    const base = Math.max(0.7, Math.min(cell * 0.10, 2.0)), rmax = Math.max(1.4, cell * 0.30);
    ctx.lineWidth = 0.8;
    for (let i = 0; i < ag.n; i++) {
      const st = colors[ag.species[i]];
      if (!st) continue;                                            // species filtered out
      const [px, py] = toPx(ag.lon[i], ag.lat[i]);
      const r = Math.min(base * (0.55 + 0.11 * Math.log10(Math.max(ag.num[i], 1))), rmax);
      ctx.beginPath();
      ctx.arc(px, py, st.ring ? Math.max(r, 1.6) : r, 0, 6.2832);
      ctx.globalAlpha = ag.stage[i] >= 2 ? 0.8 : 0.35;              // juveniles/adults solid, eggs/larvae faint
      if (st.ring) {                                                // 9th+ series: a ring of the same hue
        ctx.strokeStyle = st.color; ctx.lineWidth = 1.2; ctx.stroke();
      } else {
        ctx.fillStyle = st.color; ctx.fill();
      }
    }
    ctx.globalAlpha = 1;
  }

  function drawLegend(opt) {
    const stops = CMAPS[opt.cmap] || CMAPS.blues;
    const ticks = opt.log
      ? [opt.lo, Math.sqrt(opt.lo * opt.hi), opt.hi]
      : [opt.lo, (opt.lo + opt.hi) / 2, opt.hi];
    legendEl.innerHTML =
      `<div>${opt.label || ""}${opt.log ? " · log" : ""}</div>` +
      `<div class="bar" style="background:linear-gradient(90deg,${stops.join(",")})"></div>` +
      `<div class="ends">${ticks.map(t => Charts.fmt(t)).join("<span></span>")}</div>`;
  }

  function hover(ev) {
    if (!meta || !state.field) return;
    const r = cv.getBoundingClientRect();
    const [lon, lat] = fromPx(ev.clientX - r.left, ev.clientY - r.top);
    const dlon = meta.nx > 1 ? meta.lon[1] - meta.lon[0] : 1, dlat = meta.ny > 1 ? meta.lat[1] - meta.lat[0] : 1;
    const i = Math.round((lon - meta.lon[0]) / dlon), j = Math.round((lat - meta.lat[0]) / dlat);
    if (i < 0 || j < 0 || i >= meta.nx || j >= meta.ny) { tipEl.classList.add("hidden"); return; }
    const v = state.field[j * meta.nx + i];
    tipEl.classList.remove("hidden");
    tipEl.textContent = `${lat.toFixed(2)}°${lat >= 0 ? "N" : "S"}  ${lon.toFixed(2)}°${lon >= 0 ? "E" : "W"}\n` +
      (isFinite(v) ? `${state.label || "value"} = ${Charts.fmt(v, 4)}` : "land");
    const wrapRect = cv.parentElement.getBoundingClientRect();
    tipEl.style.left = (ev.clientX - wrapRect.left + 14) + "px";
    tipEl.style.top = (ev.clientY - wrapRect.top + 12) + "px";
  }

  return { init, setMeta, render, toPx, lut, stops: CMAPS, cmaps: Object.keys(CMAPS) };
})();
