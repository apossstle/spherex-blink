// SPHEREx Blink — browser side.
// The server sends aligned float32 frames; all display stretching, differencing,
// colour mapping and overlays happen here so the user can tune them instantly.

const DETECTOR_INFO = {
  1: { range: "0.75–1.09µm", color: "#6aa8ff" },
  2: { range: "1.10–1.62µm", color: "#52d1c4" },
  3: { range: "1.63–2.41µm", color: "#9be564" },
  4: { range: "2.42–3.82µm", color: "#ffd166" },
  5: { range: "3.83–4.41µm", color: "#ff9f5a" },
  6: { range: "4.42–5.00µm", color: "#ff6b8b" },
};
const VISIT_COLORS = ["#70d6ff", "#ffb454", "#c792ea", "#7ee787", "#ff6b8b", "#f6e05e", "#4fd1c5", "#fc8181"];

// Positions are the challenge's test targets; RA/Dec here are only used if
// the name resolver is unreachable (high proper-motion stars are given at
// ~2025.9 so they still land in the field).
const PRESETS = [
  { label: "Barnard's Star", query: "Barnard's Star", ra: 269.4463, dec: 4.7676, fov: 0.1 },
  { label: "ε Indi B", query: "eps Ind Ba", ra: 331.0959, dec: -56.8009, fov: 0.12 },
  { label: "WISE 0855", query: "WISE J085510.83-071442.5", ra: "08h55m10.831s", dec: "-7d14m42.53s", fov: 0.15 },
  { label: "RW Aurigae", query: "RW Aur", ra: 76.95642, dec: 30.40141, fov: 0.08 },
  { label: "R Draconis", query: "R Dra", ra: 248.16742, dec: 66.75497, fov: 0.08 },
];

const $ = (sel) => document.querySelector(sel);
const $$ = (sel) => [...document.querySelectorAll(sel)];
const DEG = Math.PI / 180;

const state = {
  data: null, // server response
  frames: [], // Float32Array per frame
  layers: {}, // detector -> Float32Array[]
  width: 0,
  height: 0,
  index: 0,
  playing: false,
  view: "blink",
  gaia: null,
  derivedKey: "",
  derived: null, // { frames: [...], rgb: [[r,g,b], ...], limits }
  timer: null,
};

// --------------------------------------------------------------- utilities
function toast(text) {
  const el = document.createElement("div");
  el.className = "toast";
  el.textContent = text;
  document.body.append(el);
  setTimeout(() => el.remove(), 2200);
}

function setStatus(text, kind = "") {
  const el = $("#status");
  el.textContent = text;
  el.className = `status ${kind}`;
}

async function api(path, params) {
  const url = `${path}?${new URLSearchParams(params)}`;
  const res = await fetch(url);
  let body = null;
  try { body = await res.json(); } catch { /* non-JSON error page */ }
  if (!res.ok) {
    const detail = body?.detail;
    throw new Error(typeof detail === "string" ? detail : `Server error ${res.status}`);
  }
  return body;
}

function decodeFrame(b64) {
  const bin = atob(b64);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return new Float32Array(bytes.buffer);
}

const mjdToYear = (mjd) => 2000 + (mjd - 51544.5) / 365.25;
const mjdToDate = (mjd) => new Date((mjd - 40587) * 86400000);
const fmtDate = (mjd) => mjdToDate(mjd).toISOString().slice(0, 10);

function percentile(sorted, p) {
  if (!sorted.length) return 0;
  const i = Math.min(sorted.length - 1, Math.max(0, Math.round((p / 100) * (sorted.length - 1))));
  return sorted[i];
}

function finiteSample(arrays, maxN = 60000) {
  const total = arrays.reduce((s, a) => s + a.length, 0);
  const step = Math.max(1, Math.floor(total / maxN));
  const out = [];
  let k = 0;
  for (const a of arrays) for (let i = 0; i < a.length; i++, k++) if (k % step === 0 && Number.isFinite(a[i])) out.push(a[i]);
  return Float64Array.from(out).sort();
}

function nanMean(arrays) {
  const n = arrays[0].length;
  const out = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    let s = 0, c = 0;
    for (const a of arrays) if (Number.isFinite(a[i])) { s += a[i]; c++; }
    out[i] = c ? s / c : NaN;
  }
  return out;
}

const subtract = (a, b) => a.map((v, i) => v - b[i]);

// ---------------------------------------------------------- sky projection
// Server grid: gnomonic (TAN), north up, east left, target on the centre pixel.
function skyToPixel(ra, dec) {
  const { ra: ra0, dec: dec0 } = state.data.target;
  const s = state.data.grid.pixel_scale_arcsec / 3600;
  const d = dec * DEG, d0 = dec0 * DEG, da = (ra - ra0) * DEG;
  const cosc = Math.sin(d0) * Math.sin(d) + Math.cos(d0) * Math.cos(d) * Math.cos(da);
  if (cosc <= 0) return null;
  const xi = (Math.cos(d) * Math.sin(da)) / cosc / DEG;
  const eta = (Math.cos(d0) * Math.sin(d) - Math.sin(d0) * Math.cos(d) * Math.cos(da)) / cosc / DEG;
  const cx = (state.width - 1) / 2, cy = (state.height - 1) / 2;
  return { x: cx - xi / s, y: cy + eta / s };
}

function pixelToSky(x, y) {
  const { ra: ra0, dec: dec0 } = state.data.target;
  const s = state.data.grid.pixel_scale_arcsec / 3600;
  const cx = (state.width - 1) / 2, cy = (state.height - 1) / 2;
  const xi = -(x - cx) * s * DEG, eta = (y - cy) * s * DEG;
  const rho = Math.hypot(xi, eta);
  const d0 = dec0 * DEG;
  if (rho === 0) return { ra: ra0, dec: dec0 };
  const c = Math.atan(rho);
  const dec = Math.asin(Math.cos(c) * Math.sin(d0) + (eta * Math.sin(c) * Math.cos(d0)) / rho);
  const ra = ra0 * DEG + Math.atan2(xi * Math.sin(c), rho * Math.cos(d0) * Math.cos(c) - eta * Math.sin(d0) * Math.sin(c));
  return { ra: ((ra / DEG) % 360 + 360) % 360, dec: dec / DEG };
}

function sexa(ra, dec) {
  const h = ra / 15, hh = Math.floor(h), mm = Math.floor((h - hh) * 60), ss = ((h - hh) * 60 - mm) * 60;
  const sg = dec < 0 ? "-" : "+", ad = Math.abs(dec), dd = Math.floor(ad), dm = Math.floor((ad - dd) * 60), ds = ((ad - dd) * 60 - dm) * 60;
  const p = (v, w = 2) => String(v).padStart(w, "0");
  return `${p(hh)}h${p(mm)}m${ss.toFixed(2).padStart(5, "0")}s ${sg}${p(dd)}°${p(dm)}′${ds.toFixed(1).padStart(4, "0")}″`;
}

// ------------------------------------------------------------- colour maps
function makeLut(stops) {
  const lut = new Uint8ClampedArray(256 * 3);
  for (let i = 0; i < 256; i++) {
    const t = i / 255;
    let k = 0;
    while (k < stops.length - 2 && t > stops[k + 1][0]) k++;
    const [t0, c0] = stops[k], [t1, c1] = stops[k + 1];
    const f = (t - t0) / (t1 - t0 || 1);
    for (let j = 0; j < 3; j++) lut[i * 3 + j] = c0[j] + (c1[j] - c0[j]) * f;
  }
  return lut;
}
const LUTS = {
  gray: makeLut([[0, [0, 0, 0]], [1, [255, 255, 255]]]),
  heat: makeLut([[0, [0, 0, 4]], [0.3, [87, 16, 110]], [0.55, [188, 55, 84]], [0.8, [249, 142, 9]], [1, [252, 255, 164]]]),
  diverge: makeLut([[0, [40, 120, 255]], [0.5, [12, 14, 22]], [1, [255, 170, 60]]]),
};
const STRETCH = {
  linear: (t) => t,
  sqrt: (t) => Math.sqrt(t),
  log: (t) => Math.log10(1 + 999 * t) / 3,
  asinh: (t) => Math.asinh(10 * t) / Math.asinh(10),
};

// --------------------------------------------------------- derived images
function isDiff(view) { return view === "diff-prev" || view === "diff-mean"; }

function deriveSet(source) {
  const view = state.view;
  if (view === "stack") return [nanMean(source)];
  if (view === "diff-mean") { const m = nanMean(source); return source.map((f) => subtract(f, m)); }
  if (view === "diff-prev") return source.map((f, i) => subtract(f, source[(i - 1 + source.length) % source.length]));
  return source;
}

function limitsFor(arrays) {
  const clip = Number($("#clip").value);
  const sample = finiteSample(arrays);
  if (isDiff(state.view)) {
    const abs = Float64Array.from(sample, Math.abs).sort();
    const hi = percentile(abs, clip) || 1;
    return { lo: -hi, hi };
  }
  const lo = percentile(sample, 5);
  const hi = percentile(sample, clip);
  return { lo, hi: hi > lo ? hi : lo + 1 };
}

function ensureDerived() {
  const color = $("#color").value;
  const rgbKeys = [$("#rgb-r").value, $("#rgb-g").value, $("#rgb-b").value];
  const key = [state.view, color, rgbKeys, $("#clip").value, state.frames.length].join("|");
  if (state.derivedKey === key && state.derived) return state.derived;

  const useRgb = color === "rgb" && Object.keys(state.layers).length;
  let derived;
  if (useRgb) {
    const channels = rgbKeys.map((d) => deriveSet(state.layers[d] || state.frames));
    const limits = channels.map(limitsFor);
    derived = { rgb: true, channels, limits, count: channels[0].length };
  } else {
    const frames = deriveSet(state.frames);
    derived = { rgb: false, frames, limits: limitsFor(frames), count: frames.length };
  }
  state.derived = derived;
  state.derivedKey = key;
  return derived;
}

// --------------------------------------------------------------- rendering
const imgCanvas = $("#img");
const overlay = $("#overlay");

function renderImage() {
  if (!state.data) return;
  const d = ensureDerived();
  const i = Math.min(state.index, d.count - 1);
  const { width: w, height: h } = state;
  if (imgCanvas.width !== w || imgCanvas.height !== h) { imgCanvas.width = w; imgCanvas.height = h; }
  const ctx = imgCanvas.getContext("2d");
  const out = ctx.createImageData(w, h);
  const px = out.data;
  const stretch = STRETCH[$("#stretch").value];
  const invert = $("#invert").checked;
  const diff = isDiff(state.view);

  const scale = (v, lim) => {
    if (!Number.isFinite(v)) return NaN;
    if (diff) {
      const t = Math.min(1, Math.abs(v) / lim.hi);
      return 0.5 + 0.5 * Math.sign(v) * stretch(t);
    }
    return stretch(Math.min(1, Math.max(0, (v - lim.lo) / (lim.hi - lim.lo))));
  };

  const lut = diff && $("#color").value !== "gray" ? LUTS.diverge : LUTS[$("#color").value] || LUTS.gray;
  for (let y = 0; y < h; y++) {
    const srcRow = (h - 1 - y) * w; // FITS row 0 is the bottom of the image
    for (let x = 0; x < w; x++) {
      const o = (y * w + x) * 4, s = srcRow + x;
      if (d.rgb) {
        let any = false;
        for (let c = 0; c < 3; c++) {
          let t = scale(d.channels[c][i][s], d.limits[c]);
          if (Number.isFinite(t)) any = true; else t = 0;
          if (invert) t = 1 - t;
          px[o + c] = t * 255;
        }
        if (!any) { px[o] = 10; px[o + 1] = 16; px[o + 2] = 32; }
      } else {
        let t = scale(d.frames[i][s], d.limits);
        if (!Number.isFinite(t)) { px[o] = 10; px[o + 1] = 16; px[o + 2] = 32; px[o + 3] = 255; continue; }
        if (invert) t = 1 - t;
        const k = Math.round(t * 255) * 3;
        px[o] = lut[k]; px[o + 1] = lut[k + 1]; px[o + 2] = lut[k + 2];
      }
      px[o + 3] = 255;
    }
  }
  ctx.putImageData(out, 0, 0);
  renderOverlay();
  updateBadge();
}

function currentMjd() {
  const frames = state.data.frames;
  if (state.view === "stack") return frames.reduce((s, f) => s + f.mjd, 0) / frames.length;
  return frames[Math.min(state.index, frames.length - 1)].mjd;
}

function renderOverlay(target = overlay, sizeOverride) {
  const rect = $("#viewer").getBoundingClientRect();
  const size = sizeOverride || Math.round(rect.width * devicePixelRatio);
  if (target.width !== size) { target.width = size; target.height = size; }
  const ctx = target.getContext("2d");
  ctx.clearRect(0, 0, size, size);
  if (!state.data) return;
  const k = size / state.width;
  const toScreen = (p) => ({ x: (p.x + 0.5) * k, y: (state.height - 0.5 - p.y) * k });
  const lw = Math.max(1, size / 500);

  if ($("#crosshair").checked) {
    const c = size / 2, r = 3.2 * k;
    ctx.strokeStyle = "#70d6ffcc";
    ctx.lineWidth = 1.5 * lw;
    ctx.beginPath();
    for (const [dx, dy] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) {
      ctx.moveTo(c + dx * r, c + dy * r);
      ctx.lineTo(c + dx * r * 2.2, c + dy * r * 2.2);
    }
    ctx.stroke();
  }

  if ($("#gaia").checked && state.gaia) {
    const year = mjdToYear(currentMjd());
    ctx.font = `${11 * lw}px system-ui`;
    for (const s of state.gaia) {
      const dt = year - 2016.0;
      const pm = Math.hypot(s.pmra || 0, s.pmdec || 0);
      const ra = s.ra + ((s.pmra || 0) * dt) / 3.6e6 / Math.cos(s.dec * DEG);
      const dec = s.dec + ((s.pmdec || 0) * dt) / 3.6e6;
      const p = skyToPixel(ra, dec);
      if (!p || p.x < -2 || p.y < -2 || p.x > state.width + 1 || p.y > state.height + 1) continue;
      const q = toScreen(p);
      const fast = pm > 100;
      const r = Math.max(2.2, 4.5 - (s.g ?? 18) / 6) * k * 0.6;
      ctx.strokeStyle = fast ? "#ffb454" : "#7ee787aa";
      ctx.lineWidth = (fast ? 2 : 1) * lw;
      ctx.beginPath();
      ctx.arc(q.x, q.y, r, 0, 2 * Math.PI);
      ctx.stroke();
      if (fast) {
        // Arrow: where this star will be in 5 years.
        const ahead = skyToPixel(ra + ((s.pmra || 0) * 5) / 3.6e6 / Math.cos(s.dec * DEG), dec + ((s.pmdec || 0) * 5) / 3.6e6);
        if (ahead) {
          const a = toScreen(ahead);
          ctx.beginPath(); ctx.moveTo(q.x, q.y); ctx.lineTo(a.x, a.y); ctx.stroke();
        }
        ctx.fillStyle = "#ffb454";
        ctx.fillText(`${(pm / 1000).toFixed(2)}″/yr`, q.x + r + 3, q.y - r);
      }
    }
  }

  // Compass and 1′ scale bar.
  ctx.strokeStyle = ctx.fillStyle = "#e6ecf7cc";
  ctx.lineWidth = 1.5 * lw;
  ctx.font = `${12 * lw}px system-ui`;
  const m = 14 * lw, L = 26 * lw, ox = size - m - 8 * lw, oy = size - m - 6 * lw;
  ctx.beginPath(); ctx.moveTo(ox, oy); ctx.lineTo(ox, oy - L); ctx.moveTo(ox, oy); ctx.lineTo(ox - L, oy); ctx.stroke();
  ctx.fillText("N", ox - 4 * lw, oy - L - 4 * lw);
  ctx.fillText("E", ox - L - 12 * lw, oy + 4 * lw);
  const arcminPx = (60 / state.data.grid.pixel_scale_arcsec) * k;
  if (arcminPx < size * 0.5) {
    const bx = size - m - 8 * lw - L - 30 * lw - arcminPx, by = size - m;
    ctx.beginPath(); ctx.moveTo(bx, by); ctx.lineTo(bx + arcminPx, by); ctx.stroke();
    ctx.fillText("1′", bx + arcminPx / 2 - 4 * lw, by - 5 * lw);
  }
}

function updateBadge() {
  const badge = $("#frame-badge");
  const frames = state.data.frames;
  if (state.view === "stack") {
    const n = frames.reduce((s, f) => s + f.n_exposures, 0);
    badge.innerHTML = `Deep stack · ${frames.length} frames · ${n} exposures`;
  } else {
    const f = frames[state.index];
    const unit = state.data.mode === "visit" ? "Visit" : "Exposure";
    const warn = f.flagged_fraction > 0.25 ? ` · <span class="warn">⚠ target flagged in ${Math.round(f.flagged_fraction * 100)}% of exposures</span>` : "";
    const diff = state.view === "diff-prev" ? " − previous" : state.view === "diff-mean" ? " − mean" : "";
    badge.innerHTML = `${unit} ${state.index + 1}/${frames.length}${diff} · <b>${f.date.slice(0, 10)}</b> · ${f.n_exposures} exp.${warn}`;
  }
  badge.hidden = false;
  $$("#frames-table tbody tr").forEach((tr, i) => tr.classList.toggle("current", i === state.index));
}

// ------------------------------------------------------------------ player
function setIndex(i) {
  const n = state.view === "stack" ? 1 : state.frames.length;
  state.index = ((i % n) + n) % n;
  $("#frame").value = state.index;
  renderImage();
}

function play(on) {
  state.playing = on;
  $("#btn-play").textContent = on ? "⏸" : "▶";
  $("#btn-play").setAttribute("aria-label", on ? "Pause" : "Play");
  clearTimeout(state.timer);
  if (!on) return;
  const tick = () => {
    setIndex(state.index + 1);
    state.timer = setTimeout(tick, 1000 / Number($("#speed").value));
  };
  state.timer = setTimeout(tick, 1000 / Number($("#speed").value));
}

// ------------------------------------------------------------------ charts
function setupCanvas(canvas) {
  const w = canvas.clientWidth, h = canvas.clientHeight;
  canvas.width = w * devicePixelRatio; canvas.height = h * devicePixelRatio;
  const ctx = canvas.getContext("2d");
  ctx.setTransform(devicePixelRatio, 0, 0, devicePixelRatio, 0, 0);
  ctx.clearRect(0, 0, w, h);
  return { ctx, w, h };
}

function niceTicks(lo, hi, n = 5) {
  const span = hi - lo || 1;
  const step0 = span / n, mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => span / s <= n) || 10 * mag;
  const ticks = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) ticks.push(+v.toFixed(10));
  return ticks;
}

function scatterPlot(canvas, { points, xLabel, yLabel, xFormat = (v) => v, legend = [] }) {
  const { ctx, w, h } = setupCanvas(canvas);
  const pad = { l: 62, r: 12, t: 22, b: 36 };
  if (!points.length) {
    ctx.fillStyle = "#8796b0"; ctx.font = "13px system-ui";
    ctx.fillText("No measurements at the target position.", pad.l, h / 2);
    return;
  }
  const xs = points.map((p) => p.x), ys = points.flatMap((p) => [p.y - (p.e || 0), p.y + (p.e || 0)]);
  const sortedY = Float64Array.from(ys.filter(Number.isFinite)).sort();
  let x0 = Math.min(...xs), x1 = Math.max(...xs);
  let y0 = percentile(sortedY, 1), y1 = percentile(sortedY, 99);
  if (x1 === x0) { x0 -= 0.5; x1 += 0.5; }
  if (y1 === y0) { y0 -= 1; y1 += 1; }
  const yPad = (y1 - y0) * 0.08; y0 -= yPad; y1 += yPad;
  const xPad = (x1 - x0) * 0.03; x0 -= xPad; x1 += xPad;
  const X = (v) => pad.l + ((v - x0) / (x1 - x0)) * (w - pad.l - pad.r);
  const Y = (v) => h - pad.b - ((v - y0) / (y1 - y0)) * (h - pad.t - pad.b);

  ctx.font = "11px system-ui"; ctx.strokeStyle = "#22304a"; ctx.fillStyle = "#8796b0"; ctx.lineWidth = 1;
  for (const t of niceTicks(y0, y1)) {
    ctx.beginPath(); ctx.moveTo(pad.l, Y(t)); ctx.lineTo(w - pad.r, Y(t)); ctx.stroke();
    ctx.fillText(new Intl.NumberFormat("en", { notation: "compact", maximumFractionDigits: 2 }).format(t), 8, Y(t) + 4);
  }
  for (const t of niceTicks(x0, x1, 7)) {
    ctx.beginPath(); ctx.moveTo(X(t), pad.t); ctx.lineTo(X(t), h - pad.b); ctx.stroke();
    ctx.fillText(xFormat(t), X(t) - 12, h - pad.b + 14);
  }
  if (y0 < 0 && y1 > 0) { ctx.strokeStyle = "#3a4d72"; ctx.beginPath(); ctx.moveTo(pad.l, Y(0)); ctx.lineTo(w - pad.r, Y(0)); ctx.stroke(); }
  ctx.fillStyle = "#c5d0e3";
  ctx.fillText(xLabel, w / 2 - 30, h - 6);
  ctx.save(); ctx.translate(12, h / 2 + 30); ctx.rotate(-Math.PI / 2); ctx.fillText(yLabel, 0, 0); ctx.restore();

  for (const p of points) {
    if (!Number.isFinite(p.y)) continue;
    const cx = X(p.x), cy = Y(Math.min(y1, Math.max(y0, p.y)));
    ctx.strokeStyle = ctx.fillStyle = p.c;
    ctx.globalAlpha = p.flagged ? 0.35 : 0.9;
    if (p.e) { ctx.beginPath(); ctx.moveTo(cx, Y(Math.max(y0, p.y - p.e))); ctx.lineTo(cx, Y(Math.min(y1, p.y + p.e))); ctx.stroke(); }
    ctx.beginPath(); ctx.arc(cx, cy, 2.6, 0, 2 * Math.PI); ctx.fill();
  }
  ctx.globalAlpha = 1;
  let lx = pad.l + 4;
  for (const { label, color } of legend) {
    ctx.fillStyle = color; ctx.fillRect(lx, 7, 9, 9);
    ctx.fillStyle = "#c5d0e3"; ctx.fillText(label, lx + 13, 15);
    lx += ctx.measureText(label).width + 28;
  }
}

const isDefect = (flags) => flags.some((f) => f !== "SOURCE");

function renderCharts() {
  if (!state.data) return;
  const spec = state.data.spectrum;
  const visits = [...new Set(spec.map((p) => p.visit))].sort((a, b) => a - b);
  scatterPlot($("#spectrum"), {
    points: spec.map((p) => ({ x: p.wavelength_um, y: p.flux_ujy, e: p.error_ujy, c: VISIT_COLORS[p.visit % VISIT_COLORS.length], flagged: isDefect(p.flags) })),
    xLabel: "Wavelength (µm)", yLabel: "Flux (µJy)", xFormat: (v) => v.toFixed(1),
    legend: visits.map((v) => {
      const first = spec.find((p) => p.visit === v);
      return { label: `Visit ${v + 1} · ${fmtDate(first.mjd)}`, color: VISIT_COLORS[v % VISIT_COLORS.length] };
    }),
  });
  scatterPlot($("#lightcurve"), {
    points: spec.map((p) => ({ x: mjdToYear(p.mjd), y: p.flux_ujy, e: p.error_ujy, c: DETECTOR_INFO[p.detector].color, flagged: isDefect(p.flags) })),
    xLabel: "Year", yLabel: "Flux (µJy)", xFormat: (v) => v.toFixed(2),
    legend: [...new Set(spec.map((p) => p.detector))].sort().map((d) => ({ label: `D${d}`, color: DETECTOR_INFO[d].color })),
  });
}

function renderTable() {
  const body = $("#frames-table tbody");
  body.innerHTML = "";
  state.data.frames.forEach((f, i) => {
    const tr = document.createElement("tr");
    const flags = f.flags_at_target.map((n) => `<span class="flag ${n === "SOURCE" ? "src" : ""}">${n}</span>`).join("") || "—";
    tr.innerHTML = `<td>${i + 1}</td><td>${f.date.replace("T", " ").slice(0, 16)}</td><td>${f.n_exposures}</td><td>${f.detectors.map((d) => "D" + d).join(" ")}</td><td>${Math.round(f.coverage * 100)}%</td><td>${flags}</td>`;
    tr.addEventListener("click", () => { play(false); if (state.view === "stack") setView("blink"); setIndex(i); });
    body.append(tr);
  });
}

function renderMeta() {
  const d = state.data;
  const { ra, dec } = d.target;
  const first = d.frames[0], last = d.frames[d.frames.length - 1];
  const rows = [
    ["Target", sexa(ra, dec)],
    ["Degrees", `${ra.toFixed(5)}, ${dec.toFixed(5)}`],
    ["Field", `${(d.fov_deg * 60).toFixed(1)}′ (${d.grid.width}×${d.grid.height} px at ${d.grid.pixel_scale_arcsec}″)`],
    ["Exposures", `${d.n_exposures} in ${d.n_visits} survey visit(s)`],
    ["Time span", `${fmtDate(first.mjd_start)} → ${fmtDate(last.mjd_end)}`],
    ["Arrays", d.detectors.map((x) => "D" + x).join(" ")],
    ["Releases", d.collections.join(", ")],
  ];
  if (d.warnings.length) rows.push(["Notes", d.warnings.join(" ")]);
  $("#meta").innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join("");
}

// -------------------------------------------------------------- view state
function setView(view) {
  state.view = view;
  $$(".seg button").forEach((b) => { b.classList.toggle("on", b.dataset.view === view); b.setAttribute("aria-checked", b.dataset.view === view); });
  const single = view === "stack";
  for (const id of ["#btn-prev", "#btn-next", "#btn-play", "#frame"]) $(id).disabled = single;
  if (single) { play(false); state.index = 0; }
  state.derivedKey = "";
  if (state.data) renderImage();
  syncUrl();
}

function readForm() {
  return {
    ra: $("#ra").value.trim(),
    dec: $("#dec").value.trim(),
    fov: Number($("#fov").value),
    mode: $("input[name=mode]:checked").value,
    detectors: $$("#detectors input:checked").map((i) => i.value).join(","),
  };
}

function syncUrl() {
  const f = readForm();
  if (!f.ra || !f.dec) return;
  const p = new URLSearchParams({ ra: f.ra, dec: f.dec, fov: f.fov, mode: f.mode });
  if (f.detectors !== "1,2,3,4,5,6") p.set("det", f.detectors);
  if (state.view !== "blink") p.set("view", state.view);
  if ($("#color").value !== "gray") p.set("color", $("#color").value);
  if ($("#gaia").checked) p.set("gaia", "1");
  if (document.body.classList.contains("embed")) p.set("embed", "1");
  history.replaceState(null, "", `${location.pathname}?${p}`);
}

// -------------------------------------------------------------- data load
async function load() {
  const form = readForm();
  if (!form.ra || !form.dec) { setStatus("Enter RA and Dec, or a name.", "error"); return; }
  if (!form.detectors) { setStatus("Select at least one detector array.", "error"); return; }
  play(false);
  $("#btn-load").disabled = true;
  $("#viewer").classList.add("loading");
  const started = performance.now();
  const slowNote = setTimeout(() => setStatus("Fetching cutouts from IRSA — the first look at a new field can take a minute or two…"), 4000);
  setStatus("Searching the SPHEREx archive…");
  syncUrl();
  try {
    const data = await api("/api/sequence", form);
    state.data = data;
    state.width = data.grid.width;
    state.height = data.grid.height;
    state.frames = data.frames.map((f) => decodeFrame(f.data));
    state.layers = Object.fromEntries(Object.entries(data.layers).map(([k, v]) => [k, v.map((f) => decodeFrame(f.data))]));
    state.index = 0;
    state.derivedKey = "";
    state.gaia = null;
    $("#frame").max = Math.max(0, state.frames.length - 1);
    $("#empty").hidden = true;
    $("#color option[value=rgb]").disabled = !Object.keys(state.layers).length;
    if ($("#color").value === "rgb" && !Object.keys(state.layers).length) $("#color").value = "gray";
    renderTable();
    renderMeta();
    renderImage();
    renderCharts();
    centerAladin();
    if ($("#gaia").checked) loadGaia();
    const secs = ((performance.now() - started) / 1000).toFixed(1);
    setStatus(`Loaded ${data.frames.length} frame(s) from ${data.n_exposures} exposures in ${secs} s.`, "ok");
    if (data.frames.length > 1 && state.view !== "stack") play(true);
  } catch (err) {
    setStatus(err.message, "error");
  } finally {
    clearTimeout(slowNote);
    $("#btn-load").disabled = false;
    $("#viewer").classList.remove("loading");
  }
}

async function loadGaia() {
  if (!state.data || state.gaia) { renderOverlay(); return; }
  const { ra, dec } = state.data.target;
  try {
    const res = await api("/api/gaia", { ra, dec, radius: Math.min(0.6, state.data.fov_deg * 0.75 + 0.01) });
    state.gaia = res.sources;
    renderOverlay();
  } catch (err) {
    toast(`Gaia overlay unavailable: ${err.message}`);
    $("#gaia").checked = false;
  }
}

async function resolveName(query, fallback) {
  setStatus(`Looking up “${query}”…`);
  try {
    const r = await api("/api/resolve", { name: query });
    $("#ra").value = r.ra.toFixed(6);
    $("#dec").value = r.dec.toFixed(6);
    const pm = r.pmra != null ? ` (moved for proper motion to ${r.epoch})` : "";
    setStatus(`${r.name || query}: ${sexa(r.ra, r.dec)}${pm}`, "ok");
    return true;
  } catch (err) {
    if (fallback) {
      $("#ra").value = fallback.ra;
      $("#dec").value = fallback.dec;
      setStatus(`Name service unavailable; using stored coordinates for ${query}.`);
      return true;
    }
    setStatus(err.message, "error");
    return false;
  }
}

// ------------------------------------------------------------------ Aladin
let aladin = null;
function initAladin() {
  if (!window.A || aladin) return;
  A.init.then(() => {
    aladin = A.aladin("#aladin", {
      survey: $("#hips").value, fov: 0.4, showReticle: true, showZoomControl: true,
      showLayersControl: false, showGotoControl: false, showFullscreenControl: false, showFrame: false,
      showCooGrid: false, target: "0 0",
    });
    aladin.on?.("click", (e) => {
      if (!e || e.ra == null) return;
      $("#ra").value = e.ra.toFixed(6);
      $("#dec").value = e.dec.toFixed(6);
      setStatus("New target picked on the map — press Load.");
    });
    centerAladin();
  }).catch(() => { $("#aladin").textContent = "Sky context unavailable."; });
}

let aladinOverlay = null;
function centerAladin() {
  if (!aladin) return;
  try {
    const f = readForm();
    const ra = state.data?.target.ra ?? Number(f.ra), dec = state.data?.target.dec ?? Number(f.dec);
    if (!Number.isFinite(ra) || !Number.isFinite(dec)) return;
    aladin.gotoRaDec(ra, dec);
    const fov = state.data?.fov_deg ?? f.fov;
    aladin.setFoV(Math.max(0.15, fov * 3));
    if (!aladinOverlay) { aladinOverlay = A.graphicOverlay({ color: "#ffb454", lineWidth: 2 }); aladin.addOverlay(aladinOverlay); }
    aladinOverlay.removeAll();
    const h = fov / 2, c = Math.cos(dec * DEG);
    aladinOverlay.add(A.polygon([[ra - h / c, dec - h], [ra + h / c, dec - h], [ra + h / c, dec + h], [ra - h / c, dec + h]]));
  } catch { /* Aladin API differences between versions are non-fatal */ }
}

// ------------------------------------------------------------------ export
async function recordVideo() {
  if (!state.data) return;
  if (!window.MediaRecorder) { toast("Video recording is not supported in this browser."); return; }
  const size = 720;
  const out = document.createElement("canvas");
  out.width = out.height = size;
  const ctx = out.getContext("2d");
  const ov = document.createElement("canvas");
  const stream = out.captureStream(30);
  const type = ["video/webm;codecs=vp9", "video/webm", "video/mp4"].find((t) => MediaRecorder.isTypeSupported(t));
  const rec = new MediaRecorder(stream, { mimeType: type });
  const chunks = [];
  rec.ondataavailable = (e) => e.data.size && chunks.push(e.data);
  const wasPlaying = state.playing;
  play(false);
  const n = state.view === "stack" ? 1 : state.frames.length;
  const loops = n < 4 ? 3 : 2;
  const hold = 1000 / Number($("#speed").value);
  rec.start();
  for (let k = 0; k < n * loops; k++) {
    setIndex(k % n);
    renderOverlay(ov, size);
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(imgCanvas, 0, 0, size, size);
    ctx.drawImage(ov, 0, 0);
    ctx.fillStyle = "#000b"; ctx.fillRect(0, 0, size, 30);
    ctx.fillStyle = "#e6ecf7"; ctx.font = "16px system-ui";
    ctx.fillText(`SPHEREx · ${sexa(state.data.target.ra, state.data.target.dec)} · ${fmtDate(currentMjd())}`, 10, 21);
    await new Promise((r) => setTimeout(r, hold));
  }
  rec.stop();
  await new Promise((r) => (rec.onstop = r));
  const blob = new Blob(chunks, { type: rec.mimeType });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `spherex_${state.data.target.ra.toFixed(4)}_${state.data.target.dec.toFixed(4)}.${rec.mimeType.includes("mp4") ? "mp4" : "webm"}`;
  a.click();
  if (wasPlaying) play(true);
}

// ---------------------------------------------------------------- bindings
function buildControls() {
  const presetBox = $("#presets");
  for (const p of PRESETS) {
    const b = document.createElement("button");
    b.type = "button";
    b.textContent = p.label;
    b.addEventListener("click", async () => {
      $("#name").value = p.query;
      $("#fov").value = p.fov; $("#fov-out").textContent = `${p.fov.toFixed(2)}°`;
      if (await resolveName(p.query, p)) load();
    });
    presetBox.append(b);
  }

  const detBox = $("#detectors");
  for (const [d, info] of Object.entries(DETECTOR_INFO)) {
    const l = document.createElement("label");
    l.innerHTML = `<input type="checkbox" value="${d}" checked><span class="swatch" style="background:${info.color}"></span>D${d} <small class="muted">${info.range}</small>`;
    detBox.append(l);
  }
  for (const [id, def] of [["#rgb-r", 5], ["#rgb-g", 3], ["#rgb-b", 1]]) {
    $(id).innerHTML = Object.keys(DETECTOR_INFO).map((d) => `<option value="${d}" ${+d === def ? "selected" : ""}>D${d}</option>`).join("");
  }

  $("#target-form").addEventListener("submit", (e) => { e.preventDefault(); load(); });
  $("#btn-resolve").addEventListener("click", () => { const q = $("#name").value.trim(); if (q) resolveName(q).then(centerAladin); });
  $("#name").addEventListener("keydown", (e) => {
    if (e.key === "Enter") { e.preventDefault(); const q = $("#name").value.trim(); if (q) resolveName(q).then((ok) => ok && load()); }
  });
  $("#fov").addEventListener("input", () => { $("#fov-out").textContent = `${Number($("#fov").value).toFixed(2)}°`; });

  $("#btn-play").addEventListener("click", () => play(!state.playing));
  $("#btn-prev").addEventListener("click", () => { play(false); setIndex(state.index - 1); });
  $("#btn-next").addEventListener("click", () => { play(false); setIndex(state.index + 1); });
  $("#frame").addEventListener("input", (e) => { play(false); setIndex(Number(e.target.value)); });
  $("#speed").addEventListener("input", () => state.playing && play(true));

  $$(".seg button").forEach((b) => b.addEventListener("click", () => setView(b.dataset.view)));
  for (const id of ["#stretch", "#invert", "#clip", "#rgb-r", "#rgb-g", "#rgb-b"]) {
    $(id).addEventListener("input", () => { state.derivedKey = ""; renderImage(); });
  }
  $("#color").addEventListener("change", () => {
    $("#rgb-pick").hidden = $("#color").value !== "rgb";
    state.derivedKey = ""; renderImage(); syncUrl();
  });
  $("#crosshair").addEventListener("change", renderOverlay);
  $("#gaia").addEventListener("change", () => { $("#gaia").checked ? loadGaia() : renderOverlay(); syncUrl(); });
  $("#hips").addEventListener("change", () => aladin?.setImageSurvey($("#hips").value));
  $("#btn-record").addEventListener("click", recordVideo);

  $$(".tabs button").forEach((b) => b.addEventListener("click", () => {
    $$(".tabs button").forEach((x) => x.classList.toggle("on", x === b));
    $$(".tab-body").forEach((x) => (x.hidden = x.dataset.body !== b.dataset.tab));
    renderCharts();
  }));

  $("#btn-help").addEventListener("click", () => {
    const h = $("#help"); h.hidden = !h.hidden; $("#btn-help").setAttribute("aria-expanded", !h.hidden);
  });
  $("#btn-share").addEventListener("click", async () => {
    syncUrl();
    await navigator.clipboard?.writeText(location.href).catch(() => {});
    toast("Link copied");
  });
  $("#btn-embed").addEventListener("click", async () => {
    syncUrl();
    const u = new URL(location.href); u.searchParams.set("embed", "1");
    const snippet = `<iframe src="${u}" width="520" height="560" style="border:0" title="SPHEREx Blink"></iframe>`;
    await navigator.clipboard?.writeText(snippet).catch(() => {});
    toast("Embed snippet copied");
  });

  // Hover readout and click-to-retarget on the image.
  const viewer = $("#viewer");
  const eventToPixel = (e) => {
    const r = viewer.getBoundingClientRect();
    const k = state.width / r.width;
    return { x: (e.clientX - r.left) * k - 0.5, y: state.height - 0.5 - (e.clientY - r.top) * k };
  };
  viewer.addEventListener("mousemove", (e) => {
    if (!state.data) return;
    const p = eventToPixel(e);
    const sky = pixelToSky(p.x, p.y);
    const d = ensureDerived();
    const ix = Math.round(p.x), iy = Math.round(p.y);
    let val = "";
    if (!d.rgb && ix >= 0 && iy >= 0 && ix < state.width && iy < state.height) {
      const v = d.frames[Math.min(state.index, d.count - 1)][iy * state.width + ix];
      val = Number.isFinite(v) ? ` · ${v.toPrecision(3)} MJy/sr` : " · no data";
    }
    const ro = $("#readout"); ro.hidden = false; ro.textContent = `${sexa(sky.ra, sky.dec)}${val}`;
  });
  viewer.addEventListener("mouseleave", () => ($("#readout").hidden = true));
  viewer.addEventListener("dblclick", (e) => {
    if (!state.data) return;
    const sky = pixelToSky(...Object.values(eventToPixel(e)));
    $("#ra").value = sky.ra.toFixed(6); $("#dec").value = sky.dec.toFixed(6);
    load();
  });

  document.addEventListener("keydown", (e) => {
    if (e.target.matches("input[type=text], select, textarea")) return;
    if (e.key === " ") { e.preventDefault(); play(!state.playing); }
    else if (e.key === "ArrowRight") { play(false); setIndex(state.index + 1); }
    else if (e.key === "ArrowLeft") { play(false); setIndex(state.index - 1); }
    else if (e.key.toLowerCase() === "d") {
      const views = ["blink", "diff-prev", "diff-mean", "stack"];
      setView(views[(views.indexOf(state.view) + 1) % views.length]);
    } else if (e.key.toLowerCase() === "g") { $("#gaia").checked = !$("#gaia").checked; $("#gaia").dispatchEvent(new Event("change")); }
  });

  new ResizeObserver(() => { if (state.data) { renderOverlay(); renderCharts(); } }).observe(viewer);
}

function applyUrl() {
  const p = new URLSearchParams(location.search);
  if (p.get("embed") === "1") document.body.classList.add("embed");
  if (p.has("fov")) { $("#fov").value = p.get("fov"); $("#fov-out").textContent = `${Number(p.get("fov")).toFixed(2)}°`; }
  if (p.get("mode") === "exposure") $("input[name=mode][value=exposure]").checked = true;
  if (p.has("det")) { const set = new Set(p.get("det").split(",")); $$("#detectors input").forEach((i) => (i.checked = set.has(i.value))); }
  if (p.has("color")) { $("#color").value = p.get("color"); $("#rgb-pick").hidden = p.get("color") !== "rgb"; }
  if (p.get("gaia") === "1") $("#gaia").checked = true;
  if (p.has("view")) setView(p.get("view"));
  if (p.has("ra") && p.has("dec")) {
    $("#ra").value = p.get("ra"); $("#dec").value = p.get("dec");
    load();
  } else if (p.has("target")) {
    $("#name").value = p.get("target");
    resolveName(p.get("target")).then((ok) => ok && load());
  }
}

buildControls();
applyUrl();
if (window.A) initAladin(); else window.addEventListener("load", initAladin);
