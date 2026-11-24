(() => {
  "use strict";

  const canvas = document.getElementById("flow-canvas");
  const ctx = canvas.getContext("2d");
  const chart = document.getElementById("convergence-chart");
  const chartCtx = chart.getContext("2d");
  const byId = (id) => document.getElementById(id);
  const q = [0, 1, 2, 3, 4, 5, 6, 7, 8];
  const cx = [0, 1, 0, -1, 0, 1, -1, -1, 1];
  const cy = [0, 0, 1, 0, -1, 1, 1, -1, -1];
  const opposite = [0, 3, 4, 1, 2, 7, 8, 5, 6];
  const weights = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36];
  const fields = {
    velocity: { label: "VELOCITY MAGNITUDE", unit: "lu/ts" },
    pressure: { label: "PRESSURE", unit: "lattice units" },
    temperature: { label: "TEMPERATURE", unit: "normalized" }
  };
  let state = { geometry: "cylinder", length: 4, blockage: .25, velocity: .06, viscosity: .04, diffusivity: .015, field: "velocity", streamlines: true, mesh: false, solved: false, obstacleX: .31 };
  let solver = null;
  let isRunning = false;
  let dragging = false;
  let lastResults = null;

  class LatticeSolver {
    constructor(nx, ny, options) {
      this.nx = nx;
      this.ny = ny;
      this.size = nx * ny;
      this.uIn = options.velocity;
      this.tau = .5 + 3 * options.viscosity;
      this.diffusivity = options.diffusivity;
      this.geometry = options.geometry;
      this.blockage = options.blockage;
      this.obstacleX = options.obstacleX;
      this.f = Array.from({ length: 9 }, () => new Float32Array(this.size));
      this.next = Array.from({ length: 9 }, () => new Float32Array(this.size));
      this.solid = new Uint8Array(this.size);
      this.ux = new Float32Array(this.size);
      this.uy = new Float32Array(this.size);
      this.rho = new Float32Array(this.size);
      this.temp = new Float32Array(this.size);
      this.tempNext = new Float32Array(this.size);
      this.makeGeometry();
      this.initialize();
    }

    makeGeometry() {
      const { nx, ny } = this;
      const diameter = Math.max(4, Math.round(ny * this.blockage));
      this.diameter = diameter;
      this.centerX = Math.max(diameter + 4, Math.min(nx - diameter - 5, Math.round(nx * this.obstacleX)));
      this.centerY = Math.floor(ny / 2);
      for (let y = 0; y < ny; y++) {
        for (let x = 0; x < nx; x++) {
          const index = y * nx + x;
          if (y === 0 || y === ny - 1) this.solid[index] = 1;
          if (this.geometry === "cylinder" || this.geometry === "square") {
            const dx = x - this.centerX;
            const dy = y - this.centerY;
            if (this.geometry === "cylinder" ? dx * dx + dy * dy <= (diameter / 2) ** 2 : Math.abs(dx) <= diameter / 2 && Math.abs(dy) <= diameter / 2) this.solid[index] = 1;
          }
          this.temp[index] = 0;
          if (this.solid[index] && y > 0 && y < ny - 1) this.temp[index] = 1;
        }
      }
    }

    initialize() {
      for (let i = 0; i < this.size; i++) {
        this.rho[i] = 1;
        const x = i % this.nx;
        const u = x < this.nx - 1 ? this.uIn : 0;
        this.ux[i] = this.solid[i] ? 0 : u;
        for (const k of q) this.f[k][i] = this.equilibrium(k, 1, this.ux[i], 0);
      }
    }

    equilibrium(k, rho, u, v) {
      const dot = cx[k] * u + cy[k] * v;
      return weights[k] * rho * (1 + 3 * dot + 4.5 * dot * dot - 1.5 * (u * u + v * v));
    }

    iterate() {
      const { nx, ny, size, f, next, solid, ux, uy, rho, tau } = this;
      for (let i = 0; i < size; i++) {
        if (solid[i]) continue;
        let density = 0, mx = 0, my = 0;
        for (const k of q) {
          const value = f[k][i];
          density += value;
          mx += value * cx[k];
          my += value * cy[k];
        }
        density = Math.max(.5, density);
        const u = mx / density, v = my / density;
        rho[i] = density; ux[i] = u; uy[i] = v;
        for (const k of q) {
          const equilibrium = this.equilibrium(k, density, u, v);
          next[k][i] = f[k][i] - (f[k][i] - equilibrium) / tau;
        }
      }
      for (let y = 1; y < ny - 1; y++) {
        for (let x = 1; x < nx - 1; x++) {
          const i = y * nx + x;
          if (solid[i]) continue;
          for (const k of q) {
            const sx = x - cx[k], sy = y - cy[k], source = sy * nx + sx;
            f[k][i] = solid[source] ? next[opposite[k]][i] : next[k][source];
          }
        }
      }
      for (let y = 1; y < ny - 1; y++) {
        const left = y * nx, right = left + nx - 1, leftRho = (f[0][left] + f[2][left] + f[4][left] + 2 * (f[3][left] + f[6][left] + f[7][left])) / (1 - this.uIn);
        f[1][left] = f[3][left] + 2 / 3 * leftRho * this.uIn;
        f[5][left] = f[7][left] + .5 * (f[4][left] - f[2][left]) + leftRho * this.uIn / 6;
        f[8][left] = f[6][left] + .5 * (f[2][left] - f[4][left]) + leftRho * this.uIn / 6;
        for (const k of q) f[k][right] = f[k][right - 1];
      }
      this.updateTemperature();
    }

    updateTemperature() {
      const { nx, ny, temp, tempNext, ux, uy, solid } = this;
      const alpha = Math.min(.12, this.diffusivity);
      for (let y = 1; y < ny - 1; y++) {
        for (let x = 1; x < nx - 1; x++) {
          const i = y * nx + x;
          if (solid[i]) { tempNext[i] = 1; continue; }
          const t = temp[i];
          const dx = ux[i] >= 0 ? t - temp[i - 1] : temp[i + 1] - t;
          const dy = uy[i] >= 0 ? t - temp[i - nx] : temp[i + nx] - t;
          const laplace = temp[i - 1] + temp[i + 1] + temp[i - nx] + temp[i + nx] - 4 * t;
          tempNext[i] = Math.max(0, Math.min(1, t - ux[i] * dx - uy[i] * dy + alpha * laplace));
        }
      }
      for (let y = 1; y < ny - 1; y++) {
        tempNext[y * nx] = 0;
        tempNext[y * nx + nx - 1] = tempNext[y * nx + nx - 2];
      }
      for (let x = 0; x < nx; x++) { tempNext[x] = tempNext[nx + x]; tempNext[(ny - 1) * nx + x] = tempNext[(ny - 2) * nx + x]; }
      this.temp = tempNext.slice();
      this.tempNext = temp;
    }

    solve(iterations, onProgress) {
      let step = 0;
      const batch = () => {
        const stop = Math.min(iterations, step + 18);
        for (; step < stop; step++) this.iterate();
        onProgress(step / iterations);
        if (step < iterations) requestAnimationFrame(batch);
        else onProgress(1, this.measure());
      };
      requestAnimationFrame(batch);
    }

    measure() {
      const { nx, ny, rho, ux, temp, solid } = this;
      let inlet = 0, outlet = 0, outletTemp = 0, outletCount = 0;
      const probe = Math.round(nx * .12);
      let meanU = 0, meanTemp = 0;
      for (let y = 1; y < ny - 1; y++) {
        const a = y * nx + probe, b = y * nx + nx - 2;
        inlet += rho[a];
        outlet += rho[b];
        meanU += ux[a];
        meanTemp += temp[b];
        outletTemp += temp[b];
        outletCount++;
      }
      const pressureDrop = Math.max(0, (inlet - outlet) / (ny - 2) / 3);
      return { pressureDrop, outletTemp: outletTemp / outletCount, meanU: meanU / (ny - 2) };
    }
  }

  function colorAt(value, max, field) {
    let t = Math.max(0, Math.min(1, value / (max || 1)));
    const stops = field === "temperature"
      ? [[225, 238, 225], [141, 191, 165], [71, 133, 119], [235, 174, 106], [213, 100, 72]]
      : [[224, 235, 226], [165, 203, 192], [91, 158, 157], [53, 117, 130], [228, 163, 105], [211, 105, 75]];
    if (field === "pressure") {
      const colors = [[89, 125, 153], [171, 201, 202], [242, 237, 215], [220, 159, 116], [151, 84, 76]];
      return blendStops(t, colors);
    }
    return blendStops(t, stops);
  }

  function blendStops(t, stops) {
    const pos = t * (stops.length - 1), index = Math.min(stops.length - 2, Math.floor(pos)), local = pos - index;
    const rgb = stops[index].map((v, i) => Math.round(v + (stops[index + 1][i] - v) * local));
    return `rgb(${rgb.join(",")})`;
  }

  function resizeCanvas() {
    const rect = canvas.getBoundingClientRect();
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = Math.max(1, Math.round(rect.width * dpr));
    canvas.height = Math.max(1, Math.round(rect.height * dpr));
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    render();
    drawChart(lastResults?.convergence || []);
  }

  function render() {
    if (!solver) { drawPlaceholder(); return; }
    const { nx, ny, ux, uy, rho, temp, solid } = solver;
    const width = canvas.clientWidth, height = canvas.clientHeight;
    const pad = { left: 32, right: 14, top: 12, bottom: 33 };
    const scale = Math.min((width - pad.left - pad.right) / nx, (height - pad.top - pad.bottom) / ny);
    const dw = scale * nx, dh = scale * ny, x0 = pad.left + (width - pad.left - pad.right - dw) / 2, y0 = pad.top + (height - pad.top - pad.bottom - dh) / 2;
    ctx.clearRect(0, 0, width, height);
    ctx.fillStyle = "#edf0e9"; ctx.fillRect(0, 0, width, height);
    const field = state.field;
    let min = Infinity, max = -Infinity;
    for (let i = 0; i < solver.size; i++) {
      if (solid[i]) continue;
      const value = field === "velocity" ? Math.hypot(ux[i], uy[i]) : field === "pressure" ? rho[i] - 1 : temp[i];
      min = Math.min(min, value); max = Math.max(max, value);
    }
    const span = max - min || 1;
    for (let y = 0; y < ny; y++) {
      for (let x = 0; x < nx; x++) {
        const i = y * nx + x, px = x0 + x * scale, py = y0 + (ny - 1 - y) * scale;
        if (solid[i]) { ctx.fillStyle = y === 0 || y === ny - 1 ? "#e5e6df" : "#374b40"; ctx.fillRect(px, py, scale + .35, scale + .35); continue; }
        const value = field === "velocity" ? Math.hypot(ux[i], uy[i]) : field === "pressure" ? rho[i] - 1 : temp[i];
        const normalized = (value - min) / span;
        ctx.fillStyle = colorAt(normalized, 1, field);
        ctx.fillRect(px, py, scale + .35, scale + .35);
        if (state.mesh && scale > 3) { ctx.strokeStyle = "#ffffff28"; ctx.lineWidth = .4; ctx.strokeRect(px, py, scale, scale); }
      }
    }
    if (state.streamlines) drawStreamlines(x0, y0, scale, nx, ny, ux, uy, solid);
    drawObstacle(x0, y0, scale, nx, ny, solver);
  }

  function drawObstacle(x0, y0, scale, nx, ny, current) {
    if (state.geometry === "none") return;
    const d = current.diameter * scale, cxp = x0 + current.centerX * scale + scale / 2, cyp = y0 + (ny - 1 - current.centerY) * scale + scale / 2;
    ctx.save(); ctx.fillStyle = "#34493e"; ctx.strokeStyle = "#f6f1e4"; ctx.lineWidth = 1.3;
    if (state.geometry === "cylinder") { ctx.beginPath(); ctx.arc(cxp, cyp, d / 2, 0, Math.PI * 2); ctx.fill(); ctx.stroke(); }
    else { ctx.beginPath(); ctx.roundRect(cxp - d / 2, cyp - d / 2, d, d, 2); ctx.fill(); ctx.stroke(); }
    ctx.restore();
  }

  function drawStreamlines(x0, y0, scale, nx, ny, ux, uy, solid) {
    const width = canvas.clientWidth, seeds = 17;
    ctx.save(); ctx.lineWidth = 1; ctx.strokeStyle = "rgba(255,255,247,.76)";
    for (let s = 0; s < seeds; s++) {
      let x = 2, y = 2 + s * (ny - 4) / (seeds - 1), points = 0;
      ctx.beginPath();
      while (x < nx - 2 && points < nx * 2) {
        const ix = Math.max(0, Math.min(nx - 1, Math.floor(x))), iy = Math.max(0, Math.min(ny - 1, Math.floor(y))), i = iy * nx + ix;
        if (solid[i]) break;
        const px = x0 + x * scale, py = y0 + (ny - 1 - y) * scale;
        if (points === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
        const u = ux[i], v = uy[i], speed = Math.hypot(u, v);
        if (speed < .002) break;
        const ds = Math.min(.8, .45 / (speed / .06 + .2));
        x += (u / speed) * ds; y += (v / speed) * ds;
        if (y < 1 || y > ny - 2) break;
        points++;
      }
      if (points > 1) ctx.stroke();
    }
    ctx.restore();
  }

  function drawPlaceholder() {
    const width = canvas.clientWidth, height = canvas.clientHeight;
    ctx.clearRect(0, 0, width, height); ctx.fillStyle = "#edf0e9"; ctx.fillRect(0, 0, width, height);
    ctx.fillStyle = "#e4e9e1"; ctx.fillRect(32, height * .21, width - 46, height * .58);
    ctx.fillStyle = "#dfe2d9"; ctx.fillRect(32, height * .21, width - 46, 1); ctx.fillRect(32, height * .79, width - 46, 1);
    if (state.geometry !== "none") {
      const x = width * .39, y = height * .5, r = Math.min(height * state.blockage * .29, 42);
      ctx.fillStyle = "#34493e";
      if (state.geometry === "cylinder") { ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); ctx.fill(); }
      else ctx.fillRect(x - r, y - r, 2 * r, 2 * r);
    }
    ctx.strokeStyle = "#95ab98"; ctx.lineWidth = 1;
    for (let i = 0; i < 7; i++) {
      const y = height * (.27 + i * .075);
      ctx.beginPath(); ctx.moveTo(40, y);
      if (state.geometry === "none") ctx.lineTo(width - 18, y);
      else {
        const ox = width * state.obstacleX, oy = height * .5, r = Math.min(height * state.blockage * .29, 42);
        ctx.bezierCurveTo(width * .23, y, ox - r * 1.5, y, ox - r, oy + (y - oy) * .55);
        ctx.bezierCurveTo(ox + r, oy + (y - oy) * .55, width * .58, y, width - 18, y);
      }
      ctx.stroke();
    }
  }

  function drawChart(points) {
    const width = chart.clientWidth, height = chart.clientHeight, dpr = Math.min(window.devicePixelRatio || 1, 2);
    if (!width || !height) return;
    chart.width = Math.round(width * dpr); chart.height = Math.round(height * dpr);
    chartCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
    chartCtx.clearRect(0, 0, width, height);
    const pad = { l: 28, r: 8, t: 8, b: 19 }, plotW = width - pad.l - pad.r, plotH = height - pad.t - pad.b;
    chartCtx.strokeStyle = "#e8ebe4"; chartCtx.lineWidth = 1;
    for (let i = 0; i < 4; i++) { const y = pad.t + i * plotH / 3; chartCtx.beginPath(); chartCtx.moveTo(pad.l, y); chartCtx.lineTo(width - pad.r, y); chartCtx.stroke(); }
    chartCtx.fillStyle = "#919a91"; chartCtx.font = "8px monospace"; chartCtx.fillText("Cp", 2, 13);
    if (!points.length) {
      chartCtx.setLineDash([3, 4]); chartCtx.strokeStyle = "#a9b5aa"; chartCtx.beginPath(); chartCtx.moveTo(pad.l, pad.t + plotH * .45); chartCtx.lineTo(width - pad.r, pad.t + plotH * .45); chartCtx.stroke(); chartCtx.setLineDash([]);
      ["COARSE", "MEDIUM", "FINE"].forEach((label, i) => chartCtx.fillText(label, pad.l + i * plotW / 2.5, height - 4));
      return;
    }
    const values = points.map(p => p.value), low = Math.min(...values) * .85, high = Math.max(...values) * 1.15 || 1;
    const coords = points.map((point, i) => [pad.l + i * plotW / Math.max(1, points.length - 1), pad.t + plotH * (1 - (point.value - low) / (high - low || 1))]);
    chartCtx.strokeStyle = "#5d8769"; chartCtx.lineWidth = 2; chartCtx.beginPath();
    coords.forEach(([x, y], i) => i ? chartCtx.lineTo(x, y) : chartCtx.moveTo(x, y)); chartCtx.stroke();
    coords.forEach(([x, y], i) => { chartCtx.fillStyle = "#fffefa"; chartCtx.strokeStyle = "#5d8769"; chartCtx.lineWidth = 2; chartCtx.beginPath(); chartCtx.arc(x, y, 3.2, 0, Math.PI * 2); chartCtx.fill(); chartCtx.stroke(); chartCtx.fillStyle = "#929990"; chartCtx.font = "8px monospace"; chartCtx.fillText(`${points[i].nx}×${points[i].ny}`, Math.min(width - 35, x - 13), height - 4); });
  }

  function options(nx = 96, ny = 48, scaleProperties = false) {
    const scale = scaleProperties ? ny / 48 : 1;
    return { nx, ny, velocity: state.velocity * scale, viscosity: state.viscosity * scale, diffusivity: state.diffusivity * scale, geometry: state.geometry, blockage: state.blockage, obstacleX: state.obstacleX };
  }

  function invalidateSolution() {
    state.solved = false;
    solver = null;
    lastResults = null;
    ["re-result", "pressure-result", "temp-result"].forEach(id => { byId(id).textContent = "—"; });
    const nx = Math.round(96 * state.length / 4);
    byId("grid-result").textContent = `${nx} × 48`;
    byId("grid-label").textContent = `${nx} × 48 CELLS`;
    byId("canvas-status").textContent = "SOLVER IDLE";
    byId("convergence-status").textContent = "NO STUDY RUN YET";
    drawChart([]);
    updateReadouts();
  }

  function updateReadouts() {
    byId("length-value").textContent = Number(state.length).toFixed(1);
    byId("blockage-value").textContent = Number(state.blockage).toFixed(2);
    byId("velocity-value").textContent = Number(state.velocity).toFixed(2);
    byId("viscosity-value").textContent = Number(state.viscosity).toFixed(3).replace(/0$/, "");
    byId("diffusivity-value").textContent = Number(state.diffusivity).toFixed(3);
    byId("field-badge").lastChild.textContent = ` ${fields[state.field].label}`;
    byId("canvas-hint").style.display = state.solved ? "none" : "block";
    if (solver) render();
    else drawPlaceholder();
  }

  function setProgress(value) { byId("progress-bar").style.width = `${Math.max(0, Math.min(100, value * 100))}%`; }

  function runSimulation() {
    if (isRunning) return;
    isRunning = true; byId("solve-button").disabled = true; byId("run-state").textContent = "SOLVING FLOW FIELD"; byId("canvas-status").textContent = "SOLVER RUNNING";
    byId("solve-time").textContent = "ITERATING"; setProgress(0);
    const nx = Math.round(96 * state.length / 4);
    solver = new LatticeSolver(nx, 48, options(nx, 48));
    byId("grid-label").textContent = `${nx} × 48 CELLS`;
    solver.solve(1800, (progress, result) => {
      setProgress(progress);
      if (progress === 1 && result) {
        isRunning = false; state.solved = true; lastResults = result;
        byId("solve-button").disabled = false; byId("run-state").textContent = "SOLUTION CONVERGED"; byId("canvas-status").textContent = "FIELD SOLVED";
        byId("solve-time").textContent = "1800 ITERATIONS";
        const re = state.velocity * solver.diameter / state.viscosity;
        byId("re-result").textContent = re.toFixed(1);
        byId("pressure-result").textContent = result.pressureDrop.toFixed(4);
        byId("temp-result").textContent = result.outletTemp.toFixed(3);
        byId("grid-result").textContent = `${nx} × 48`;
        render();
      }
    });
  }

  async function runConvergence() {
    if (isRunning) return;
    isRunning = true;
    const button = byId("convergence-button");
    button.disabled = true; button.textContent = "Solving…";
    const grids = [[40, 20], [64, 32], [96, 48]].map(([nx, ny]) => [Math.round(nx * state.length / 4), ny]), results = [];
    byId("convergence-status").textContent = "RUNNING THREE MESHES";
    for (const [nx, ny] of grids) {
      byId("convergence-status").textContent = `SOLVING ${nx} × ${ny}`;
      const candidate = new LatticeSolver(nx, ny, options(nx, ny, true));
      const result = await new Promise(resolve => candidate.solve(350, (progress, value) => { if (value) resolve(value); }));
      results.push({ nx, ny, value: 2 * result.pressureDrop / (candidate.uIn * candidate.uIn) });
      drawChart(results);
    }
    const difference = results.length > 1 ? Math.abs(results[2].value - results[1].value) / Math.max(Math.abs(results[2].value), 1e-8) * 100 : 0;
    byId("convergence-status").textContent = `FINE GRID CHANGE: ${difference.toFixed(1)}%`;
    lastResults = { ...(lastResults || {}), convergence: results };
    isRunning = false; button.disabled = false; button.innerHTML = 'Run study <span>↗</span>';
    drawChart(results);
  }

  document.querySelectorAll("input[type=range]").forEach(input => {
    input.addEventListener("input", () => {
      state[input.id] = Number(input.value);
      input.style.background = `linear-gradient(90deg,#7d9d87 ${((input.value - input.min) / (input.max - input.min)) * 100}%,#e6e9e2 50%)`;
      invalidateSolution();
    });
  });
  byId("geometry").addEventListener("change", event => { state.geometry = event.target.value; invalidateSolution(); });
  document.querySelectorAll(".field-tab").forEach(button => button.addEventListener("click", () => {
    document.querySelectorAll(".field-tab").forEach(tab => tab.classList.toggle("active", tab === button));
    state.field = button.dataset.field; updateReadouts();
  }));
  byId("streamlines").addEventListener("change", event => { state.streamlines = event.target.checked; render(); });
  byId("mesh-overlay").addEventListener("change", event => { state.mesh = event.target.checked; render(); });
  byId("solve-button").addEventListener("click", runSimulation);
  byId("convergence-button").addEventListener("click", runConvergence);
  byId("reset-button").addEventListener("click", () => { if (isRunning) return; window.location.reload(); });
  byId("fullscreen-button").addEventListener("click", () => {
    const viewport = byId("viewport"); viewport.classList.toggle("expanded"); byId("fullscreen-button").textContent = viewport.classList.contains("expanded") ? "×" : "⛶"; resizeCanvas();
  });
  canvas.addEventListener("pointerdown", event => {
    if (!solver || state.geometry === "none") return;
    const rect = canvas.getBoundingClientRect(), nx = solver.nx;
    const scale = Math.min((rect.width - 46) / nx, (rect.height - 45) / solver.ny);
    const x0 = 32 + (rect.width - 46 - scale * nx) / 2;
    const cxp = x0 + solver.centerX * scale;
    if (Math.abs(event.offsetX - cxp) < solver.diameter * scale && Math.abs(event.offsetY - rect.height / 2) < solver.diameter * scale) { dragging = true; canvas.setPointerCapture(event.pointerId); }
  });
  canvas.addEventListener("pointermove", event => {
    if (!dragging || !solver) return;
    const rect = canvas.getBoundingClientRect(), scale = Math.min((rect.width - 46) / solver.nx, (rect.height - 45) / solver.ny);
    const x0 = 32 + (rect.width - 46 - scale * solver.nx) / 2;
    state.obstacleX = Math.max(.18, Math.min(.76, (event.offsetX - x0) / (scale * solver.nx)));
    drawPlaceholder();
  });
  canvas.addEventListener("pointerup", () => { if (dragging) { dragging = false; invalidateSolution(); } });
  window.addEventListener("resize", resizeCanvas);
  updateReadouts(); resizeCanvas();
})();
