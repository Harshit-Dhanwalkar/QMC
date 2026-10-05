// Physics check: run the playground's double-slit scenario and compare the
// detector fringe spacing with the textbook estimate  dy = lambda * L / d.
import fs from "node:fs";
const N = 128, BOX = 40, DT = 0.02, DX = BOX / N;
const { instance } = await WebAssembly.instantiate(
  fs.readFileSync(new URL("../../docs/src/playground/qmc.wasm", import.meta.url)), {});
const w = instance.exports; if (w._initialize) w._initialize();
const k = 3.0, sep = 6.0, slitW = 1.6, wallX = 2.0, screenX = 14.0;
w.qmc_init(N, BOX, DT);
w.qmc_build_slits(wallX, 1.0, 2, sep, slitW, 60.0);
w.qmc_wavepacket(-12, 0, k, 0, 2.0);
const SX = Math.round(N / 2 + screenX / DX);
const acc = new Float64Array(N);
for (let f = 0; f < 200; f++) {          // 200 frames x 4 steps = t = 16
  w.qmc_step(4);
  const d = new Float32Array(w.memory.buffer, w.qmc_density(), N * N);
  for (let iy = 0; iy < N; iy++) acc[iy] += d[SX * N + iy];
}
// local maxima of the accumulated pattern (smoothed)
const sm = acc.map((_, i) => { let s = 0, c = 0; for (let j = -1; j <= 1; j++) { const a = acc[i + j]; if (a !== undefined) { s += a; c++; } } return s / c; });
const max = Math.max(...sm), peaks = [];
for (let i = 2; i < N - 2; i++) if (sm[i] > sm[i - 1] && sm[i] >= sm[i + 1] && sm[i] > 0.12 * max) peaks.push(i);
const ys = peaks.map((i) => (i - N / 2) * DX);
console.log("peak y positions:", ys.map((y) => y.toFixed(2)).join(", "));
const gaps = ys.slice(1).map((y, i) => y - ys[i]);
console.log("gaps:", gaps.map((g) => g.toFixed(2)).join(", "));
const lambda = 2 * Math.PI / k, L = screenX - wallX, pred = lambda * L / sep;
// Compare only the two fringe gaps around the optical axis: far from it the
// small-angle formula (and the far-field assumption) stops being accurate.
const centre = ys.findIndex((y) => Math.abs(y) < DX);
const cgaps = [ys[centre + 1] - ys[centre], ys[centre] - ys[centre - 1]];
const meas = (cgaps[0] + cgaps[1]) / 2, err = Math.abs(meas - pred) / pred;
console.log(`predicted fringe spacing lambda*L/d = ${pred.toFixed(2)}, ` +
  `measured (central gaps) = ${meas.toFixed(2)}, error ${(100 * err).toFixed(1)}%`);
console.log("norm at end:", w.qmc_norm().toFixed(4));
if (!(err < 0.1)) { console.error("FAIL: central fringe spacing off by more than 10%"); process.exit(1); }
