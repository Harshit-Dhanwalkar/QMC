import fs from "node:fs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_haldane.wasm", import.meta.url),
);

// complex eigensolver carries stdio error paths, so module imports three WASI calls; none runs on valid input.
const wasi = { fd_write: () => 8, fd_seek: () => 8, fd_close: () => 8 };
const { instance } = await WebAssembly.instantiate(bytes, {
  wasi_snapshot_preview1: wasi,
});
const w = instance.exports;
if (w._initialize) w._initialize();

const f = (x) => x.toFixed(5);
// memory can grow, so make views after calls that touch it
const E = () =>
  new Float64Array(w.memory.buffer, w.qmc_hd_energy_buffer(), 121 * 28);
const M = () =>
  new Float64Array(w.memory.buffer, w.qmc_hd_mean_buffer(), 121 * 28);

w.qmc_hd_set(0.1, Math.PI / 2, 0);
w.qmc_hd_compute(25);

let e = E();
let m = M();
console.log(
  `chern=${w.qmc_hd_chern()} limit=${f(w.qmc_hd_limit())} gapK=${f(w.qmc_hd_valley_gap(1))} gapKp=${f(w.qmc_hd_valley_gap(-1))} vmax=${f(w.qmc_hd_vmax())} cmin=${f(w.qmc_hd_cmin())}`,
);
console.log(
  `k=pi E13=${f(Math.abs(e[12 * 28 + 13]))} E14=${f(Math.abs(e[12 * 28 + 14]))}`,
);
console.log(
  `k=3pi/4 E13=${f(e[9 * 28 + 13])} E14=${f(e[9 * 28 + 14])} m13=${m[9 * 28 + 13].toFixed(4)} m14=${m[9 * 28 + 14].toFixed(4)}`,
);
console.log(`k=0 E0=${f(e[0])} E27=${f(e[27])}`);

w.qmc_hd_profile(9, 25, 13);
const p = new Float64Array(w.memory.buffer, w.qmc_hd_prob_buffer(), 14);
console.log(`profile p0=${f(p[0])} p1=${f(p[1])} p13=${p[13].toFixed(6)}`);
w.qmc_hd_set(0.1, -Math.PI / 2, 1);
console.log(`trivial chern=${w.qmc_hd_chern()}`);
if (w.qmc_hd_set(0.5, 0, 0) !== -1 || w.qmc_hd_compute(500) !== -1) {
  console.error("bad args accepted");
  process.exit(1);
}
