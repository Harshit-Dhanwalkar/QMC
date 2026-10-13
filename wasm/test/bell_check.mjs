import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_bell.wasm", import.meta.url),
);

let memory = null;
const { instance } = await WebAssembly.instantiate(
  bytes,
  makeWasiImports(() => memory),
);
memory = instance.exports.memory;
const w = instance.exports;
if (w._initialize) w._initialize();

const f = (x) => x.toFixed(5);
const pi = Math.PI;

// memory can grow, so make views after calls that touch it
const Q = () =>
  new Float64Array(w.memory.buffer, w.qmc_bl_quantum_buffer(), 361);
const L = () => new Float64Array(w.memory.buffer, w.qmc_bl_local_buffer(), 361);
w.qmc_bl_set(0, pi / 2, pi / 4, (3 * pi) / 4, 1.0);
console.log(
  `S=${f(w.qmc_bl_chsh())} local=${f(w.qmc_bl_local_chsh())} max=${f(w.qmc_bl_max())} e0=${f(w.qmc_bl_e(0))} e1=${f(w.qmc_bl_e(1))} el0=${f(w.qmc_bl_e_local(0))}`,
);

w.qmc_bl_set(0, pi / 2, pi / 4, (3 * pi) / 4, 0.6);
console.log(`V=0.6 S=${f(w.qmc_bl_chsh())} max=${f(w.qmc_bl_max())}`);
w.qmc_bl_set(0.1, 1.2, 0.7, 2.0, 0.9);
console.log(`odd S=${f(w.qmc_bl_chsh())} local=${f(w.qmc_bl_local_chsh())}`);
w.qmc_bl_curve(181);

const q = Q();
const l = L();
console.log(
  `curve q0=${f(q[0])} q40=${f(q[40])} q90=${f(q[90])} l40=${f(l[40])} l60=${f(l[60])}`,
);

w.qmc_bl_set(0, pi / 2, pi / 4, (3 * pi) / 4, 1.0);
const sEst = w.qmc_bl_sample(5000, 42);
console.log(
  `sample S=${f(sEst)} e0=${f(w.qmc_bl_est(0))} sigma=${f(w.qmc_bl_sigma())}`,
);

if (
  w.qmc_bl_set(0, 0, 0, 0, 1.5) !== -1 ||
  w.qmc_bl_curve(1) !== -1 ||
  !Number.isNaN(w.qmc_bl_e(7))
) {
  console.log("bad args accepted");
  process.exit(1);
}
console.log("bad args rejected");
