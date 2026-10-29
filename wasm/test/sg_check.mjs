import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_sg.wasm", import.meta.url),
);

let memory = null;
const { instance } = await WebAssembly.instantiate(
  bytes,
  makeWasiImports(() => memory),
);
memory = instance.exports.memory;
const w = instance.exports;
if (w._initialize) w._initialize();

const pr = (label, x) => `${label}=${(Math.abs(x) < 5e-6 ? 0 : x).toFixed(5)} `;

// memory can grow, so make views after calls that touch it
const buf = (fn, n) => new Float64Array(w.memory.buffer, fn(), n);
w.qmc_sg_set(0, 0, 0, 0, 1);
w.qmc_sg_set(1, Math.PI / 2, 0, 0, 1);
w.qmc_sg_set(2, 0, 0, 0, 1);
w.qmc_sg_set(3, 0.8, 0.4, 2, 0.35);
w.qmc_sg_run(0, 0, 1, 3);

let inn = buf(w.qmc_sg_intensity_buffer, 5);
let r = buf(w.qmc_sg_r_buffer, 15);
let s = "zxz ";
for (let k = 0; k <= 3; k++) {
  s += pr("i", inn[k]);
}
console.log(s);

w.qmc_sg_run(1.0, 0.5, 0.7, 4);
inn = buf(w.qmc_sg_intensity_buffer, 5);
r = buf(w.qmc_sg_r_buffer, 15);
s = "mixed ";
for (let k = 0; k <= 4; k++) {
  s += pr("i", inn[k]) + pr("rz", r[3 * k + 2]);
}
console.log(s);

const ov = w.qmc_sg_screen(0.4, 3.0, 0.7, 101, 6.0);
const q = buf(w.qmc_sg_quantum_buffer, 101);
const c = buf(w.qmc_sg_classical_buffer, 101);
console.log(
  "screen " +
    pr("ov", ov) +
    pr("q20", q[20]) +
    pr("q50", q[50]) +
    pr("q85", q[85]) +
    pr("c50", c[50]) +
    pr("c85", c[85]),
);

const ups = w.qmc_sg_shoot(500, 0.3, 3.0, 0.7, 7);
const h = buf(w.qmc_sg_hits_buffer, 1000);
console.log(
  `shots ups=${ups} ` + pr("z0", h[0]) + pr("z1", h[2]) + pr("z499", h[998]),
);

if (
  w.qmc_sg_set(4, 0, 0, 0, 1) !== -1 ||
  w.qmc_sg_set(0, 0, 0, 3, 1) !== -1 ||
  w.qmc_sg_set(0, 0, 0, 0, 2) !== -1 ||
  w.qmc_sg_run(0, 0, 1.5, 1) !== -1 ||
  w.qmc_sg_run(0, 0, 1, 5) !== -1 ||
  w.qmc_sg_screen(0.5, 1, 0.0, 10, 5) !== -1 ||
  w.qmc_sg_shoot(99999, 0.5, 1, 1, 1) !== -1
) {
  console.log("bad args accepted");
  process.exit(1);
}
console.log("bad args rejected");
