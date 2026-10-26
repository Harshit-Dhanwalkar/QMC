import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_doublewell.wasm", import.meta.url),
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

// memory can grow, so make the views after the calls that touch it
const buf = (fn, n) => new Float64Array(w.memory.buffer, fn(), n);
w.qmc_dw_set(8.0, 1.5);
let e = buf(w.qmc_dw_e_buffer, 8);
let v = buf(w.qmc_dw_v_buffer, 500);
let s = "levels ";

for (let k = 0; k < 6; k++) {
  s += pr("e", e[k]);
}
console.log(s);
console.log(
  "split " +
    pr("d", w.qmc_dw_splitting()) +
    pr("wkb", w.qmc_dw_wkb()) +
    pr("v0", v[0]) +
    pr("v250", v[250]),
);

const cap = w.qmc_dw_launch(-1.5, 0.35);
s = "launch " + pr("cap", cap);
const d = w.qmc_dw_splitting();
for (let i = 0; i < 4; i++) {
  const t = (i * 0.6 * Math.PI) / d / 3.0;
  const p = w.qmc_dw_evolve(t);
  s += pr("P", p) + pr("x", w.qmc_dw_mean_x());
}
console.log(s);

w.qmc_dw_evolve(1.0);
const dn = buf(w.qmc_dw_density_buffer, 500);
console.log(
  "dens " + pr("d100", dn[100]) + pr("d150", dn[150]) + pr("d350", dn[350]),
);

w.qmc_dw_split_curve(20.0);
const sp = buf(w.qmc_dw_split_buffer, 61);
const wk = buf(w.qmc_dw_wkb_buffer, 61);
console.log(
  "curve " +
    pr("s10", sp[10]) +
    pr("s30", sp[30]) +
    pr("s60", sp[60]) +
    pr("w30", wk[30]) +
    pr("w60", wk[60]),
);

if (
  w.qmc_dw_set(0, 1.5) !== -1 ||
  w.qmc_dw_set(8, 9) !== -1 ||
  w.qmc_dw_split_curve(1) !== -1 ||
  w.qmc_dw_launch(-1.5, 5) !== -1 ||
  w.qmc_dw_launch(9, 0.3) !== -1
) {
  console.log("bad args accepted");
  process.exit(1);
}
console.log("bad args rejected");
