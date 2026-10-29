import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_zeeman.wasm", import.meta.url),
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
w.qmc_zm_curves(1.0, 5.0, 101);
const c = buf(w.qmc_zm_curves_buffer, 808);
const mj = buf(w.qmc_zm_mj_buffer, 8);

let s = "mj ";

for (let k = 0; k < 8; k++) {
  s += pr("m", mj[k]);
}
console.log(s);

for (let k = 0; k < 8; k++)
  console.log(
    `curve${k} ` +
      pr("b0", c[k * 101]) +
      pr("b1", c[k * 101 + 20]) +
      pr("b5", c[k * 101 + 100]),
  );

const n = w.qmc_zm_lines(1.0, 0.7);
const l = buf(w.qmc_zm_lines_buffer, 60);
console.log(`lines n=${n}`);

for (let j = 0; j < 12; j++)
  console.log(
    `row${j} ` +
      pr("e", l[j * 5]) +
      pr("i", l[j * 5 + 1]) +
      pr("q", l[j * 5 + 2]),
  );
console.log(
  "moment " +
    pr("m0", w.qmc_zm_moment(1, 0.7, 0)) +
    pr("m4", w.qmc_zm_moment(1, 0.7, 4)) +
    pr("m5", w.qmc_zm_moment(1, 0.7, 5)),
);

if (
  w.qmc_zm_curves(0, 5, 101) !== -1 ||
  w.qmc_zm_curves(1, 5, 1) !== -1 ||
  w.qmc_zm_curves(1, 5, 9999) !== -1 ||
  w.qmc_zm_lines(-1, 1) !== -1 ||
  !Number.isNaN(w.qmc_zm_moment(1, 1, 7))
) {
  console.log("bad args accepted");
  process.exit(1);
}
console.log("bad args rejected");
