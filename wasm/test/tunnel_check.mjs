import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_tunnel.wasm", import.meta.url),
);

let memory = null;
const { instance } = await WebAssembly.instantiate(
  bytes,
  makeWasiImports(() => memory),
);
memory = instance.exports.memory;
const w = instance.exports;
if (w._initialize) w._initialize();

const f = (x, d = 5) => x.toFixed(d);

// memory can grow, so make views after calls that touch it
const buf = (fn, n) => new Float64Array(w.memory.buffer, fn(), n);
w.qmc_dt_set(30.0, 0.5, 3.0, 0.5, 30.0, 0.0);
console.log(
  `length ${f(w.qmc_dt_length())} well1 ${f(w.qmc_dt_well_level(1))} well2 ${f(w.qmc_dt_well_level(2))}`,
);

w.qmc_dt_curves(20.0);
const te = buf(w.qmc_dt_te_buffer, 401);
const dw = buf(w.qmc_dt_dwell_buffer, 401);
console.log(
  `T i0=${f(te[0])} i20=${f(te[20])} i100=${f(te[100])} i300=${f(te[300])} i400=${f(te[400])}`,
);

console.log(`dwell i20=${f(dw[20])} i100=${f(dw[100])} i300=${f(dw[300])}`);
let t = w.qmc_dt_point(2.2);
console.log(
  `point T=${f(t)} R=${f(w.qmc_dt_last_reflection())} dwell=${f(w.qmc_dt_dwell(2.2))}`,
);

w.qmc_dt_profile(2.2, -2.0, 6.0);
const ps = buf(w.qmc_dt_psi_buffer, 401);
const pt = buf(w.qmc_dt_pot_buffer, 401);
console.log(
  `\\psi ${f(ps[0])} ${f(ps[100])} ${f(ps[200])} ${f(ps[300])} ${f(ps[400])} | pot ${f(pt[0])} ${f(pt[150])} ${f(pt[200])}`,
);

w.qmc_dt_set(20.0, 0.3, 2.0, 0.6, 35.0, -3.0);
t = w.qmc_dt_point(4.0);
console.log(
  `asym T=${f(t)} R=${f(w.qmc_dt_last_reflection())} length ${f(w.qmc_dt_length())}`,
);

if (
  w.qmc_dt_set(-1, 0.5, 3, 0.5, 30, 0) !== -1 ||
  w.qmc_dt_set(30, 0, 3, 0.5, 30, 0) !== -1 ||
  w.qmc_dt_curves(0) !== -1 ||
  !Number.isNaN(w.qmc_dt_point(-1)) ||
  w.qmc_dt_profile(1, 2, 1) !== -1
) {
  console.log("bad args accepted");
  process.exit(1);
}
console.log("bad args rejected");
