import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_josephson.wasm", import.meta.url),
);

let memory = null;
const { instance } = await WebAssembly.instantiate(
  bytes,
  makeWasiImports(() => memory),
);
memory = instance.exports.memory;
const w = instance.exports;
if (w._initialize) w._initialize();

// memory can grow, so make views after calls that touch it
const buf = (fn, n) => new Float64Array(w.memory.buffer, fn(), n);
w.qmc_jj_set(3.0, 1.4, 0.0, 0.5);
w.qmc_jj_reset(0.2, 0.0);
w.qmc_jj_advance(10.0);
w.qmc_jj_advance(10.0);
w.qmc_jj_advance(5.0);
console.log(
  `phi=${f(w.qmc_jj_phi())} v=${f(w.qmc_jj_v())} t=${f(w.qmc_jj_time())} `,
);

w.qmc_jj_set(3.0, 0.0, 0.0, 0.5);
w.qmc_jj_iv(2.0);
let up = buf(w.qmc_jj_up_buffer, 161);
let dn = buf(w.qmc_jj_down_buffer, 161);
console.log(
  `hyst up40=${f(up[40])} up100=${f(up[100])} up160=${f(up[160])} dn40=${f(dn[40])} dn80=${f(dn[80])} dn160=${f(dn[160])} `,
);

w.qmc_jj_set(0.7, 0.0, 1.2, 0.8);
w.qmc_jj_iv(2.2);
up = buf(w.qmc_jj_up_buffer, 161);
dn = buf(w.qmc_jj_down_buffer, 161);
console.log(
  `shapiro up60=${f(up[60])} up90=${f(up[90])} up120=${f(up[120])} dn90=${f(dn[90])} `,
);

if (
  w.qmc_jj_set(0, 0, 0, 0) !== -1 ||
  w.qmc_jj_set(1, 9, 0, 0) !== -1 ||
  w.qmc_jj_iv(0) !== -1 ||
  w.qmc_jj_advance(-1) !== -1 ||
  w.qmc_jj_reset(NaN, 0) !== -1
) {
  console.log("bad args accepted");
  process.exit(1);
}
console.log("bad args rejected");
