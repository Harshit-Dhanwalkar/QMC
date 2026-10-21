import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_ring.wasm", import.meta.url),
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
// memory can grow, so make views after calls that touch it
const buf = (fn, n) => new Float64Array(w.memory.buffer, fn(), n);
w.qmc_ar_spectrum(12);
const sp = buf(w.qmc_ar_spec_buffer, 32 * 201);
console.log(
  `spec l0 i0=${f(sp[0])} i50=${f(sp[50])} i100=${f(sp[100])} l11 i75=${f(sp[11 * 201 + 75])}`,
);

w.qmc_ar_current(11, 5);
const cu = buf(w.qmc_ar_cur_buffer, 201);
console.log(
  `cur i20=${f(cu[20])} i60=${f(cu[60])} i120=${f(cu[120])} i140=${f(cu[140])}`,
);

const e = w.qmc_ar_levels(12, 6, 0.21);
const lv = buf(w.qmc_ar_levels_buffer, 32);
console.log(
  `levels E6=${f(e)} l0=${f(lv[0])} l5=${f(lv[5])} l11=${f(lv[11])} I=${f(w.qmc_ar_current_at(12, 6, 0.21))}`,
);

const t = w.qmc_ar_point(16, 8, 0.15, 0.6);
console.log(`point T=${f(t)} R=${f(w.qmc_ar_last_reflection())}`);

w.qmc_ar_trans_flux(16, 8, 0.4);
const tf = buf(w.qmc_ar_tf_buffer, 201);
console.log(
  `tflux t0=${f(tf[0])} t50=${f(tf[50])} t100=${f(tf[100])} t200=${f(tf[200])}`,
);

w.qmc_ar_trans_energy(16, 5, 0.3);
const te = buf(w.qmc_ar_te_buffer, 201);
console.log(
  `tenergy t0=${f(te[0])} t60=${f(te[60])} t100=${f(te[100])} t180=${f(te[180])}`,
);

w.qmc_ar_map(10, 5);
const mp = buf(w.qmc_ar_map_buffer, 4096);
console.log(
  `map ${f(mp[0])} ${f(mp[20 * 64 + 10])} ${f(mp[32 * 64 + 32])} ${f(mp[63 * 64 + 40])} ${f(mp[40 * 64 + 63])}`,
);
if (
  w.qmc_ar_spectrum(3) !== -1 ||
  w.qmc_ar_current(8, 9) !== -1 ||
  w.qmc_ar_trans_flux(8, 0, 0) !== -1 ||
  w.qmc_ar_trans_flux(8, 4, 2.5) !== -1 ||
  w.qmc_ar_map(40, 3) !== -1 ||
  !Number.isNaN(w.qmc_ar_point(8, 8, 0, 0)) ||
  !Number.isNaN(w.qmc_ar_levels(8, 9, 0))
) {
  console.log("bad args accepted");
  process.exit(1);
}
console.log("bad args rejected");
