import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_tfim.wasm", import.meta.url),
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
// memory can grow, so make the views after the calls that touch it
const Z = () => new Float64Array(w.memory.buffer, w.qmc_tq_zz_buffer(), 41);
w.qmc_tq_set(0.3, 2.0);
w.qmc_tq_eval(0);
let z = Z();
console.log(
  `t=0 mx=${f(w.qmc_tq_last_mx())} rate=${f(w.qmc_tq_last_rate())} zz1=${f(z[1])} zz10=${f(z[10])}`,
);

w.qmc_tq_eval(0.8);
z = Z();
console.log(
  `t=0.8 mx=${f(w.qmc_tq_last_mx())} rate=${f(w.qmc_tq_last_rate())} zz1=${f(z[1])} zz2=${f(z[2])} zz5=${f(z[5])}`,
);
console.log(
  `tcrit ${f(w.qmc_tq_tcrit(0))} ${f(w.qmc_tq_tcrit(1))} front ${w.qmc_tq_front_speed().toFixed(3)}`,
);

w.qmc_tq_curve(3.0, 61);
const r = new Float64Array(w.memory.buffer, w.qmc_tq_rate_buffer(), 400),
  m = new Float64Array(w.memory.buffer, w.qmc_tq_mx_buffer(), 400);
console.log(
  `curve r[10]=${f(r[10])} r[30]=${f(r[30])} r[60]=${f(r[60])} m[10]=${f(m[10])} m[60]=${f(m[60])}`,
);

w.qmc_tq_set(3.0, 0.5);
w.qmc_tq_eval(3.0);
z = Z();
console.log(
  `para->ordered zz3=${f(z[3])} zz4=${f(z[4])} zz30_ok=${Math.abs(z[30]) < 1e-6 ? 1 : 0} tcrit_nan=${Number.isNaN(w.qmc_tq_tcrit(0)) ? 0 : 1}`,
);

const d = new Float64Array(w.memory.buffer, w.qmc_tq_density_buffer(), 128);
console.log(`density d[0]=${f(d[0])} d[60]=${f(d[60])}`);
if (
  w.qmc_tq_set(3.5, 1.0) !== -1 ||
  w.qmc_tq_curve(1.0, 1000) !== -1 ||
  w.qmc_tq_eval(-1) !== -1
) {
  console.error("bad args accepted");
  process.exit(1);
}
