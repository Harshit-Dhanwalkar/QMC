import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_zeno.wasm", import.meta.url),
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
console.log(
  `pulsed N=1 ${f(w.qmc_zn_pulsed(1))} N=4 ${f(w.qmc_zn_pulsed(4))} N=50 ${f(w.qmc_zn_pulsed(50))}`,
);

w.qmc_zn_pulsed_curves(8);
let t = buf(w.qmc_zn_time_buffer, 401);
let nc = buf(w.qmc_zn_n_buffer, 400);
console.log(
  `time c0=${f(t[0])} c100=${f(t[100])} c250=${f(t[250])} c400=${f(t[400])}`,
);
console.log(
  `ncurve n1=${f(nc[0])} n8=${f(nc[7])} n100=${f(nc[99])} n400=${f(nc[399])}`,
);

const end = w.qmc_zn_continuous(5.0);
const cc = buf(w.qmc_zn_cont_buffer, 401);
console.log(`continuous g5 end=${f(end)} c200=${f(cc[200])}`);
w.qmc_zn_rates(10.0);

const r = buf(w.qmc_zn_rate_buffer, 401);
console.log(
  `rates r0=${f(r[0])} r80=${f(r[80])} r400=${f(r[400])} peak ${f(w.qmc_zn_rate(2.0))}`,
);

const s = w.qmc_zn_sample(10, 5000, 7);
const al = buf(w.qmc_zn_alive_buffer, 401);
console.log(`sample final=${f(s)} a1=${f(al[1])} a5=${f(al[5])}`);
if (
  w.qmc_zn_pulsed_curves(0) !== -1 ||
  w.qmc_zn_rates(-1) !== -1 ||
  !Number.isNaN(w.qmc_zn_continuous(-1)) ||
  !Number.isNaN(w.qmc_zn_sample(0, 10, 1))
) {
  console.log("bad args accepted");
  process.exit(1);
}
console.log("bad args rejected");
