import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_corral.wasm", import.meta.url),
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
const F = () =>
  new Float64Array(w.memory.buffer, w.qmc_cb_field_buffer(), 6400);
console.log(`modes=${w.qmc_cb_mode_count() > 150 ? 1 : 0}`);
console.log(
  `mode0 m=${w.qmc_cb_mode_m(0)} n=${w.qmc_cb_mode_n(0)} k=${f(w.qmc_cb_mode_k(0))}  mode1 m=${w.qmc_cb_mode_m(1)} n=${w.qmc_cb_mode_n(1)} s=${w.qmc_cb_mode_sine(1)} k=${f(w.qmc_cb_mode_k(1))}  mode9 m=${w.qmc_cb_mode_m(9)} n=${w.qmc_cb_mode_n(9)} k=${f(w.qmc_cb_mode_k(9))}`,
);

w.qmc_cb_show_mode(0);
let a = F();
console.log(
  `ground f[40,40]=${f(a[40 * 80 + 40])} f[10,40]=${f(a[40 * 80 + 10])} f[0,0]=${f(a[0])}`,
);

w.qmc_cb_show_mode(5);
a = F();
console.log(
  `mode5 f[30,50]=${f(a[30 * 80 + 50])} f[50,30]=${f(a[50 * 80 + 30])}`,
);

const kept = w.qmc_cb_launch(-0.4, 0.0, 0.1, 8.0, 3.0);
const energy = w.qmc_cb_energy();
console.log(`launch kept=${f(kept, 4)} energy=${f(energy, 3)}`);

const p0 = w.qmc_cb_evolve(0.0);
a = F();
console.log(`t=0 total=${f(p0, 4)} f[40,24]=${f(a[40 * 80 + 24])}`);

const p1 = w.qmc_cb_evolve(0.05);
a = F();
console.log(
  `t=0.05 total=${f(p1, 4)} f[40,24]=${f(a[40 * 80 + 24])} f[44,45]=${f(a[44 * 80 + 45])}`,
);

if (
  w.qmc_cb_show_mode(-1) !== -1 ||
  w.qmc_cb_launch(2, 0, 0.1, 0, 0) !== -1 ||
  w.qmc_cb_launch(0, 0, 0.1, 50, 0) !== -1
) {
  console.log("bad args accepted");
  process.exit(1);
}
console.log("bad args rejected");
