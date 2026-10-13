import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_grover.wasm", import.meta.url),
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
const A = () => new Float64Array(w.memory.buffer, w.qmc_gv_amp_buffer(), 1024);
const C = () => new Float64Array(w.memory.buffer, w.qmc_gv_curve_buffer(), 400);
const E = () => new Float64Array(w.memory.buffer, w.qmc_gv_exact_buffer(), 400);

w.qmc_gv_set_n(10);
w.qmc_gv_mark(37, 1);
console.log(`M=1 theta=${f(w.qmc_gv_theta())} kopt=${w.qmc_gv_kopt()}`);
w.qmc_gv_goto(10);

let a = A();
console.log(
  `k=${w.qmc_gv_iter()} p=${f(w.qmc_gv_success())} exact=${f(w.qmc_gv_exact())} a37=${f(a[37])} a0=${f(a[0])}`,
);
w.qmc_gv_goto(25);

a = A();
console.log(
  `k=${w.qmc_gv_iter()} p=${f(w.qmc_gv_success())} exact=${f(w.qmc_gv_exact())} a37=${f(a[37])} a0=${f(a[0])}`,
);

w.qmc_gv_mark(100, 1);
w.qmc_gv_mark(500, 1);
w.qmc_gv_mark(900, 1);
w.qmc_gv_goto(w.qmc_gv_kopt());
console.log(
  `M=${w.qmc_gv_marks()} kopt=${w.qmc_gv_kopt()} p=${f(w.qmc_gv_success())}`,
);

w.qmc_gv_curve(60);
const c = C();
const e = E();
console.log(
  `curve c5=${f(c[5])} c12=${f(c[12])} c40=${f(c[40])} e12=${f(e[12])}`,
);

w.qmc_gv_set_n(4);
w.qmc_gv_mark(3, 1);
console.log(`step_one_mark=${w.qmc_gv_step()}`);
for (let i = 0; i < 8; i++) {
  w.qmc_gv_mark(i, 1);
}
console.log(`marks=${w.qmc_gv_marks()} (capped below N/2)`);
if (w.qmc_gv_set_n(1) !== -1 || w.qmc_gv_mark(99, 1) !== -1) {
  console.log("bad args accepted");
  process.exit(1);
}
console.log("bad args rejected");
