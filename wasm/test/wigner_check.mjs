import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_wigner.wasm", import.meta.url),
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
const W = () =>
  new Float64Array(w.memory.buffer, w.qmc_wg_w_buffer(), 128 * 128);
const MX = () => new Float64Array(w.memory.buffer, w.qmc_wg_x_buffer(), 128);
const MP = () => new Float64Array(w.memory.buffer, w.qmc_wg_p_buffer(), 128);
const PN = () => new Float64Array(w.memory.buffer, w.qmc_wg_pn_buffer(), 101);

console.log(`extent ${w.qmc_wg_extent().toFixed(4)}`);
w.qmc_wg_state(0, 2.5, 0, 0);
console.log(
  `coherent nbar=${f(w.qmc_wg_nbar())} neg=${f(w.qmc_wg_negativity())} wmax=${f(w.qmc_wg_wmax())} mx64=${f(MX()[64])}`,
);
w.qmc_wg_state(0, 2.5, 0, Math.PI / 2);
console.log(
  `kerr cat nbar=${f(w.qmc_wg_nbar())} neg=${f(w.qmc_wg_negativity())} wmin=${f(w.qmc_wg_wmin())} w[64,64]=${f(W()[64 * 128 + 64])} p64=${f(MP()[64])}`,
);
w.qmc_wg_state(1, 3, 0, 0);
console.log(
  `fock3 nbar=${f(w.qmc_wg_nbar())} wmin=${f(w.qmc_wg_wmin())} w[63,63]=${f(W()[63 * 128 + 63])} pn3=${f(PN()[3])}`,
);
w.qmc_wg_state(2, 2, 0, 0);

const p = PN();
console.log(
  `even cat neg=${f(w.qmc_wg_negativity())} pn0=${f(p[0])} pn1=${f(p[1])} pn2=${f(p[2])}`,
);
w.qmc_wg_state(4, 0.8, 0, 0);
console.log(`squeezed nbar=${f(w.qmc_wg_nbar())} wmax=${f(w.qmc_wg_wmax())}`);
if (
  w.qmc_wg_state(0, 5, 0, 0) !== -1 ||
  w.qmc_wg_state(1, 20, 0, 0) !== -1 ||
  w.qmc_wg_state(9, 1, 0, 0) !== -1
) {
  console.log("bad args accepted");
  process.exit(1);
}
console.log("bad args rejected");
