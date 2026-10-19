import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_carpet.wasm", import.meta.url),
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
const CAR = () =>
  new Float64Array(w.memory.buffer, w.qmc_qc_carpet_buffer(), 384 * 192);
const SL = () =>
  new Float64Array(w.memory.buffer, w.qmc_qc_slice_buffer(), 192);
const CV = () =>
  new Float64Array(w.memory.buffer, w.qmc_qc_curve_buffer(), 513);

w.qmc_qc_set(0.3, 0.05, 20);
let car = CAR();
let cv = CV();
console.log(
  `weight=${f(w.qmc_qc_weight())} energy=${w.qmc_qc_energy().toFixed(4)}`,
);
console.log(
  `carpet t0 x57=${f(car[57])} x134=${f(car[134])} t192 x40=${f(car[192 * 192 + 40])} t383 x57=${f(car[383 * 192 + 57])}`,
);
console.log(
  `curve c0=${f(cv[0])} c128=${f(cv[128])} c256=${f(cv[256])} c512=${f(cv[512])}`,
);

let r = w.qmc_qc_slice(0.5);
let sl = SL();
console.log(`slice half r=${f(r)} s57=${f(sl[57])} s134=${f(sl[134])}`);

r = w.qmc_qc_slice(0.25);
sl = SL();
console.log(`slice quarter r=${f(r)} s57=${f(sl[57])} s134=${f(sl[134])}`);

w.qmc_qc_set(0.5, 0.1, 0);
const half = w.qmc_qc_slice(0.5);
console.log(
  `centred weight=${f(w.qmc_qc_weight())} energy=${w.qmc_qc_energy().toFixed(4)} half=${f(half)}`,
);

if (
  w.qmc_qc_set(1.5, 0.05, 1) !== -1 ||
  w.qmc_qc_set(0.5, 0.001, 1) !== -1 ||
  !Number.isNaN(w.qmc_qc_slice(NaN))
) {
  console.log("bad args accepted");
  process.exit(1);
}
console.log("bad args rejected");
