import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_bloch.wasm", import.meta.url),
);

let memory = null;
const { instance } = await WebAssembly.instantiate(
  bytes,
  makeWasiImports(() => memory),
);
memory = instance.exports.memory;
const w = instance.exports;
if (w._initialize) w._initialize();

const v = () => new Float64Array(w.memory.buffer, w.qmc_bl_vec(), 3);
w.qmc_bl_reset(0.2, 0.9, -0.3);
w.qmc_bl_step(1.3, 0.7, 0.3, 0.2, 240);
console.log(
  `damped v=(${v()[0].toFixed(9)}, ${v()[1].toFixed(9)}, ${v()[2].toFixed(9)}) t=${w.qmc_bl_time().toFixed(3)} purity=${w.qmc_bl_purity().toFixed(9)}`,
);

w.qmc_bl_reset(0, 0, 1);
w.qmc_bl_step(2.0, 0, 0, 0, 157);
console.log(
  `pulse z=${v()[2].toFixed(9)} rabi=${(1 - 2 * w.qmc_bl_rabi(w.qmc_bl_time(), 2.0, 0)).toFixed(9)}`,
);
if (w.qmc_bl_reset(2, 0, 0) !== -1 || w.qmc_bl_step(NaN, 0, 0, 0, 1) !== -1) {
  console.error("bad args accepted");
  process.exit(1);
}
