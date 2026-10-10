import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_rotor.wasm", import.meta.url),
);

let memory = null;
const { instance } = await WebAssembly.instantiate(
  bytes,
  makeWasiImports(() => memory),
);
memory = instance.exports.memory;
const w = instance.exports;
if (w._initialize) w._initialize();

const f = (x) => x.toFixed(9);
w.qmc_kr_reset(5, 1, 7);
w.qmc_kr_step(100);
w.qmc_kr_step(100);

// memory can grow, so make views after calls that touch it
const q = new Float64Array(w.memory.buffer, w.qmc_kr_qbuffer(), 2048);
console.log(
  `local t=${w.qmc_kr_time()} qm2=${f(w.qmc_kr_qm2())} norm=${w.qmc_kr_norm().toFixed(10)} D=${w.qmc_kr_diffusion().toFixed(6)} edge_ok=${w.qmc_kr_edge() < 1e-20 ? 1 : 0}`,
);
console.log(`q[1024]=${f(q[1024])} q[1030]=${f(q[1030])}`);

// classical map is chaotic: compare it only after a few kicks
w.qmc_kr_reset(5, 1, 7);
w.qmc_kr_step(8);
const c = new Float64Array(w.memory.buffer, w.qmc_kr_cbuffer(), 6000);
console.log(
  `classical t=${w.qmc_kr_time()} cm2=${w.qmc_kr_cm2().toFixed(5)} c0=(${c[0].toFixed(5)},${c[1].toFixed(5)})`,
);

w.qmc_kr_reset(5, 4 * Math.PI, 1);
w.qmc_kr_step(10);
const x = (10 * 5) / (4 * Math.PI);
console.log(`resonance qm2=${f(w.qmc_kr_qm2())} exact=${f((x * x) / 2)}`);
if (w.qmc_kr_reset(13, 1, 1) !== -1 || w.qmc_kr_step(1000) !== -1) {
  console.error("bad args accepted");
  process.exit(1);
}
