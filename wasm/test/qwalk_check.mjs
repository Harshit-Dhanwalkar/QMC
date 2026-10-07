import fs from "node:fs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_qwalk.wasm", import.meta.url),
);
const { instance } = await WebAssembly.instantiate(bytes, {});
const w = instance.exports;
if (w._initialize) w._initialize();
const f = (x) => x.toFixed(9);
w.qmc_qw_reset(0.7853981633974483, 2);
w.qmc_qw_step(300);
const b = new Float64Array(w.memory.buffer, w.qmc_qw_buffer(), 2048);

console.log(
  `hadamard t=${w.qmc_qw_time()} var=${f(w.qmc_qw_variance())} ratio=${f(w.qmc_qw_variance() / 90000)} asym=${f(w.qmc_qw_asymptote())} norm=${w.qmc_qw_norm().toFixed(12)}`,
);
console.log(
  `P(0)=${f(b[512])} P(212)=${f(b[512 + 212])} classical P(0)=${f(b[1024 + 512])} P(2)=${f(b[1024 + 514])}`,
);

w.qmc_qw_reset(0.3, 0);
w.qmc_qw_step(120);
console.log(`biased mean=${f(w.qmc_qw_mean())} var=${f(w.qmc_qw_variance())}`);
if (w.qmc_qw_reset(2, 0) !== -1 || w.qmc_qw_step(1000) !== -1) {
  console.error("bad args accepted");
  process.exit(1);
}
