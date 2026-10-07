import fs from "node:fs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_anderson.wasm", import.meta.url),
);
const { instance } = await WebAssembly.instantiate(bytes, {});
const w = instance.exports;
if (w._initialize) w._initialize();

const f = (x) => x.toFixed(9);
w.qmc_an_reset(0, 1);
w.qmc_an_step(1000);
console.log(
  `clean t=${w.qmc_an_time().toFixed(2)} width=${f(w.qmc_an_width())} exact=${f(Math.SQRT2 * w.qmc_an_time())} norm=${f(w.qmc_an_norm())}`,
);
w.qmc_an_reset(6, 2024);
w.qmc_an_step(3000);

const b = new Float64Array(w.memory.buffer, w.qmc_an_buffer(), 2 * 801);
console.log(
  `disordered t=${w.qmc_an_time().toFixed(2)} width=${f(w.qmc_an_width())} ipr=${f(w.qmc_an_ipr())} energy=${f(w.qmc_an_energy())} norm=${f(w.qmc_an_norm())}`,
);
console.log(`centre=${f(b[400])} eps0=${f(b[801])}`);
console.log(`xi(W=1,E=1)=${w.qmc_an_xi(1, 1).toFixed(6)}`);
if (w.qmc_an_reset(-1, 1) !== -1 || !Number.isNaN(w.qmc_an_xi(1, 2.5))) {
  console.error("bad args accepted");
  process.exit(1);
}
