import fs from "node:fs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_dirac.wasm", import.meta.url),
);
const { instance } = await WebAssembly.instantiate(bytes, {});
const w = instance.exports;
if (w._initialize) w._initialize();
for (const v0 of [0.0, 0.5, 4.0, 10.0]) {
  w.qmc_dc_klein(Math.sqrt(3), v0);
  w.qmc_dc_step(2500);
  console.log(
    `klein V0=${v0.toFixed(1)} norm=${w.qmc_dc_norm().toFixed(10)} T=${w.qmc_dc_right().toFixed(10)}`,
  );
}
w.qmc_dc_zitter();
w.qmc_dc_step(39);
console.log(`zitter x=${w.qmc_dc_position().toFixed(10)}`);

// invalid input must be rejected, not crash
if (w.qmc_dc_klein(-1, 0) !== -1 || w.qmc_dc_klein(1, 99) !== -1) {
  console.error("bad args accepted");
  process.exit(1);
}
