import fs from "node:fs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_lz.wasm", import.meta.url),
);
const { instance } = await WebAssembly.instantiate(bytes, {});
const w = instance.exports;
if (w._initialize) w._initialize();

const f = (x) => x.toFixed(9);
const  g = (x) => x.toFixed(6);
w.qmc_lz_reset(1, 1, 12, 1);
let rc = 0;
while (rc === 0) rc = w.qmc_lz_step(0.5);
console.log(
  `one pass t=${w.qmc_lz_time().toFixed(3)} delta=${g(w.qmc_lz_delta())} upper=${f(w.qmc_lz_upper())} exact=${f(w.qmc_lz_exact())} rc=${rc}`,
);

w.qmc_lz_reset(1, 1, 12, 2);
rc = 0;
while (rc === 0) rc = w.qmc_lz_step(0.7);
console.log(
  `two pass t=${w.qmc_lz_time().toFixed(3)} delta=${g(w.qmc_lz_delta())} upper=${f(w.qmc_lz_upper())}`,
);

w.qmc_lz_scan_rate(1, 10, 5, 0.1, 10);
const s = new Float64Array(w.memory.buffer, w.qmc_lz_scan_buffer(), 256);
console.log(`scan rate ${g(s[0])} ${g(s[1])} ${g(s[2])} ${g(s[3])} ${g(s[4])}`);
w.qmc_lz_scan_amp(1, 1, 4, 10, 11);
console.log(`scan amp ${g(s[0])} ${g(s[1])} ${g(s[2])} ${g(s[3])}`);
if (
  w.qmc_lz_reset(5, 1, 10, 1) !== -1 ||
  w.qmc_lz_scan_rate(1, 10, 1, 0.1, 1) !== -1
) {
  console.error("bad args accepted");
  process.exit(1);
}
