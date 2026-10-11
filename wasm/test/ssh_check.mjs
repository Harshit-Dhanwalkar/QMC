import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_ssh.wasm", import.meta.url),
);

let memory = null;
const { instance } = await WebAssembly.instantiate(
  bytes,
  makeWasiImports(() => memory),
);
memory = instance.exports.memory;
const w = instance.exports;
if (w._initialize) w._initialize();

const f6 = (x) => x.toFixed(6);
const f9 = (x) => x.toFixed(9);

// memory can grow, so make views after calls that touch it
const E = () => new Float64Array(w.memory.buffer, w.qmc_ssh_energies(), 40);
w.qmc_ssh_set(2, 0, 0, 1);
let e = E();
console.log(
  `topo w=2 wind=${w.qmc_ssh_winding()} zak=${f6(w.qmc_ssh_zak())} xi=${f6(w.qmc_ssh_xi())} gap=${f6(w.qmc_ssh_gap())}`,
);
console.log(
  `E19=${f9(e[19])} E20=${f9(e[20])} E0=${f6(e[0])} E39=${f6(e[39])} edge=${f6(w.qmc_ssh_edge_weight(19, 3))}`,
);

w.qmc_ssh_set(0.5, 0, 0, 1);
e = E();
console.log(
  `trivial w=0.5 wind=${w.qmc_ssh_winding()} zak=${f6(w.qmc_ssh_zak())} xi_nan=${Number.isNaN(w.qmc_ssh_xi()) ? 1 : 0} E19=${f6(e[19])} E20=${f6(e[20])}`,
);

w.qmc_ssh_set(2, 0.4, 0, 5);
e = E();
console.log(`noisy chiral E0+E39=${f6(e[0] + e[39])} E19=${f6(e[19])}`);

w.qmc_ssh_set(2, 0, 0.3, 1);
e = E();
console.log(`staggered E19=${f6(e[19])} E20=${f6(e[20])}`);

w.qmc_ssh_set(1.5, 0, 0, 1);
w.qmc_ssh_sweep(0, 3, 61);
const s = new Float64Array(w.memory.buffer, w.qmc_ssh_sweep_buffer(), 160 * 40);
console.log(
  `sweep w=0 E0=${f6(s[0])} w=1.5 s[30*40+19]=${f6(s[30 * 40 + 19])} w=3 E39=${f6(s[60 * 40 + 39])}`,
);

if (w.qmc_ssh_set(3.5, 0, 0, 1) !== -1 || w.qmc_ssh_sweep(0, 3, 1000) !== -1) {
  console.error("bad args accepted");
  process.exit(1);
}
