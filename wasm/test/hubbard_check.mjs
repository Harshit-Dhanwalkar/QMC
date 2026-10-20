import fs from "node:fs";
import { makeWasiImports } from "../wasi-shim.mjs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_hubbard.wasm", import.meta.url),
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
const buf = (fn, n) => new Float64Array(w.memory.buffer, fn(), n);

w.qmc_hb_curves(12.0);
const sp = buf(w.qmc_hb_spec_buffer, 6 * 201);
const d = buf(w.qmc_hb_gs_d_buffer, 201);
const ss = buf(w.qmc_hb_gs_ss_buffer, 201);
const se = buf(w.qmc_hb_gs_s_buffer, 201);
console.log(
  `spec l0 i0=${f(sp[0])} i100=${f(sp[100])} i200=${f(sp[200])} l5 i100=${f(sp[5 * 201 + 100])}`,
);
console.log(`gs d i0=${f(d[0])} i50=${f(d[50])} i200=${f(d[200])}`);
console.log(`gs ss i0=${f(ss[0])} i50=${f(ss[50])} i200=${f(ss[200])}`);
console.log(`gs s i0=${f(se[0])} i50=${f(se[50])} i200=${f(se[200])}`);

const e0 = w.qmc_hb_levels(4.0),
  lv = buf(w.qmc_hb_levels_buffer, 6);
console.log(
  `levels e0=${f(e0)} ${f(lv[1])} ${f(lv[2])} ${f(lv[3])} ${f(lv[4])} ${f(lv[5])}`,
);

w.qmc_hb_quench(8.0, 20.0);
const qd = buf(w.qmc_hb_q_d_buffer, 401);
const qm = buf(w.qmc_hb_q_m_buffer, 401);
const qs = buf(w.qmc_hb_q_s_buffer, 401);
const qss = buf(w.qmc_hb_q_ss_buffer, 401);
console.log(`quench d q0=${f(qd[0])} q100=${f(qd[100])} q400=${f(qd[400])}`);
console.log(`quench m q0=${f(qm[0])} q100=${f(qm[100])} q400=${f(qm[400])}`);
console.log(
  `quench s q100=${f(qs[100])} q300=${f(qs[300])} ss q100=${f(qss[100])}`,
);

w.qmc_hb_probs(4.0, 0, 0.0);
let p = buf(w.qmc_hb_probs_buffer, 8);
console.log(
  `probs gs ${f(p[0])} ${f(p[1])} ${f(p[2])} ${f(p[3])} | ${f(p[4])} ${f(p[5])} ${f(p[6])} ${f(p[7])}`,
);

w.qmc_hb_probs(4.0, 1, 1.3);
p = buf(w.qmc_hb_probs_buffer, 8);
console.log(
  `probs neel ${f(p[0])} ${f(p[1])} ${f(p[2])} ${f(p[3])} | ${f(p[4])} ${f(p[5])} ${f(p[6])} ${f(p[7])}`,
);

w.qmc_hb_thermal(4.0, 5.0);
const td = buf(w.qmc_hb_th_d_buffer, 201);
const tss = buf(w.qmc_hb_th_ss_buffer, 201);
const ts = buf(w.qmc_hb_th_s_buffer, 201);
const te = buf(w.qmc_hb_th_e_buffer, 201);
console.log(
  `thermal d ${f(td[0])} ${f(td[100])} ${f(td[200])} ss ${f(tss[0])} ${f(tss[100])} ${f(tss[200])}`,
);
console.log(
  `thermal s ${f(ts[0])} ${f(ts[100])} ${f(ts[200])} e ${f(te[0])} ${f(te[100])} ${f(te[200])}`,
);

if (
  w.qmc_hb_curves(0) !== -1 ||
  w.qmc_hb_quench(-1, 5) !== -1 ||
  w.qmc_hb_quench(2, 0) !== -1 ||
  w.qmc_hb_probs(2, 2, 0) !== -1 ||
  w.qmc_hb_thermal(2, 0) !== -1 ||
  !Number.isNaN(w.qmc_hb_levels(99))
) {
  console.log("bad args accepted");
  process.exit(1);
}
console.log("bad args rejected");
