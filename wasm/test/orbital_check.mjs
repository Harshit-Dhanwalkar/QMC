import fs from "node:fs";

const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_orbital.wasm", import.meta.url),
);
const { instance } = await WebAssembly.instantiate(bytes, {});
const w = instance.exports;
if (w._initialize) w._initialize();
const orbs = [
  [1, 0, 0, 0],
  [2, 1, 1, 0],
  [3, 2, -1, 1],
  [4, 3, 2, 1],
  [6, 2, 0, 0],
];
orbs.forEach(([n, l, m, r], k) => {
  if (w.qmc_orb_sample(n, l, m, r, 2000, 7 + k) !== 0) {
    process.exit(1);
  }

  const p = new Float64Array(w.memory.buffer, w.qmc_orb_buffer(), 8000);
  let sr = 0;
  let sz = 0;
  let sp = 0;
  for (let i = 0; i < 2000; i++) {
    sr += Math.hypot(p[4 * i], p[4 * i + 1], p[4 * i + 2]);
    sz += p[4 * i + 2];
    sp += p[4 * i + 3];
  }

  const sg = m < 0 ? "-" : "+";
  console.log(
    `orb ${n}${l}${sg}${Math.abs(m)} r=${sr.toFixed(6)} z=${sz.toFixed(6)} ph=${sp.toFixed(6)}`,
  );
});

console.log(`E2=${w.qmc_orb_energy_ev(2).toFixed(6)}`);
if (
  w.qmc_orb_sample(1, 0, 0, 0, 0, 1) !== -1 ||
  w.qmc_orb_sample(9, 0, 0, 0, 10, 1) === 0
) {
  console.error("bad args accepted");
  process.exit(1);
}
