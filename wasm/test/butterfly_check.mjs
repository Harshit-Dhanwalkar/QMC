import fs from "node:fs";
const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_butterfly.wasm", import.meta.url),
);
let mem;
const wasi = {
  fd_close: () => 0,
  fd_seek: () => 70,
  fd_write: (fd, iov, n, out) => {
    // swallow stdio, report bytes "written"
    const v = new DataView(mem.buffer);
    let tot = 0;
    for (let i = 0; i < n; i++) tot += v.getUint32(iov + 8 * i + 4, true);
    v.setUint32(out, tot, true);
    return 0;
  },
};
const { instance } = await WebAssembly.instantiate(bytes, {
  wasi_snapshot_preview1: wasi,
});
const w = instance.exports;
mem = w.memory;
if (w._initialize) w._initialize();
const gcd = (a, b) => (b ? gcd(b, a % b) : a);
const INT_MIN = -2147483648;
let sum = 0,
  tsum = 0;
for (let q = 2; q <= 24; q++)
  for (let p = 1; p < q; p++) {
    if (gcd(p, q) !== 1) continue;
    if (w.qmc_bf_edges(p, q) !== 0) throw new Error(`edges ${p}/${q}`);
    const g = new Float64Array(w.memory.buffer, w.qmc_bf_buffer(), 2 * q);
    for (let b = 0; b < q; b++)
      sum +=
        (g[2 * b + 1] - g[2 * b]) * (b + 1) + g[2 * b] * g[2 * b] * (q - b);
    for (let r = 1; r < q; r++) {
      const t = w.qmc_bf_tknn(p, q, r);
      if (t !== INT_MIN) tsum += Math.abs(t) * r;
    }
  }
console.log(`digest edges=${sum.toFixed(9)} tknn=${tsum}`);
// live spot checks (same calls the page makes on click)
const c = w.qmc_bf_chern(1, 3, 1, 16);
console.error(
  `spot check 1/3, gap 1: numeric C = ${c.toFixed(3)}, TKNN = ${w.qmc_bf_tknn(1, 3, 1)}`,
);
if (Math.abs(c - w.qmc_bf_tknn(1, 3, 1)) > 0.05) process.exit(1);
if (!Number.isNaN(w.qmc_bf_chern(1, 4, 2, 16))) {
  console.error("expected closed middle gap at 1/4");
  process.exit(1);
}
