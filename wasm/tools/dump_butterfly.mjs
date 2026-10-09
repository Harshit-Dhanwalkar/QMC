// Compute the Hofstadter butterfly (band edges + TKNN gap labels) with the
// WASM module and dump it as JSON for wasm/tools/make_butterfly_png.py.
//   node wasm/tools/dump_butterfly.mjs 50 /tmp/butterfly.json
import fs from "node:fs";
const QMAX = parseInt(process.argv[2] || "50", 10);
const out = process.argv[3] || "/tmp/butterfly.json";
const bytes = fs.readFileSync(
  new URL("../../docs/src/playground/qmc_butterfly.wasm", import.meta.url),
);
let mem;
const wasi = {
  fd_close: () => 0,
  fd_seek: () => 70,
  fd_write: (fd, iov, n, o) => {
    const v = new DataView(mem.buffer);
    let t = 0;
    for (let i = 0; i < n; i++) {
      t += v.getUint32(iov + 8 * i + 4, true);
    }
    v.setUint32(o, t, true);
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
const rows = [];
const t0 = performance.now();
for (let q = 2; q <= QMAX; q++) {
  for (let p = 1; p < q; p++) {
    if (gcd(p, q) !== 1) {
      continue;
    }
    w.qmc_bf_edges(p, q);
    const g = Array.from(
      new Float64Array(w.memory.buffer, w.qmc_bf_buffer(), 2 * q),
    );
    const t = [];
    for (let r = 1; r < q; r++) {
      const v = w.qmc_bf_tknn(p, q, r);
      t.push(v === INT_MIN ? null : v);
    }
    rows.push({ p, q, edges: g, tknn: t });
  }
}

fs.writeFileSync(out, JSON.stringify(rows));
console.log(
  `q<=${QMAX}: ${rows.length} fluxes in ${((performance.now() - t0) / 1000).toFixed(2)} s -> ${out}`,
);
