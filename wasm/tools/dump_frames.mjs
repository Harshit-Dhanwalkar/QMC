// Dump density frames from the WASM solver for wasm/tools/make_gif.py.
//   node wasm/tools/dump_frames.mjs /tmp/qmc_frames.bin
import fs from "node:fs";
const out = process.argv[2] || "/tmp/qmc_frames.bin";
const N = 128,
  BOX = 40,
  DT = 0.02,
  FRAMES = 100,
  STEPS = 6;
const { instance } = await WebAssembly.instantiate(
  fs.readFileSync(
    new URL("../../docs/src/playground/qmc.wasm", import.meta.url),
  ),
  {},
);
const w = instance.exports;
if (w._initialize) w._initialize();
w.qmc_init(N, BOX, DT);
w.qmc_build_slits(2.0, 1.0, 2, 6.0, 1.6, 60.0);
w.qmc_wavepacket(-12, 0, 3.0, 0, 2.0);
const parts = [
  Buffer.from(
    new Float64Array(w.memory.buffer, w.qmc_potential(), N * N).slice().buffer,
  ),
];
for (let f = 0; f < FRAMES; f++) {
  w.qmc_step(STEPS);
  parts.push(
    Buffer.from(
      new Float32Array(w.memory.buffer, w.qmc_density(), N * N).slice().buffer,
    ),
  );
}
fs.writeFileSync(out, Buffer.concat(parts));
console.log(`wrote ${FRAMES} frames to ${out}`);
