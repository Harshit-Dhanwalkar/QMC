// Shared WASI shim for the QMC playground demos
//
//   <script src="wasi-shim.js"></script>
//   ...
//   let wasmMemory = null;
//   const { instance } = await WebAssembly.instantiate(
//     bytes, makeWasiImports(() => wasmMemory));
//   wasmMemory = instance.exports.memory;
//
(function () {
  "use strict";

  function makeWasiImports(getMemory) {
    const decoder = new TextDecoder();
    const view = () => new DataView(getMemory().buffer);

    function writeToFd(fd, bytes) {
      const text = decoder.decode(bytes);
      if (fd === 2) console.warn("[wasm]", text);
      else console.log("[wasm]", text);
    }

    return {
      wasi_snapshot_preview1: {
        fd_write: (fd, iovsPtr, iovsLen, nwrittenPtr) => {
          const mem = getMemory();
          if (!mem) return 0;
          const dv = view();
          let written = 0;
          for (let i = 0; i < iovsLen; i++) {
            const ptr = dv.getUint32(iovsPtr + i * 8, true);
            const len = dv.getUint32(iovsPtr + i * 8 + 4, true);
            writeToFd(fd, new Uint8Array(mem.buffer, ptr, len));
            written += len;
          }
          if (nwrittenPtr) dv.setUint32(nwrittenPtr, written, true);
          return 0;
        },
        fd_close: () => 0,
        fd_seek: () => 0,
        fd_fdstat_get: () => 0,
        fd_prestat_get: () => 8,
        fd_prestat_dir_name: () => 8,

        args_sizes_get: (argcPtr, bufPtr) => {
          if (!getMemory()) return 0;
          const dv = view();
          dv.setUint32(argcPtr, 0, true);
          dv.setUint32(bufPtr, 0, true);
          return 0;
        },

        args_get: () => 0,
        environ_sizes_get: (countPtr, sizePtr) => {
          if (!getMemory()) return 0;
          const dv = view();
          dv.setUint32(countPtr, 0, true);
          dv.setUint32(sizePtr, 0, true);
          return 0;
        },

        environ_get: () => 0,

        random_get: (ptr, len) => {
          const mem = getMemory();
          if (!mem) return 0;
          const buf = new Uint8Array(mem.buffer, ptr, len);
          if (typeof crypto !== "undefined" && crypto.getRandomValues) {
            crypto.getRandomValues(buf);
          } else {
            for (let i = 0; i < len; i++) buf[i] = (Math.random() * 256) | 0;
          }
          return 0;
        },

        clock_time_get: (clockId, precision, timePtr) => {
          const mem = getMemory();
          if (!mem) return 0;
          const ns = BigInt(Date.now()) * 1000000n;
          new DataView(mem.buffer).setBigUint64(timePtr, ns, true);
          return 0;
        },

        proc_exit: (code) => {
          throw new Error("WASI proc_exit(" + code + ")");
        },

        proc_raise: () => 0,
        sched_yield: () => 0,
      },
    };
  }

  window.makeWasiImports = makeWasiImports;
})();
