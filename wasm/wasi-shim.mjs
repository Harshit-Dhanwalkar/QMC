// Usage:
//   import { makeWasiImports } from "./wasi-shim.mjs";
//   let memory = null;
//   const { instance } = await WebAssembly.instantiate(
//     bytes,
//     makeWasiImports(() => memory),
//   );
//   memory = instance.exports.memory;

export async function instantiate(bytes) {
  let memory = null;
  const { instance } = await WebAssembly.instantiate(
    bytes,
    makeWasiImports(() => memory),
  );
  memory = instance.exports.memory;
  if (instance.exports._initialize) {
    instance.exports._initialize();
  }
  return instance.exports;
}

export function makeWasiImports(getMemory) {
  const decoder = new TextDecoder();

  const view = () => new DataView(getMemory().buffer);

  function writeToFd(fd, bytes) {
    // fd 1 = stdout, fd 2 = stderr
    const text = decoder.decode(bytes);
    if (fd === 2) {
      console.error(text);
    } else {
      console.log(text);
    }
  }

  return {
    wasi_snapshot_preview1: {
      // stdio
      fd_write: (fd, iovsPtr, iovsLen, nwrittenPtr) => {
        const mem = getMemory();
        if (!mem) {
          return 0;
        }

        const dv = view();
        let written = 0;
        for (let i = 0; i < iovsLen; i++) {
          const ptr = dv.getUint32(iovsPtr + i * 8, true);
          const len = dv.getUint32(iovsPtr + i * 8 + 4, true);
          writeToFd(fd, new Uint8Array(mem.buffer, ptr, len));
          written += len;
        }

        if (nwrittenPtr) {
          dv.setUint32(nwrittenPtr, written, true);
        }

        return 0;
      },

      fd_close: () => 0,
      fd_seek: () => 0,
      fd_fdstat_get: () => 0,
      fd_prestat_get: () => 8,
      fd_prestat_dir_name: () => 8,

      // args / environ (reactor: none)
      args_sizes_get: (argcPtr, bufPtr) => {
        if (!getMemory()) return 0;
        const dv = view();
        dv.setUint32(argcPtr, 0, true);
        dv.setUint32(bufPtr, 0, true);

        return 0;
      },

      args_get: () => 0,
      environ_sizes_get: (countPtr, sizePtr) => {
        if (!getMemory()) {
          return 0;
        }
        const dv = view();
        dv.setUint32(countPtr, 0, true);
        dv.setUint32(sizePtr, 0, true);

        return 0;
      },
      environ_get: () => 0,

      // random: WASI expects real entropy; browsers have crypto
      random_get: (ptr, len) => {
        const mem = getMemory();
        if (!mem) {
          return 0;
        }

        const buf = new Uint8Array(mem.buffer, ptr, len);
        if (typeof globalThis.crypto?.getRandomValues === "function") {
          globalThis.crypto.getRandomValues(buf);
        } else {
          for (let i = 0; i < len; i++) {
            buf[i] = (Math.random() * 256) | 0;
          }
        }

        return 0;
      },

      // clock: nanosecond monotonic + realtime
      clock_time_get: (clockId, precision, timePtr) => {
        const mem = getMemory();
        if (!mem) {
          return 0;
        }
        const ns = BigInt(Date.now()) * 1_000_000n;
        new DataView(mem.buffer).setBigUint64(timePtr, ns, true);

        return 0;
      },

      // process control
      proc_exit: (code) => {
        throw new Error(`WASI proc_exit(${code})`);
      },
      proc_raise: () => 0,
      sched_yield: () => 0,
    },
  };
}
