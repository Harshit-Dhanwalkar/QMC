// Growable GPU batches. Each class owns both CPU-side data and a GL
// buffer and exposes push / reset / flush
//
// flush(gl, program, uniforms) draws everything pushed since last reset()
// `uniforms` is a plain object: each key is set with matching gl.uniform*
// call inferred from value's length (2 -> 2f, 4 -> 4f, etc)

function setUniform(gl, loc, v) {
  if (!loc) {
    return;
  }

  if (Array.isArray(v) || ArrayBuffer.isView(v)) {
    switch (v.length) {
      case 1:
        gl.uniform1f(loc, v[0]);
        break;
      case 2:
        gl.uniform2f(loc, v[0], v[1]);
        break;
      case 3:
        gl.uniform3f(loc, v[0], v[1], v[2]);
        break;
      case 4:
        gl.uniform4f(loc, v[0], v[1], v[2], v[3]);
        break;
      default:
        console.warn("Uniform array of length", v.length, "unsupported");
    }
  } else if (typeof v === "number") {
    gl.uniform1f(loc, v);
  }
}

class InstanceBatchBase {
  constructor(gl, floatsPerInstance, attributes) {
    this.gl = gl;
    this.floatsPerInstance = floatsPerInstance;
    this.attributes = attributes;
    this.program = null;
    // Cache uniform locations per program using a Map to prevent
    // cache-loss on program switches
    this.programUniformLocs = new Map();
    this.buffer = gl.createBuffer();
    this.vao = gl.createVertexArray();
    this.data = new Float32Array(4096 * floatsPerInstance);
    this.count = 0;
    // GPU buffer is re-uploaded only when CPU data changed
    // so pan / zoom / hover redraws cost zero buffer traffic
    this.dirty = true;
    this.gpuBytes = 0;
  }

  reset() {
    this.count = 0;
    this.dirty = true;
  }

  _ensure(n) {
    const need = n * this.floatsPerInstance;
    if (this.data.length < need) {
      let cap = Math.max(this.data.length, 1024);
      while (cap < need) {
        cap *= 2;
      }

      const grown = new Float32Array(cap);
      grown.set(this.data.subarray(0, this.count * this.floatsPerInstance));
      this.data = grown;
    }
  }

  _write(vals) {
    this._ensure(this.count + 1);
    const o = this.count * this.floatsPerInstance;
    for (let i = 0; i < this.floatsPerInstance; i++) {
      this.data[o + i] = vals[i];
    }
    this.count++;
  }

  // Binds VAO, uploads data, sets uniforms, draws. `uniforms` is a plain
  // object, `setup` is optional (gl, program) -> void for extra state
  flush(gl, program, uniforms = {}, setup = null) {
    if (this.count === 0) {
      return;
    }

    gl.useProgram(program);
    gl.bindVertexArray(this.vao);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.buffer);

    if (this.dirty) {
      const view = this.data.subarray(0, this.count * this.floatsPerInstance);
      if (view.byteLength > this.gpuBytes) {
        // Grow-only: reallocate (bufferData) when it no longer fits
        this.gpuBytes = this.data.byteLength;
        gl.bufferData(gl.ARRAY_BUFFER, this.data, gl.DYNAMIC_DRAW);
      } else {
        gl.bufferSubData(gl.ARRAY_BUFFER, 0, view);
      }

      this.dirty = false;
    }

    // Attribute setup on first use per program
    if (this.program !== program) {
      for (const attr of this.attributes) {
        const loc = gl.getAttribLocation(program, attr.name);
        if (loc < 0) {
          continue;
        }

        gl.enableVertexAttribArray(loc);
        gl.vertexAttribPointer(
          loc,
          attr.size,
          gl.FLOAT,
          false,
          this.floatsPerInstance * 4,
          attr.offset * 4,
        );

        if (attr.divisor !== undefined) {
          gl.vertexAttribDivisor(loc, attr.divisor);
        }
      }

      this.program = program;
    }

    // Retrieve or initialize uniform location cache for active program
    let uniformLocs = this.programUniformLocs.get(program);
    if (!uniformLocs) {
      uniformLocs = {};
      this.programUniformLocs.set(program, uniformLocs);
    }

    for (const [k, v] of Object.entries(uniforms)) {
      if (!(k in uniformLocs)) {
        uniformLocs[k] = gl.getUniformLocation(program, k);
      }

      setUniform(gl, uniformLocs[k], v);
    }

    if (setup) {
      setup(gl, program);
    }

    this._draw(gl);
    gl.bindVertexArray(null);
  }

  dispose() {
    this.gl.deleteBuffer(this.buffer);
    this.gl.deleteVertexArray(this.vao);
  }
}

// RectBatch: instanced unit quad. Each instance carries
//   (x0, y0, x1, y1) in data coords + (r, g, b, a)
// Vertex shader mixes a_quad into a_rect and maps through u_view
export class RectBatch extends InstanceBatchBase {
  constructor(gl) {
    super(gl, 8, [
      { name: "a_rect", size: 4, offset: 0, divisor: 1 },
      { name: "a_color", size: 4, offset: 4, divisor: 1 },
    ]);

    // Unit quad: TRIANGLE_STRIP (0,0) (1,0) (0,1) (1,1)
    this._quadBuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this._quadBuf);
    gl.bufferData(
      gl.ARRAY_BUFFER,
      new Float32Array([0, 0, 1, 0, 0, 1, 1, 1]),
      gl.STATIC_DRAW,
    );
    gl.bindBuffer(gl.ARRAY_BUFFER, null);
  }

  push(x0, y0, x1, y1, rgba) {
    this._write([x0, y0, x1, y1, rgba[0], rgba[1], rgba[2], rgba[3]]);
  }

  flush(gl, program, uniforms = {}, setup = null) {
    // a_quad lives at whatever location the linker picks; look it up once
    // per program and bind it into the shared VAO
    if (this.program !== program) {
      gl.bindVertexArray(this.vao);
      const qLoc = gl.getAttribLocation(program, "a_quad");
      if (qLoc >= 0) {
        gl.bindBuffer(gl.ARRAY_BUFFER, this._quadBuf);
        gl.enableVertexAttribArray(qLoc);
        gl.vertexAttribPointer(qLoc, 2, gl.FLOAT, false, 0, 0);
        gl.vertexAttribDivisor(qLoc, 0);
      }

      gl.bindVertexArray(null);
    }

    super.flush(gl, program, uniforms, setup);
  }

  _draw(gl) {
    gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, this.count);
  }

  dispose() {
    this.gl.deleteBuffer(this._quadBuf);
    super.dispose();
  }
}

// PointBatch: plain gl.POINTS. Each vertex is
//   (x, y, size, r, g, b, a)
// Size is a multiplier on gl_PointSize in shader
export class PointBatch extends InstanceBatchBase {
  constructor(gl) {
    super(gl, 7, [
      { name: "a_pos", size: 2, offset: 0 },
      { name: "a_size", size: 1, offset: 2 },
      { name: "a_color", size: 4, offset: 3 },
    ]);
  }

  push(x, y, size, rgba) {
    this._write([x, y, size, rgba[0], rgba[1], rgba[2], rgba[3]]);
  }

  _draw(gl) {
    gl.drawArrays(gl.POINTS, 0, this.count);
  }
}

// LineBatch: gl.LINES, each vertex is (x, y, r, g, b, a)
// Caller pushes vertices in pairs; each pair is one segment
// Line width in WebGL2 is capped at 1 on most drivers; use RectBatch if
// thicker strokes are needed
export class LineBatch extends InstanceBatchBase {
  constructor(gl) {
    super(gl, 6, [
      { name: "a_pos", size: 2, offset: 0 },
      { name: "a_color", size: 4, offset: 2 },
    ]);
  }

  pushVertex(x, y, rgba) {
    this._write([x, y, rgba[0], rgba[1], rgba[2], rgba[3]]);
  }

  pushSegment(x0, y0, x1, y1, rgba) {
    this.pushVertex(x0, y0, rgba);
    this.pushVertex(x1, y1, rgba);
  }

  _draw(gl) {
    gl.drawArrays(gl.LINES, 0, this.count);
  }
}
