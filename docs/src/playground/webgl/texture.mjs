// 2D image texture renderer. Uploads RGBA8 pixels and draws them into
// plot area with a shader that maps visible view (u_view) into full data
// extent (u_data) to compute UVs
//
// Use this for scalar-field demos: quantum carpet, Wigner function, Anderson
// density, Ising correlations, Dirac density
//
// Colours are baked in on CPU before upload, so no LUT texture is
// needed. If you later want shader-side colormapping, add a u_lut
// sampler and remap s.r instead of sampling s directly

export const IMAGE_VS = `#version 300 es
in vec2 a_quad;
uniform vec4 u_view;
uniform vec4 u_data;
uniform vec4 u_plot;
uniform vec2 u_canvas;
out vec2 v_uv;
void main() {
  float dx = mix(u_view.x, u_view.z, a_quad.x);
  float dy = mix(u_view.y, u_view.w, a_quad.y);
  v_uv = (vec2(dx, dy) - u_data.xy) / (u_data.zw - u_data.xy);
  float px = u_plot.x + a_quad.x * u_plot.z;
  float py = u_plot.y + (1.0 - a_quad.y) * u_plot.w;
  gl_Position = vec4(
    (px / u_canvas.x) * 2.0 - 1.0,
    1.0 - (py / u_canvas.y) * 2.0,
    0.0, 1.0);
}`;

export const IMAGE_FS = `#version 300 es
precision highp float;
in vec2 v_uv;
uniform sampler2D u_tex;
out vec4 outColor;
void main() {
  outColor = texture(u_tex, v_uv);
}`;

export class Image2D {
  constructor(gl) {
    this.gl = gl;
    this.texture = gl.createTexture();
    this.quadBuf = gl.createBuffer();
    this.vao = gl.createVertexArray();
    this.width = 0;
    this.height = 0;
    this.program = null;
    this.uniformLocs = {};

    gl.bindVertexArray(this.vao);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.quadBuf);
    gl.bufferData(
      gl.ARRAY_BUFFER,
      new Float32Array([0, 0, 1, 0, 0, 1, 1, 1]),
      gl.STATIC_DRAW,
    );
    gl.bindVertexArray(null);
  }

  // pixels: Uint8ClampedArray of length width * height * 4 (RGBA)
  // Row 0 of input ends up at TOP of drawn image because we set
  // UNPACK_FLIP_Y_WEBGL below
  upload(pixels, width, height) {
    const gl = this.gl;
    this.width = width;
    this.height = height;
    gl.bindTexture(gl.TEXTURE_2D, this.texture);
    gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
    gl.texImage2D(
      gl.TEXTURE_2D,
      0,
      gl.RGBA8,
      width,
      height,
      0,
      gl.RGBA,
      gl.UNSIGNED_BYTE,
      pixels,
    );
    gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  }

  flush(gl, program, uniforms = {}, setup = null) {
    gl.useProgram(program);
    gl.bindVertexArray(this.vao);

    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.texture);

    if (this.program !== program) {
      const qLoc = gl.getAttribLocation(program, "a_quad");
      if (qLoc >= 0) {
        gl.bindBuffer(gl.ARRAY_BUFFER, this.quadBuf);
        gl.enableVertexAttribArray(qLoc);
        gl.vertexAttribPointer(qLoc, 2, gl.FLOAT, false, 0, 0);
        gl.vertexAttribDivisor(qLoc, 0);
      }

      this.program = program;
      this.uniformLocs = {};
    }

    if (!("u_tex" in this.uniformLocs)) {
      this.uniformLocs.u_tex = gl.getUniformLocation(program, "u_tex");
    }

    gl.uniform1i(this.uniformLocs.u_tex, 0);

    for (const [k, v] of Object.entries(uniforms)) {
      if (!(k in this.uniformLocs)) {
        this.uniformLocs[k] = gl.getUniformLocation(program, k);
      }

      const loc = this.uniformLocs[k];
      if (!loc) {
        continue;
      }

      if (typeof v === "number") {
        gl.uniform1f(loc, v);
      } else if (v.length === 2) {
        gl.uniform2f(loc, v[0], v[1]);
      } else if (v.length === 3) {
        gl.uniform3f(loc, v[0], v[1], v[2]);
      } else if (v.length === 4) {
        gl.uniform4f(loc, v[0], v[1], v[2], v[3]);
      }
    }

    if (setup) {
      setup(gl, program);
    }
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    gl.bindVertexArray(null);
  }

  dispose() {
    this.gl.deleteTexture(this.texture);
    this.gl.deleteBuffer(this.quadBuf);
    this.gl.deleteVertexArray(this.vao);
  }
}
