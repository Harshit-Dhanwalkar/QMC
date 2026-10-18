// 2D-canvas overlay for axes, ticks, crosshairs and text on top of
// WebGL canvas
// Uses a separate <canvas> element; no WebGL state here
// with High-DPI (Retina) support

export class Overlay2D {
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx = canvas.getContext("2d");
  }

  setSize(W, H) {
    if (this.canvas.width !== W || this.canvas.height !== H) {
      this.canvas.width = W;
      this.canvas.height = H;
      this.ctx.resetTransform();
    }
  }

  clear() {
    this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
  }

  // Draws ticks and labels around plot area
  //
  //   camera:  OrthoCamera
  //   xTicks:  [{ value, label }]  in data coords
  //   yTicks:  same, on y axis
  //   opts:    { font, color, tickLen, xLabel, yLabel, topRight }
  drawAxes(camera, xTicks, yTicks, opts = {}) {
    const c = this.ctx;
    const W = this.canvas.width;
    const H = this.canvas.height;
    const { ml, mt, pw, ph } = camera.plotRect;
    const v = camera.view;

    const fs = opts.font ?? Math.max(11, Math.round(H * 0.026));
    const fam =
      opts.family ?? "system-ui,-apple-system,Segoe UI,Roboto,sans-serif";
    const col = opts.color ?? "#8993a8";
    const tickLen = opts.tickLen ?? 5;

    c.font = `${fs}px ${fam}`;
    c.fillStyle = col;
    c.strokeStyle = col;
    c.lineWidth = 1;

    // X axis (bottom)
    c.textAlign = "center";
    c.textBaseline = "top";
    for (const t of xTicks) {
      if (t.value < v.x0 || t.value > v.x1) {
        continue;
      }

      const [px] = camera.dataToPix(t.value, v.y0);
      c.beginPath();
      c.moveTo(px + 0.5, mt + ph);
      c.lineTo(px + 0.5, mt + ph - tickLen);
      c.stroke();
      c.fillText(t.label, px, mt + ph + 3);
    }

    // Y axis (left)
    c.textAlign = "right";
    c.textBaseline = "middle";
    for (const t of yTicks) {
      if (t.value < v.y0 || t.value > v.y1) {
        continue;
      }

      const [, py] = camera.dataToPix(v.x0, t.value);
      c.beginPath();
      c.moveTo(ml - tickLen, py + 0.5);
      c.lineTo(ml, py + 0.5);
      c.stroke();
      c.fillText(t.label, ml - tickLen - 4, py);
    }

    // Axis titles
    if (opts.xLabel) {
      c.fillStyle = "#d7dce5";
      c.textAlign = "center";
      c.textBaseline = "alphabetic";
      c.fillText(opts.xLabel, ml + pw / 2, H - 6);
    }
    if (opts.yLabel) {
      c.fillStyle = "#d7dce5";
      c.textAlign = "left";
      c.textBaseline = "top";
      c.fillText(opts.yLabel, ml + 6, mt + 4);
    }
    if (opts.topRight) {
      c.fillStyle = "#d7dce5";
      c.textAlign = "right";
      c.textBaseline = "top";
      c.fillText(opts.topRight, ml + pw - 6, mt + 4);
    }
  }

  // Crosshair drawn at data coords (x, y); solid or dashed
  drawCrosshair(camera, x, y, { solid = false, radius = 4 } = {}) {
    const c = this.ctx;
    const { ml, mt, pw, ph } = camera.plotRect;
    const [px, py] = camera.dataToPix(x, y);

    c.strokeStyle = solid ? "rgba(255,255,255,.95)" : "rgba(255,255,255,.55)";
    c.lineWidth = Math.max(
      1,
      this.canvas.height / window.devicePixelRatio / 500,
    );
    c.setLineDash(solid ? [] : [5, 4]);
    c.beginPath();
    c.moveTo(px, mt);
    c.lineTo(px, mt + ph);
    c.moveTo(ml, py);
    c.lineTo(ml + pw, py);
    c.stroke();
    c.setLineDash([]);
    c.beginPath();
    c.arc(px, py, radius, 0, 2 * Math.PI);
    c.stroke();
  }

  // Free text, in device pixels
  drawText(str, x, y, opts = {}) {
    const c = this.ctx;

    c.font = `${opts.size ?? 13}px ${opts.family ?? "system-ui,sans-serif"}`;
    c.fillStyle = opts.color ?? "#d7dce5";
    c.textAlign = opts.align ?? "left";
    c.textBaseline = opts.baseline ?? "top";
    c.fillText(str, x, y);
  }
}

// Common tick helpers
export function uniformTicks(min, max, n, fmt = (x) => x.toFixed(2)) {
  const out = [];
  for (let i = 0; i <= n; i++) {
    const v = min + (max - min) * (i / n);
    out.push({ value: v, label: fmt(v) });
  }

  return out;
}

export function customTicks(pairs) {
  return pairs.map(([value, label]) => ({ value, label: String(label) }));
}
