// 2D orthographic camera for data-space <-> screen-space mapping
//
//   view:      { x0, y0, x1, y1 }  in data coordinates
//   plotRect:  { ml, mt, pw, ph }  in device pixels
//   canvas:    { W, H }            in device pixels
//
// y0 is BOTTOM of view (math convention: y up)

export class OrthoCamera {
  constructor(view, plotRect, canvas) {
    this._initial = { ...view };
    this.view = { ...view };
    this.plotRect = { ...plotRect };
    this.canvas = { ...canvas };
  }

  reset() {
    this.view = { ...this._initial };
  }

  setPlotRect({ ml, mt, pw, ph }) {
    this.plotRect = { ml, mt, pw, ph };
  }

  setCanvas({ W, H }) {
    this.canvas = { W, H };
  }

  // Data (x, y) -> device pixels
  dataToPix(x, y) {
    const { view, plotRect } = this;
    const u = (x - view.x0) / (view.x1 - view.x0);
    const v = (y - view.y0) / (view.y1 - view.y0);

    return [plotRect.ml + u * plotRect.pw, plotRect.mt + (1 - v) * plotRect.ph];
  }

  // Device pixels -> data (x, y)
  pixToData(px, py) {
    const { view, plotRect } = this;
    const u = (px - plotRect.ml) / plotRect.pw;
    const v = 1 - (py - plotRect.mt) / plotRect.ph;

    return [
      view.x0 + u * (view.x1 - view.x0),
      view.y0 + v * (view.y1 - view.y0),
    ];
  }

  // Data (x, y) -> clip space NDC in [-1, 1]
  dataToNDC(x, y) {
    const [px, py] = this.dataToPix(x, y);

    return [(px / this.canvas.W) * 2 - 1, 1 - (py / this.canvas.H) * 2];
  }

  // Zoom by factor keeping device-pixel anchor (px, py) fixed
  // factor < 1 = zoom in
  zoomAt(px, py, factor) {
    const [ax, ay] = this.pixToData(px, py);
    const v = this.view;
    const w = (v.x1 - v.x0) * factor;
    const h = (v.y1 - v.y0) * factor;
    // Keep anchor at same relative position
    const u = (ax - v.x0) / (v.x1 - v.x0);
    const s = (ay - v.y0) / (v.y1 - v.y0);
    v.x0 = ax - u * w;
    v.x1 = v.x0 + w;
    v.y0 = ay - s * h;
    v.y1 = v.y0 + h;
  }

  // Pan by device-pixel \delta (positive dx moves view right)
  pan(dxPix, dyPix) {
    const v = this.view;
    const sx = (v.x1 - v.x0) / this.plotRect.pw;
    const sy = (v.y1 - v.y0) / this.plotRect.ph;
    v.x0 -= dxPix * sx;
    v.x1 -= dxPix * sx;
    v.y0 += dyPix * sy;
    v.y1 += dyPix * sy;
  }

  // Clamp both axes to stay inside [minW, maxW] x [minH, maxH] and
  // inside given data bounds (if provided)
  clamp({
    minW = 1e-5,
    maxW = 1e9,
    minH = 1e-5,
    maxH = 1e9,
    bounds = null,
  } = {}) {
    const v = this.view;
    let w = v.x1 - v.x0;
    if (w < minW) {
      const c = 0.5 * (v.x0 + v.x1);
      v.x0 = c - minW / 2;
      v.x1 = c + minW / 2;
    } else if (w > maxW) {
      const c = 0.5 * (v.x0 + v.x1);
      v.x0 = c - maxW / 2;
      v.x1 = c + maxW / 2;
    }

    let h = v.y1 - v.y0;
    if (h < minH) {
      const c = 0.5 * (v.y0 + v.y1);
      v.y0 = c - minH / 2;
      v.y1 = c + minH / 2;
    } else if (h > maxH) {
      const c = 0.5 * (v.y0 + v.y1);
      v.y0 = c - maxH / 2;
      v.y1 = c + maxH / 2;
    }
    if (bounds) {
      const w2 = v.x1 - v.x0,
        h2 = v.y1 - v.y0;
      if (v.x0 < bounds.x0) {
        v.x0 = bounds.x0;
        v.x1 = v.x0 + w2;
      }
      if (v.x1 > bounds.x1) {
        v.x1 = bounds.x1;
        v.x0 = v.x1 - w2;
      }
      if (v.y0 < bounds.y0) {
        v.y0 = bounds.y0;
        v.y1 = v.y0 + h2;
      }
      if (v.y1 > bounds.y1) {
        v.y1 = bounds.y1;
        v.y0 = v.y1 - h2;
      }
    }
  }

  // Uniform payload for a shader
  //   uniform vec4 u_view;   // x0, y0, x1, y1
  //   uniform vec4 u_plot;   // ml, mt, pw, ph
  //   uniform vec2 u_canvas; // W, H
  uniforms() {
    const v = this.view;
    const r = this.plotRect;
    const c = this.canvas;

    return {
      u_view: [v.x0, v.y0, v.x1, v.y1],
      u_plot: [r.ml, r.mt, r.pw, r.ph],
      u_canvas: [c.W, c.H],
    };
  }
}

// Derive a plot rect that leaves room for axes
export function defaultPlotRect(W, H) {
  const ml = Math.round(H * 0.075);
  const mr = Math.round(H * 0.02);
  const mt = Math.round(H * 0.03);
  const mb = Math.max(Math.round(H * 0.095), 46);

  return { ml, mt, pw: W - ml - mr, ph: H - mt - mb };
}
