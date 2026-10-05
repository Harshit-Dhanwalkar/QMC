#!/usr/bin/env python3
"""Build docs/src/playground/preview.gif from frames dumped by dump_frames.mjs.

    node wasm/tools/dump_frames.mjs /tmp/qmc_frames.bin
    python3 wasm/tools/make_gif.py /tmp/qmc_frames.bin docs/src/playground/preview.gif

Needs: pip install pillow numpy
"""
import sys
import numpy as np
from PIL import Image, ImageDraw

N, SX, SCALE, SIDE = 128, 109, 3, 96   # SX = detector column (x = +14 on a 40-wide box)
src, dst = sys.argv[1], sys.argv[2]
raw = np.fromfile(src, dtype=np.uint8)
V = np.frombuffer(raw[: N * N * 8].tobytes(), dtype=np.float64).reshape(N, N)
frames = np.frombuffer(raw[N * N * 8 :].tobytes(), dtype=np.float32).reshape(-1, N, N)

stops = np.array([[0,0,4],[40,11,84],[101,21,110],[159,42,99],[212,72,66],[245,125,21],[250,193,39],[252,255,164]], float)
xs = np.linspace(0, len(stops) - 1, 256)
lut = np.stack([np.interp(xs, np.arange(len(stops)), stops[:, c]) for c in range(3)], axis=1)

wall = np.where(V > 20, 0.8, np.where(V > 0.5, 0.45, 0.0))[..., None]
acc = np.zeros(N)
vmax = 0.02
out = []
for d in frames:
    vmax = max(vmax * 0.992, float(d.max()), 0.003)
    idx = np.minimum(255, (np.sqrt(d / vmax) * 255).astype(int))
    rgb = lut[idx]
    rgb = rgb * (1 - wall) + np.array([130, 160, 210]) * wall
    img = Image.fromarray(rgb.transpose(1, 0, 2).astype(np.uint8)).resize((N * SCALE, N * SCALE), Image.BICUBIC)
    acc += d[SX]
    side = Image.new("RGB", (SIDE, N * SCALE), (5, 7, 11))
    dr = ImageDraw.Draw(side)
    m = acc.max() or 1.0

    for iy in range(N):
        dr.rectangle([0, iy * SCALE, int(acc[iy] / m * (SIDE - 6)), iy * SCALE + SCALE - 1], fill=(61, 220, 151))

    dr2 = ImageDraw.Draw(img)
    x = int((SX + 0.5) * SCALE)

    for y in range(0, N * SCALE, 9):
        dr2.line([x, y, x, y + 4], fill=(61, 220, 151))

    canvas = Image.new("RGB", (N * SCALE + 6 + SIDE, N * SCALE), (11, 14, 20))
    canvas.paste(img, (0, 0)); canvas.paste(side, (N * SCALE + 6, 0))
    out.append(canvas.quantize(colors=96, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE))

out[0].save(dst, save_all=True, append_images=out[1:], duration=40, loop=0, optimize=True)
print(f"wrote {dst}: {len(out)} frames")
