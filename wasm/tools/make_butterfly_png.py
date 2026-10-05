#!/usr/bin/env python3
"""Render the Hofstadter butterfly (gaps coloured by Chern number) from JSON
written by dump_butterfly.mjs.

    node wasm/tools/dump_butterfly.mjs 50 /tmp/butterfly.json
    python3 wasm/tools/make_butterfly_png.py /tmp/butterfly.json \
        docs/src/playground/butterfly.png docs/src/playground/butterfly_grow.gif

Needs: pip install pillow
"""
import json
import sys
from PIL import Image, ImageDraw, ImageFont

src, png_out = sys.argv[1], sys.argv[2]
gif_out = sys.argv[3] if len(sys.argv) > 3 else None
rows = json.load(open(src))

BG = (11, 14, 20)
BAND = (22, 28, 42)
E_MAX = 4.15


def font(size):
    for path in ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
                 "/usr/share/fonts/dejavu/DejaVuSans.ttf",
                 "/usr/share/fonts/TTF/DejaVuSans.ttf",
                 "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"):
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            continue
    return ImageFont.load_default()


def gap_color(t):
    if t is None:
        return (255, 255, 255)
    if t == 0:
        return (96, 110, 138)
    s = min(abs(t), 8) / 8.0
    if t > 0:   # warm: pale yellow -> red
        a, b = (255, 214, 120), (232, 40, 62)
    else:       # cool: pale cyan -> blue
        a, b = (130, 228, 255), (48, 84, 255)
    return tuple(int(a[i] + (b[i] - a[i]) * s) for i in range(3))


def render(qmax, W, H, margin_l=70, margin_b=72, margin_t=18, margin_r=18, labels=True):
    img = Image.new("RGB", (W, H), BG)
    d = ImageDraw.Draw(img)
    pw, ph = W - margin_l - margin_r, H - margin_t - margin_b

    def X(alpha):
        return margin_l + alpha * pw

    def Y(e):
        return margin_t + (E_MAX - e) / (2 * E_MAX) * ph

    # each flux paints the horizontal span nearest to it (nearest-neighbour
    # fill), so low-denominator fluxes do not leave empty stripes
    sel = sorted((r for r in rows if r["q"] <= qmax), key=lambda r: r["p"] / r["q"])
    alphas = [r["p"] / r["q"] for r in sel]
    for i, row in enumerate(sel):
        p, q, ed, tk = row["p"], row["q"], row["edges"], row["tknn"]
        left = 0.0 if i == 0 else 0.5 * (alphas[i - 1] + alphas[i])
        right = 1.0 if i == len(sel) - 1 else 0.5 * (alphas[i] + alphas[i + 1])
        x0, x1 = int(X(left)), max(int(X(left)), int(X(right)) - 1)
        for b in range(q):
            d.rectangle([x0, Y(ed[2 * b + 1]), x1, Y(ed[2 * b])], fill=BAND)
        for r in range(1, q):
            lo, hi = ed[2 * (r - 1) + 1], ed[2 * r]
            if hi - lo > 1e-9:
                d.rectangle([x0, Y(hi), x1, Y(lo)], fill=gap_color(tk[r - 1]))
    if labels:
        f = font(max(11, H // 52))
        for a, txt in [(0, "0"), (0.25, "1/4"), (1 / 3, "1/3"), (0.5, "1/2"), (2 / 3, "2/3"), (0.75, "3/4"), (1, "1")]:
            d.line([X(a), Y(-E_MAX), X(a), Y(-E_MAX) + 5], fill=(137, 147, 168))
            d.text((X(a), Y(-E_MAX) + 8), txt, fill=(137, 147, 168), font=f, anchor="mt")
        for e in (-4, -2, 0, 2, 4):
            d.line([margin_l - 5, Y(e), margin_l, Y(e)], fill=(137, 147, 168))
            d.text((margin_l - 9, Y(e)), str(e), fill=(137, 147, 168), font=f, anchor="rm")
        d.text((margin_l + pw / 2, H - 10), "magnetic flux per plaquette  α = p/q", fill=(215, 220, 229), font=f, anchor="ms")
        d.text((margin_l + 6, margin_t + 4), "E / t", fill=(215, 220, 229), font=f, anchor="la")
        d.text((W - margin_r - 6, margin_t + 4), f"q ≤ {qmax}", fill=(215, 220, 229), font=f, anchor="ra")
    return img


qmax_all = max(r["q"] for r in rows)
big = render(qmax_all, 2400, 1500)
big.resize((1600, 1000), Image.LANCZOS).save(png_out, optimize=True)
print("wrote", png_out)

if gif_out:
    frames = []
    for qm in list(range(3, 21)) + list(range(22, qmax_all + 1, 2)):
        frames.append(render(qm, 900, 560).quantize(colors=128, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE))
    last = render(qmax_all, 900, 560).quantize(colors=128, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    frames += [last] * 20
    durations = [120] * (len(frames) - 20) + [100] * 20
    frames[0].save(gif_out, save_all=True, append_images=frames[1:], duration=durations, loop=0, optimize=True)
    print("wrote", gif_out, len(frames), "frames")
