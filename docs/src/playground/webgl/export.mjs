// PNG export: renders GL content and renders overlay text and
// exportPNG composites them into a single image and downloads it

export async function exportPNG({
  width,
  height,
  renderGL,
  renderOverlay,
  filename = "export.png",
}) {
  // Offscreen GL canvas
  const glCanvas = document.createElement("canvas");
  glCanvas.width = width;
  glCanvas.height = height;

  await renderGL(glCanvas, width, height);

  // Composite canvas (2D) for final image
  const out = document.createElement("canvas");
  out.width = width;
  out.height = height;
  const ctx = out.getContext("2d");
  ctx.drawImage(glCanvas, 0, 0);

  // Browsers cap live WebGL contexts (~16); release this one now rather than
  // waiting for GC, or repeated exports can evict page's main context
  glCanvas
    .getContext("webgl2")
    ?.getExtension("WEBGL_lose_context")
    ?.loseContext();

  if (renderOverlay) {
    renderOverlay(ctx, width, height);
  }

  return new Promise((resolve) => {
    out.toBlob((blob) => {
      if (!blob) {
        resolve(null);

        return;
      }

      const a = document.createElement("a");
      a.href = URL.createObjectURL(blob);
      a.download = filename;
      document.body.appendChild(a); // Safari/old Firefox need it in DOM
      a.click();
      a.remove();

      setTimeout(() => URL.revokeObjectURL(a.href), 2000);

      resolve(blob);
    });
  });
}
