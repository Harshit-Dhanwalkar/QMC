(() => {
  const DEMOS = [
    ["index", "Double slit", "2D wave interference, split-operator FFT"],
    [
      "butterfly",
      "Hofstadter butterfly",
      "Electrons on a lattice in a magnetic field",
    ],
    ["dirac", "Klein paradox", "Relativistic wavepacket on a potential step"],
    ["orbitals", "Hydrogen orbitals", "Rotatable 3D point clouds"],
    ["bloch", "Bloch sphere", "A driven, damped qubit"],
    ["anderson", "Anderson localisation", "Why disorder stops a wave"],
    ["quantum_walk", "Quantum walk", "Linear spreading instead of √t"],
    ["landau_zener", "Landau-Zener", "Jumping across an avoided crossing"],
    ["kicked_rotor", "Kicked rotor", "Classical chaos, quantum freeze"],
    ["ssh_chain", "SSH chain", "Topological edge states in one dimension"],
    ["ising_quench", "Ising quench", "Light cone and dynamical phase transitions"],
  ];
  const GH =
    '<svg viewBox="0 0 16 16" width="14" height="14" fill="currentColor" aria-hidden="true" style="vertical-align:-2px;margin-right:5px"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82a7.6 7.6 0 0 1 4 0c1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.01 8.01 0 0 0 16 8c0-4.42-3.58-8-8-8z"/></svg>';

  const css = `
  .qnav { display:flex; flex-wrap:wrap; gap:8px; align-items:center; }
  .qnav a, .qnav summary { margin:0 !important; color:#d7dce5; text-decoration:none; border:1px solid #232b3d; padding:5px 11px; border-radius:8px;
    font-size:.9rem; display:inline-block; cursor:pointer; background:transparent; list-style:none; line-height:1.5; }
  .qnav a:hover, .qnav summary:hover, .qnav details[open] > summary { border-color:#3ddc97; color:#3ddc97; }
  .qnav summary::-webkit-details-marker { display:none; }
  .qnav summary::after { content:" \\25BE"; font-size:.8em; }
  .qnav details { position:relative; }
  .qnav .menu { position:absolute; left:0; top:calc(100% + 6px); z-index:50; width:min(330px, 88vw); padding:6px; border-radius:12px;
    background:#131824; border:1px solid #232b3d; box-shadow:0 14px 40px rgba(0,0,0,.55); }
  @media (min-width:700px) { .qnav .menu { left:auto; right:0; } }
  .qnav .menu a { display:block; border:0; border-radius:8px; padding:7px 10px; color:#d7dce5; }
  .qnav .menu a small { display:block; color:#8993a8; font-size:.78rem; margin-top:1px; }
  .qnav .menu a:hover { background:#1b2234; color:#3ddc97; }
  .qnav .menu a.cur { background:#10241d; color:#3ddc97; }
  .qnav .menu a.cur small { color:#6fae95; }
  .qnav a:focus-visible, .qnav summary:focus-visible { outline:2px solid #3ddc97; outline-offset:2px; }
  `;

  const style = document.createElement("style");
  style.textContent = css;
  document.head.appendChild(style);
  const here =
    (location.pathname.split("/").pop() || "index.html").replace(
      /\.html$/,
      "",
    ) || "index";
  const mount = document.querySelector(".links");
  if (!mount) {
    return;
  }
  mount.classList.add("qnav");
  mount.removeAttribute("style");
  const items = DEMOS.map(
    ([id, name, blurb]) =>
      `<a href="${id}.html"${id === here ? ' class="cur" aria-current="page"' : ""}>${name}<small>${blurb}</small></a>`,
  ).join("");

  const cur = DEMOS.find((d) => d[0] === here);
  mount.innerHTML =
    `<details><summary>${cur ? "More demos" : "Demos"}</summary><div class="menu">${items}</div></details>` +
    `<a href="https://github.com/Harshit-Dhanwalkar/QMC">${GH}GitHub</a>` +
    `<a href="../showcase.html">Showcase</a><a href="../introduction.html">Docs</a>`;
  document.addEventListener("click", (e) => {
    const d = mount.querySelector("details");
    if (d && !d.contains(e.target)) d.open = false;
  });
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") {
      const d = mount.querySelector("details");
      if (d) d.open = false;
    }
  });
})();
