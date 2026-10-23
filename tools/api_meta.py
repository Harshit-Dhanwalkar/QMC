"""Hand-curated metadata for tools/gen_api.py.

Every public header (core/, physics/, export/, latex/) must appear in
SUMMARIES; the generator fails otherwise, so a new module cannot ship without a
one-line description and a category. Edit here, then run `make api-docs`.
"""

CORE_TYPES = "Core: data types & utilities"
CORE_LA = "Core: linear algebra"
CORE_FFT = "Core: FFT"
CORE_ODE = "Core: ODE & time evolution"
CORE_SPECIAL = "Core: special functions"
SP = "Single-particle quantum mechanics"
ATOMS = "Atomic structure & perturbation theory"
OPEN = "Two-level, open systems & quantum optics"
CHEM = "Quantum chemistry (molecules)"
QMCARLO = "Quantum Monte Carlo"
CM = "Condensed matter & many-body"
QC = "Quantum computing & information"
EXPORT = "Export, plotting & LaTeX"

# stem -> (category, one-line summary)
SUMMARIES = {
    # ---- core ----
    "complex": (CORE_TYPES, "Complex number type complex_t and inline arithmetic (c_add, c_mul, c_exp, c_abs, ...)."),
    "constants": (CORE_TYPES, "Fundamental physical constants in SI units (HBAR, M_ELECTRON, E_CHARGE, EPSILON_0, ...)."),
    "matrix": (CORE_TYPES, "Dense complex matrix type cmatrix_t: allocation, copy, multiply, transpose and related operations."),
    "vector": (CORE_TYPES, "Complex vector type cvector_t: allocation, copy, dot product, norm and related operations."),
    "sparse": (CORE_TYPES, "CSR sparse matrices: construction from dense, matrix-vector products (including Hermitian)."),
    "shift_invert": (CORE_LA, "Shift-invert eigensolver for k interior eigenvalues and eigenvectors of a Hermitian sparse matrix."),
    "random": (CORE_TYPES, "xoshiro256** random number generator with uniform and Gaussian draws and jump-ahead parallel streams."),
    "utils": (CORE_TYPES, "Grid and array utilities: linspace, logspace, range, wavefunction normalisation and norms."),
    "qmc_version": (CORE_TYPES, "Library version macros (compile time) and qmc_version()/qmc_version_number() (run time)."),
    "catalog": (CORE_TYPES, "Module catalog and help: discover modules and print where and how to use them (qmc_help)."),
    "linalg": (CORE_LA, "Umbrella header for the linear-algebra routines (LU, QR, SVD, eigensolvers)."),
    "lu": (CORE_LA, "LU decomposition with partial pivoting: solve, determinant and inverse."),
    "qr": (CORE_LA, "QR decomposition with solve and inverse."),
    "svd": (CORE_LA, "Singular value decomposition, SVD solve and pseudoinverse."),
    "tridiag_eigh": (CORE_LA, "Eigenvalues and eigenvectors of real symmetric tridiagonal matrices (implicit QL)."),
    "complex_eigh": (CORE_LA, "Eigen-decomposition of complex Hermitian matrices via the real-embedding trick (cmatrix_eigh_complex)."),
    "eigen_generic": (CORE_LA, "Generic Hermitian eigensolver: LAPACK when built with USE_LAPACK, the built-in solver otherwise."),
    "fft": (CORE_FFT, "1D radix-2 Cooley-Tukey FFT and inverse for power-of-two lengths."),
    "fft2d": (CORE_FFT, "2D FFT by the separable row-column algorithm."),
    "fft3d": (CORE_FFT, "3D FFT by separable row-column-depth decomposition."),
    "fft_wrapper": (CORE_FFT, "FFT wrapper that dispatches to the built-in implementation or to FFTW."),
    "rk4": (CORE_ODE, "Fourth-order Runge-Kutta integrator for complex ODE systems."),
    "numerov": (CORE_ODE, "Numerov integrator for the 1D time-independent Schrodinger equation."),
    "crank_nicolson": (CORE_ODE, "Crank-Nicolson time evolution of the 1D Schrodinger equation, with complex absorbing potentials."),
    "special": (CORE_SPECIAL, "Special functions: Hermite, Laguerre, Legendre, Bessel, spherical harmonics and related polynomials."),
    # ---- export ----
    "csv_writer": (EXPORT, "Streaming CSV writer: configurable delimiter and precision, RFC 4180 escaping, append mode, column descriptors."),
    "json_writer": (EXPORT, "Streaming JSON writer for results and metadata: objects, arrays, quantities with units, schema header."),
    "hdf5_writer": (EXPORT, "HDF5 export: typed datasets, groups, attributes, chunking and compression, appendable datasets (needs USE_HDF5=1)."),
    "netcdf_writer": (EXPORT, "NetCDF export of 1D series and matrices (needs USE_NETCDF=1)."),
    "vtk_writer": (EXPORT, "Legacy VTK writer for scalar fields on structured grids (open in ParaView or VisIt)."),
    "plot": (EXPORT, "Plotting front end (lines, multi-series) over GR, gnuplot or matplotlib back ends, chosen with PLOT_BACKEND."),
    "export_schema": (EXPORT, "Conventions shared by all exporters: schema version, library name and version, UTC timestamps, column descriptors."),
    "latex_gen": (EXPORT, "LaTeX generation: escaping, tables, matrices and documents, rendered to PDF or PNG when pdflatex is available."),
    # ---- single-particle QM ----
    "wavefn": (SP, "Wavefunction container on a 1D grid: allocation, normalisation and probability density."),
    "potentials": (SP, "Standard 1D potentials (infinite and finite well, harmonic, step, barrier, Morse, Coulomb, Yukawa) as V(x) functions."),
    "schrodinger": (SP, "1D Schrodinger solvers: matrix diagonalisation, shooting, Crank-Nicolson and split-step time evolution."),
    "uncertainty": (SP, "Position and momentum uncertainties and energy expectation values for wavefunctions."),
    "central_potential": (SP, "General 3D central-potential radial solver."),
    "hydrogen": (SP, "Hydrogen atom: radial Schrodinger solver, energy levels and radial wavefunctions."),
    "orbital_sample": (SP, "Monte Carlo sampling of hydrogen orbitals psi_nlm, for visualisation."),
    "box_revival": (SP, "Particle in a box: Gaussian packet evolution, wavefunction revivals and density maps."),
    "circular_billiard": (SP, "Quantum billiard in a disk: Bessel functions and zeros, modes and Weyl counting."),
    "ab_ring": (SP, "Aharonov-Bohm ring: energy levels versus flux, persistent current and flux-dependent transmission."),
    "double_barrier": (SP, "Resonant tunnelling: transfer-matrix transmission through stacked square barriers, wavefunction profile and dwell time."),
    "wkb": (SP, "WKB approximation: Bohr-Sommerfeld quantisation and tunnelling transmission."),
    "scattering": (SP, "Scattering theory: phase shifts, Born approximation, partial waves and the S-matrix."),
    "soft": (SP, "Split-operator Fourier method for 2D and 3D time-dependent Schrodinger evolution."),
    "relativistic": (SP, "Relativistic wave equations: 1D Klein-Gordon and Dirac solvers and the radial Dirac equation."),
    "dirac_evolve": (SP, "1D Dirac equation time evolution of a two-component spinor (split-operator Fourier method)."),
    # ---- atoms & perturbation theory ----
    "helium": (ATOMS, "Two-electron atom or ion ground-state energy by the effective-nuclear-charge variational method."),
    "angular": (ATOMS, "Angular momentum and spin: Pauli matrices, L+/L-/L^2/Lz, Clebsch-Gordan and Gaunt coefficients."),
    "fine_structure": (ATOMS, "Hydrogen fine structure: spin-orbit coupling and the closed-form total fine-structure formula."),
    "zeeman": (ATOMS, "Weak-field anomalous Zeeman effect for a single valence electron."),
    "perturbation": (ATOMS, "Time-independent perturbation theory (non-degenerate and degenerate) and Fermi's golden rule."),
    "td_perturbation": (ATOMS, "Time-dependent perturbation theory: Dyson series to arbitrary order and first-order closed forms."),
    "variational": (ATOMS, "Variational method: generic golden-section minimiser and variational energy minimisation."),
    "identical": (ATOMS, "Identical particles: Slater determinants for fermions and permanents for bosons."),
    "hartree_fock": (ATOMS, "Restricted Hartree-Fock SCF for closed-shell atoms and ions."),
    "dft": (ATOMS, "Kohn-Sham LDA density functional theory for closed-shell atoms and ions."),
    # ---- two-level / open systems ----
    "rabi": (OPEN, "Rabi oscillations of a two-level system in the rotating-wave approximation, exact and numerical."),
    "driven": (OPEN, "Driven two-level systems: general time-dependent Hamiltonians, Landau-Zener sweeps, lab-frame driving beyond the RWA."),
    "landau_zener": (OPEN, "Landau-Zener sweeps of a two-level system on the Bloch sphere."),
    "floquet": (OPEN, "Floquet theory for periodically driven systems: quasienergies from time-domain and Sambe-space solvers."),
    "bloch": (OPEN, "Single-qubit Bloch-vector picture on top of the Lindblad solver."),
    "lindblad": (OPEN, "Lindblad (GKSL) master equation for open quantum systems."),
    "zeno": (OPEN, "Quantum Zeno effect: pulsed and continuous measurement survival probabilities and Monte Carlo trajectories."),
    "kicked_rotor": (OPEN, "Quantum kicked rotor and the classical standard map."),
    "phase_space": (OPEN, "Phase-space picture of a bosonic mode: Fock states, Kerr nonlinearity and the Wigner function."),
    # ---- quantum chemistry ----
    "basis_parser": (CHEM, "Parser for Gaussian94-format basis-set files (as exported by Basis Set Exchange)."),
    "molecular_integrals": (CHEM, "Gaussian-type-orbital molecular integrals by McMurchie-Davidson, with built-in STO-3G data."),
    "molecular_hf": (CHEM, "Restricted Hartree-Fock for molecules in a general basis, plus AO-to-MO integral transforms."),
    "mp2": (CHEM, "Second-order Moller-Plesset perturbation theory for closed-shell molecules."),
    "cisd": (CHEM, "Configuration interaction singles and doubles (optionally with a frozen core)."),
    "fci": (CHEM, "Full configuration interaction: exact diagonalisation of the molecular electronic Hamiltonian."),
    "ccsd": (CHEM, "Coupled cluster singles and doubles (CCSD) for closed-shell molecules."),
    "ccsd_t": (CHEM, "Perturbative triples correction CCSD(T) on top of CCSD."),
    "casscf": (CHEM, "Complete-active-space SCF: joint optimisation of CI coefficients and orbitals."),
    "molecular_dft": (CHEM, "Molecular Kohn-Sham LDA DFT on a Becke fuzzy-Voronoi integration grid."),
    "tddft": (CHEM, "Linear-response TDDFT excitation energies in the Tamm-Dancoff approximation."),
    "geometry_optimization": (CHEM, "Generic Cartesian geometry optimiser driven by an energy-and-gradient callback."),
    "vibrational": (CHEM, "Harmonic vibrational frequencies from numerical Hessians of analytic gradients."),
    # ---- QMC ----
    "vmc": (QMCARLO, "Variational Monte Carlo for two-electron atoms and ions with a Slater-Jastrow ansatz."),
    "dmc": (QMCARLO, "Diffusion Monte Carlo for two-electron atoms and ions with importance sampling."),
    "pimc": (QMCARLO, "Path integral Monte Carlo for two-electron atoms at finite temperature."),
    # ---- condensed matter / many-body ----
    "anderson": (CM, "Anderson localisation on a 1D tight-binding chain with random on-site disorder."),
    "lattice": (CM, "Tight-binding lattice models on 1D and 2D lattices, including Hofstadter and Chern-number routines."),
    "tight_binding": (CM, "Generic multi-band tight-binding models: bands, DOS, Berry curvature, Chern numbers, Zak and Wilson loops."),
    "ssh_chain": (CM, "Su-Schrieffer-Heeger chain, the simplest topological insulator."),
    "ising_chain": (CM, "Transverse-field Ising chain: exact diagonalisation on a ring and the exact Jordan-Wigner solution."),
    "tfim_quench": (CM, "Exact quench dynamics of the transverse-field Ising chain via free fermions."),
    "spin_chain": (CM, "Spin-1/2 XXZ chain on a ring: symmetry-adapted exact diagonalisation and Lanczos."),
    "hubbard_dimer": (CM, "Two-site Hubbard model: exact diagonalisation, eigenstates and time evolution."),
    "dmrg": (CM, "Infinite-system density matrix renormalisation group for spin chains."),
    "finite_dmrg": (CM, "Finite-system DMRG sweeps."),
    "quantum_walk": (CM, "Discrete-time quantum walk on a line with a two-state coin."),
    "second_quant": (CM, "Second quantisation and the Jordan-Wigner transform: fermionic operators and molecular Hamiltonians in qubit form."),
    # ---- quantum computing & information ----
    "qubits": (QC, "Minimal multi-qubit state-vector simulator: single-qubit gates and CNOT."),
    "qft": (QC, "Quantum Fourier transform and its inverse on a chosen set of qubits."),
    "qpe": (QC, "Quantum phase estimation: circuit simulation, phase estimate and output distribution."),
    "quantum_algorithms": (QC, "Deutsch-Jozsa and Grover search simulations."),
    "grover_multi": (QC, "Grover search with several marked items: state evolution and success probability."),
    "quantum_info": (QC, "Quantum information protocols: teleportation, superdense coding and CHSH."),
    "bell_chsh": (QC, "Bell and CHSH correlations: quantum values and local-hidden-variable and deterministic bounds."),
    "boson_sampling": (QC, "Boson sampling: amplitudes and probabilities for N photons through an interferometer."),
    "qec": (QC, "Three-qubit bit-flip and phase-flip error-correcting codes with ancilla syndrome extraction."),
    "qec5": (QC, "The five-qubit perfect [[5,1,3]] error-correcting code."),
    "qec_surface17": (QC, "Distance-3 rotated surface code ([[9,1,3]])."),
    "vqe": (QC, "Variational quantum eigensolver with a hardware-efficient ansatz."),
    "vqe_noisy": (QC, "Noisy VQE: amplitude-damping and dephasing noise after every gate."),
    "qaoa": (QC, "QAOA for MaxCut and Ising cost functions: simulator, closed-form p=1 expectation and parameter optimisation."),
    "ucc": (QC, "Unitary coupled-cluster ansatz built from converged CCSD amplitudes."),
    "physics": (CORE_TYPES, "Umbrella header including the core headers and every physics module."),
}

# stem -> explicit "start here" functions (default: first declared functions)
START = {
    "csv_writer": ["csv_options_default", "csv_open", "csv_write_row", "csv_close"],
    "hdf5_writer": ["hdf5_open", "hdf5_write_dataset", "hdf5_close"],
    "json_writer": ["json_open", "json_begin_object", "json_write_field", "json_close"],
    "netcdf_writer": ["netcdf_open", "netcdf_writer_write_1d", "netcdf_close"],
    "latex_gen": ["latex_tools_available", "latex_render_to_pdf", "latex_document_create"],
    "plot": ["plot_line", "plot_lines", "plot_backend_name"],
    "angular": ["lz_matrix", "l2_matrix", "clebsch_gordan"],
    "potentials": ["V_harmonic", "V_finite_well", "V_barrier"],
    "schrodinger": ["solve_tise_matrix", "solve_tise_shoot", "evolve_tdse_crank"],
    "wavefn": ["wavefunction_alloc", "wavefunction_normalize", "wavefunction_free"],
    "special": ["hermite", "laguerre", "legendre"],
}

# stem -> existing hand-written guide page (relative to docs/src), when the
# name does not match docs/src/physics/<stem>.md
GUIDES = {
    "angular": "physics/angular_momentum.md",
    "wavefn": "physics/wavefunctions.md",
    "potentials": "physics/potentials_1d.md",
    "molecular_integrals": "physics/moecular_integrals.md",
    "crank_nicolson": "internals/crank_nicolson.md",
    "numerov": "internals/numerov.md",
    "special": "internals/special.md",
    "sparse": "internals/sparse_solvers.md",
    "shift_invert": "internals/sparse_solvers.md",
    "plot": "internals/plotting.md",
    "fft": "internals/fft.md",
    "fft2d": "internals/fft.md",
    "fft3d": "internals/fft.md",
    "fft_wrapper": "internals/fft.md",
    "linalg": "internals/linalg.md",
    "lu": "internals/linalg.md",
    "qr": "internals/linalg.md",
    "svd": "internals/linalg.md",
    "tridiag_eigh": "internals/linalg.md",
    "complex_eigh": "internals/linalg.md",
    "eigen_generic": "internals/linalg.md",
    "infinite_well": None,
}
