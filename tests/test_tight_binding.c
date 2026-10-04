/*
 * Test: generic tight-binding band structure and band topology
 *
 * 1. Dispersions of chain, square lattice and graphene against closed
 *    forms; Dirac points of graphene; band-line and DOS helpers
 * 2. Haldane model: valley gaps against 2|M -+ 3 \sqrt(3) t2 \sin(\phi)|,
 *    Chern-number phase diagram, a trivial massive-graphene limit
 * 3. SSH chain: Zak phase 0 / \pi, cross-checked against edge states of
 *    open chain from lattice.h
 * 4. Kagome lattice: flat band, sum rule, Dirac point
 * 5. Hofstadter: Bloch matrix identical to lattice_hofstadter_bloch(), band
 *    Chern numbers equal to lattice_hofstadter_chern_numbers(), band-group
 *    handling for touching bands of flux 1/4
 * 6. Wilson-loop winding equals Chern number
 * 7. Invalid-input handling
 */

#include "../physics/lattice.h"
#include "../physics/tight_binding.h"
#include "complex.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.10g expected=%.10g err=%.2e\n", label, got, expected,
         err);
  if (!(err <= tol)) {
    printf("  FAIL: %s\n", label);
    failures++;
  }
}

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAIL");
  if (!cond) {
    failures++;
  }
}

static double chern(const tb_model_t *m, int first, int count, int n_k) {
  double c = NAN;
  if (tb_chern_number(m, first, count, n_k, &c) != TB_OK) {
    return NAN;
  }

  return c;
}

static void test_dispersions(void) {
  printf("  === Dispersions, DOS, band lines ===\n");

  tb_model_t *chain = tb_model_chain(0.7);
  tb_model_t *sq = tb_model_square(-0.9);
  tb_model_t *gr = tb_model_graphene(1.3);
  check_true(chain && sq && gr, "models built");
  if (!chain || !sq || !gr) {
    tb_model_free(chain);
    tb_model_free(sq);
    tb_model_free(gr);

    failures++;

    return;
  }

  check_true(tb_model_dim(chain) == 1 && tb_model_orbitals(chain) == 1,
             "chain dim/orbitals");
  check_true(tb_model_dim(gr) == 2 && tb_model_orbitals(gr) == 2,
             "graphene dim/orbitals");

  double worst_chain = 0.0;
  double worst_sq = 0.0;
  double worst_gr = 0.0;
  for (int i = 0; i < 17; i++) {
    double k = 2.0 * M_PI * i / 17.0 + 0.123;
    double e1 = 0.0;
    double e2[2];

    tb_bands(chain, k, 0.0, &e1);
    worst_chain = fmax(worst_chain, fabs(e1 - tb_chain_energy(0.7, k)));
    for (int j = 0; j < 13; j++) {
      double q = 2.0 * M_PI * j / 13.0 + 0.31;
      double es = 0.0;

      tb_bands(sq, k, q, &es);
      worst_sq = fmax(worst_sq, fabs(es - tb_square_energy(-0.9, k, q)));

      tb_bands(gr, k, q, e2);
      double eg = tb_graphene_energy(1.3, k, q);
      worst_gr = fmax(worst_gr, fabs(e2[1] - eg) + fabs(e2[0] + eg));
    }
  }

  check_close(worst_chain, 0.0, 1e-12, "chain E(k) vs 2t \\cos(k)");
  check_close(worst_sq, 0.0, 1e-12, "square E(k) vs closed form");
  check_close(worst_gr, 0.0, 1e-12, "graphene +-|f(k)| vs closed form");

  double e[2];
  tb_bands(gr, 2.0 * M_PI / 3.0, -2.0 * M_PI / 3.0, e);
  check_close(e[1] - e[0], 0.0, 1e-12, "Dirac gap at K");
  tb_bands(gr, -2.0 * M_PI / 3.0, 2.0 * M_PI / 3.0, e);
  check_close(e[1] - e[0], 0.0, 1e-12, "Dirac gap at K'");
  tb_bands(gr, 0.0, 0.0, e);
  check_close(e[1], 3.0 * 1.3, 1e-12, "band top 3t at \\Gamma");

  // Hermiticity of Bloch matrix
  cmatrix_t *h = tb_hamiltonian(gr, 0.4, 1.1);
  check_true(h != NULL, "tb_hamiltonian returns matrix");
  if (h) {
    double herm = 0.0;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        complex_t a = h->data[i * 2 + j];
        complex_t b = c_conj(h->data[j * 2 + i]);

        herm = fmax(herm, hypot(a.re - b.re, a.im - b.im));
      }
    }

    check_close(herm, 0.0, 1e-14, "H(k) Hermitian");

    cmatrix_free(h);
  }

  // Band line Gamma -> K for graphene
  enum { NP = 9 };
  double line[NP * 2];
  check_true(tb_bands_along_line(gr, 0.0, 0.0, 2.0 * M_PI / 3.0,
                                 -2.0 * M_PI / 3.0, NP, line) == TB_OK,
             "bands_along_line OK");
  check_close(line[(NP - 1) * 2 + 1] - line[(NP - 1) * 2], 0.0, 1e-12,
              "line ends at Dirac point");
  check_close(line[1], 3.0 * 1.3, 1e-12, "line starts at band top");

  // DOS normalization, symmetry and van Hove peak
  enum { NE = 801 };
  static double dos[NE];
  double emin = -5.0;
  double emax = 5.0;

  check_true(tb_dos(gr, 96, emin, emax, NE, 0.05, dos) == TB_OK, "DOS OK");

  double de = (emax - emin) / (NE - 1);
  double total = 0.0;
  double asym = 0.0;
  for (int i = 0; i < NE; i++) {
    total += dos[i] * de;
    asym = fmax(asym, fabs(dos[i] - dos[NE - 1 - i]));
  }

  check_close(total, 2.0, 1e-3, "graphene DOS integrates to n_orb");
  check_close(asym, 0.0, 1e-9, "graphene DOS particle-hole symmetric");

  int imax = 0;
  for (int i = 0; i < NE; i++) {
    if (emin + i * de > 0.0 && dos[i] > dos[imax]) {
      imax = i;
    }
  }

  check_close(emin + imax * de, 1.3, 0.1, "van Hove peak at |E| = t");
  check_true(dos[NE / 2] < 0.5 * dos[imax], "DOS dips at Dirac point");

  tb_model_free(chain);
  tb_model_free(sq);
  tb_model_free(gr);
}

static void test_haldane(void) {
  printf("  === Haldane model ===\n");
  const double t1 = 1.0;
  const double t2 = 0.1;
  const double phis[3] = {M_PI / 2.0, -M_PI / 2.0, 0.6};
  const double masses[3] = {0.0, 0.3, 0.7};

  // Valley gaps
  double worst = 0.0;
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      tb_model_t *m = tb_model_haldane(t1, t2, phis[a], masses[b]);
      if (!m) {
        failures++;

        continue;
      }

      for (int v = -1; v <= 1; v += 2) {
        double k = v * 2.0 * M_PI / 3.0;
        double e[2];
        tb_bands(m, k, -k, e);

        double expect = tb_haldane_valley_gap(t2, phis[a], masses[b], v);
        worst = fmax(worst, fabs((e[1] - e[0]) - expect));
      }

      tb_model_free(m);
    }
  }

  check_close(worst, 0.0, 1e-10,
              "valley gaps = 2|M -+ 3 \\sqrt(3) t2 \\sin(\\phi)|");
  check_true(isnan(tb_haldane_valley_gap(t2, 0.5, 0.1, 0)),
             "valley_gap NaN for bad valley");

  // Phase diagram: threshold 3 \sqrt(3) t2 \sin(\phi) = 0.5196 at \phi = \pi/2
  struct {
    double phi;
    double t2;
    double mass;
  } cases[] = {
      {M_PI / 2.0, 0.1, 0.0},  {M_PI / 2.0, 0.1, 0.3},  {M_PI / 2.0, 0.1, -0.3},
      {M_PI / 2.0, 0.1, 0.7},  {M_PI / 2.0, 0.1, -0.7}, {-M_PI / 2.0, 0.1, 0.0},
      {-M_PI / 2.0, 0.1, 0.3}, {-M_PI / 2.0, 0.1, 0.7}, {0.0, 0.1, 0.0},
      {M_PI / 2.0, -0.1, 0.0}, {M_PI / 3.0, 0.1, 0.2},  {M_PI, 0.1, 0.0},
  };

  for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
    tb_model_t *m =
        tb_model_haldane(t1, cases[i].t2, cases[i].phi, cases[i].mass);
    char label[96];
    int expect = tb_haldane_chern(t1, cases[i].t2, cases[i].phi, cases[i].mass);
    snprintf(label, sizeof label, "Chern: phi=%.3f t2=%.2f M=%.2f -> %d",
             cases[i].phi, cases[i].t2, cases[i].mass, expect);
    if (!m) {
      printf("  %s: FAIL (model)\n", label);
      failures++;

      continue;
    }

    if (expect == TB_CHERN_UNDEFINED) {
      double c = 0.0;
      check_true(tb_chern_number(m, 0, 1, 48, &c) == TB_ERR_DEGENERATE, label);
    } else {
      // lower band's Chern number
      check_close(chern(m, 0, 1, 48), expect, 1e-6, label);
      // filled lower band and empty upper band cancel
      check_close(chern(m, 0, 1, 48) + chern(m, 1, 1, 48), 0.0, 1e-6,
                  "  C_lower + C_upper = 0");
    }

    tb_model_free(m);
  }

  // Berry curvature sits at Dirac points of massive graphene
  tb_model_t *mg = tb_model_haldane(t1, 0.0, 0.0, 0.4);
  check_true(mg != NULL, "massive graphene built");
  if (mg) {
    enum { N = 36 };
    static double flux[N * N];
    double gap = 0.0;

    check_true(tb_berry_flux_grid(mg, 0, 1, N, flux, &gap) == TB_OK,
               "flux grid OK");
    check_close(gap, 0.8, 0.05, "minimum gap ~ 2M");

    double total = 0.0;
    double kflux = 0.0;
    double kpflux = 0.0;
    for (int a = 0; a < N; a++) {
      for (int b = 0; b < N; b++) {
        total += flux[a * N + b];
        // near K = (N/3, 2N/3), K' = (2N/3, N/3)
        if (abs(a - N / 3) <= 5 && abs(b - 2 * N / 3) <= 5) {
          kflux += flux[a * N + b];
        }
        if (abs(a - 2 * N / 3) <= 5 && abs(b - N / 3) <= 5) {
          kpflux += flux[a * N + b];
        }
      }
    }

    check_close(total, 0.0, 1e-9, "massive graphene: C = 0");

    printf("  valley flux: K=%.4f K'=%.4f\n", kflux, kpflux);
    check_true(fabs(kflux) > 1.5 && fabs(kpflux) > 1.5 && kflux * kpflux < 0.0,
               "opposite-sign flux at two valleys");

    tb_model_free(mg);
  }
}

static void test_ssh(void) {
  printf("  === SSH chain ===\n");

  const double vs[3] = {1.0, 0.4, 0.8};
  const double ws[3] = {0.4, 1.0, -1.1};

  for (int i = 0; i < 3; i++) {
    tb_model_t *m = tb_model_ssh(vs[i], ws[i]);

    char label[64];
    snprintf(label, sizeof label, "Zak phase v=%.1f w=%.1f", vs[i], ws[i]);
    if (!m) {
      failures++;
      continue;
    }

    double z = 0.0;
    check_true(tb_wilson_phase(m, 0, 1, 0.0, 200, &z) == TB_OK, label);
    double expect = tb_ssh_zak_phase(vs[i], ws[i]);

    // phase is defined mod 2\pi; compare on circle
    check_close(cos(z), cos(expect), 1e-9, "  \\cos(Z_ak) matches reference");
    check_close(fabs(sin(z)), 0.0, 1e-9, "  Zak phase quantized to 0 or \\pi");
    tb_model_free(m);
  }

  check_true(isnan(tb_ssh_zak_phase(0.7, 0.7)), "reference NaN when gapless");

  // Bulk-boundary: edge states of open chain appear exactly when Zak
  // phase is pi
  const double pairs[2][2] = {{0.5, 1.0}, {1.0, 0.5}};
  for (int i = 0; i < 2; i++) {
    cmatrix_t *h =
        lattice_build_ssh(40, pairs[i][0], pairs[i][1], LATTICE_OPEN);
    check_true(h != NULL, "open SSH chain built");
    if (!h) {
      continue;
    }

    eigen_t *eig = cmatrix_eigh(h);
    int near_zero = 0;
    if (eig) {
      for (int n = 0; n < eig->n; n++) {
        if (fabs(eig->eigenvalues[n]) < 1e-3) {
          near_zero++;
        }
      }
    }

    // NOTE: lattice.h's chain has t1 intracell, t2 intercell: edge states iff
    // t2 > t1
    double z = tb_ssh_zak_phase(pairs[i][0], pairs[i][1]);
    check_true((near_zero == 2) == (fabs(z) > 1.0),
               "edge states <=> Zak phase pi");

    eigen_free(eig);
    cmatrix_free(h);
  }
}

static void test_kagome(void) {
  printf("  === Kagome lattice ===\n");

  const double ts[2] = {1.0, -0.6};
  for (int s = 0; s < 2; s++) {
    double t = ts[s];
    tb_model_t *m = tb_model_kagome(t);

    check_true(m != NULL && tb_model_orbitals(m) == 3, "kagome built");
    if (!m) {
      failures++;

      continue;
    }

    double flat = tb_kagome_flat_band_energy(t);
    double worst_flat = 0.0;
    double worst_sum = 0.0;
    for (int i = 0; i < 11; i++) {
      for (int j = 0; j < 11; j++) {
        double k1 = 2.0 * M_PI * i / 11.0 + 0.2;
        double k2 = 2.0 * M_PI * j / 11.0 + 0.5;
        double e[3];

        tb_bands(m, k1, k2, e);

        // flat band = lowest for t > 0, highest for t < 0
        double fb = t > 0.0 ? e[0] : e[2];
        worst_flat = fmax(worst_flat, fabs(fb - flat));

        double others = t > 0.0 ? e[1] + e[2] : e[0] + e[1];
        worst_sum = fmax(worst_sum, fabs(others - 2.0 * t));
      }
    }

    check_close(worst_flat, 0.0, 1e-10, "flat band at -2t");
    check_close(worst_sum, 0.0, 1e-10, "dispersive bands sum to 2t");
    if (t > 0.0) {
      double e[3];

      tb_bands(m, 2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0, e);

      check_close(e[1], t, 1e-10, "Dirac point energy t (lower)");
      check_close(e[2], t, 1e-10, "Dirac point energy t (upper)");
    }

    tb_model_free(m);
  }
}

static void test_hofstadter(void) {
  printf("  === Hofstadter ===\n");

  // Bloch matrix equals lattice.h's
  const int pq[3][2] = {{1, 3}, {2, 5}, {1, 4}};
  double worst = 0.0;
  for (int c = 0; c < 3; c++) {
    int p = pq[c][0];
    int q = pq[c][1];
    tb_model_t *m = tb_model_hofstadter(p, q, 1.0);
    if (!m) {
      failures++;

      continue;
    }

    for (int i = 0; i < 6; i++) {
      double kx = 0.37 + 0.9 * i;
      double ky = 1.1 + 0.55 * i;

      cmatrix_t *a = tb_hamiltonian(m, kx, ky);
      cmatrix_t *b = lattice_hofstadter_bloch(kx, ky, p, q, 1.0);
      if (!a || !b) {
        failures++;
      } else {
        for (int e = 0; e < q * q; e++) {
          worst = fmax(worst, hypot(a->data[e].re - b->data[e].re,
                                    a->data[e].im - b->data[e].im));
        }
      }

      cmatrix_free(a);
      cmatrix_free(b);
    }

    tb_model_free(m);
  }

  check_close(worst, 0.0, 1e-12, "H(k) identical to lattice_hofstadter_bloch");

  // Chern numbers against lattice.h and TKNN values
  const int cases[3][2] = {{1, 3}, {2, 5}, {1, 5}};
  for (int c = 0; c < 3; c++) {
    int p = cases[c][0];
    int q = cases[c][1];

    tb_model_t *m = tb_model_hofstadter(p, q, 1.0);
    double ref[8] = {0};

    int ok = lattice_hofstadter_chern_numbers(p, q, 1.0, 40, ref);
    check_true(m != NULL && ok != 0, "reference Chern numbers computed");
    if (!m) {
      failures++;

      continue;
    }

    double sum = 0.0;
    for (int b = 0; b < q; b++) {
      char label[64];
      snprintf(label, sizeof label, "flux %d/%d band %d Chern", p, q, b);
      double got = chern(m, b, 1, 40);

      check_close(got, ref[b], 1e-6, label);

      sum += got;
    }

    check_close(sum, 0.0, 1e-6, "  sum of band Chern numbers = 0");
    // a group of all bands is whole spectrum: trivial
    check_close(chern(m, 0, q, 24), 0.0, 1e-9, "  all-band group C = 0");

    tb_model_free(m);
  }

  // Known values for 1/3: (+1, -2, +1)
  tb_model_t *m3 = tb_model_hofstadter(1, 3, 1.0);
  if (m3) {
    check_close(chern(m3, 0, 1, 40), 1.0, 1e-6, "1/3 lowest band C = +1");
    check_close(chern(m3, 1, 1, 40), -2.0, 1e-6, "1/3 middle band C = -2");
    check_close(chern(m3, 0, 2, 40), -1.0, 1e-6, "1/3 bands 0+1 C = -1");

    tb_model_free(m3);
  } else {
    failures++;
  }

  // Flux 1/4: middle bands touch. Band 1 alone is ill-defined, group
  // of lower two bands is well defined
  tb_model_t *m4 = tb_model_hofstadter(1, 4, 1.0);
  if (m4) {
    double c = 0.0;

    check_true(tb_chern_number(m4, 1, 1, 24, &c) == TB_ERR_DEGENERATE,
               "1/4 band 1 alone: TB_ERR_DEGENERATE");
    check_close(chern(m4, 1, 2, 32), -2.0, 1e-6,
                "1/4 touching pair (1+2) C = -2");
    check_close(chern(m4, 0, 1, 32), 1.0, 1e-6, "1/4 lowest band C = +1");
    check_close(chern(m4, 3, 1, 32), 1.0, 1e-6, "1/4 top band C = +1");

    tb_model_free(m4);
  } else {
    failures++;
  }
}

static void test_wilson_winding(void) {
  printf("  ===  Wilson-loop winding ===\n");

  tb_model_t *h = tb_model_haldane(1.0, 0.1, M_PI / 2.0, 0.2);
  tb_model_t *hm = tb_model_haldane(1.0, 0.1, -M_PI / 2.0, 0.2);
  tb_model_t *triv = tb_model_haldane(1.0, 0.1, M_PI / 2.0, 0.8);
  tb_model_t *hof = tb_model_hofstadter(2, 5, 1.0);
  tb_model_t *k3 = tb_model_hofstadter(1, 3, 1.0);
  if (!h || !hm || !triv || !hof || !k3) {
    failures++;
  } else {
    double w = NAN;

    check_true(tb_wilson_winding(h, 0, 1, 48, 48, &w) == TB_OK, "Haldane OK");
    check_close(w, chern(h, 0, 1, 48), 1e-6,
                "winding = Chern (\\phi = \\pi/2)");

    check_true(tb_wilson_winding(hm, 0, 1, 48, 48, &w) == TB_OK,
               "Haldane (-\\phi) OK");
    check_close(w, chern(hm, 0, 1, 48), 1e-6,
                "winding = Chern (\\phi = -\\pi/2)");

    check_true(tb_wilson_winding(triv, 0, 1, 48, 48, &w) == TB_OK,
               "trivial OK");
    check_close(w, 0.0, 1e-6, "trivial phase: winding 0");
    for (int b = 0; b < 5; b++) {
      char label[48];
      snprintf(label, sizeof label, "flux 2/5 band %d winding", b);

      check_true(tb_wilson_winding(hof, b, 1, 40, 80, &w) == TB_OK, label);
      check_close(w, chern(hof, b, 1, 40), 1e-6, "  = Chern");
    }

    check_true(tb_wilson_winding(k3, 0, 2, 40, 60, &w) == TB_OK, "1/3 group");
    check_close(w, -1.0, 1e-6, "1/3 bands 0+1 winding = -1");
  }

  tb_model_free(h);
  tb_model_free(hm);
  tb_model_free(triv);
  tb_model_free(hof);
  tb_model_free(k3);
}

static void test_invalid(void) {
  printf("  === Invalid input ===\n");

  check_true(tb_model_alloc(0, 2) == NULL, "dim 0 rejected");
  check_true(tb_model_alloc(3, 2) == NULL, "dim 3 rejected");
  check_true(tb_model_alloc(2, 0) == NULL, "0 orbitals rejected");
  check_true(tb_model_alloc(2, TB_MAX_ORBITALS + 1) == NULL,
             "too many orbitals rejected");
  check_true(tb_model_dim(NULL) == 0 && tb_model_orbitals(NULL) == 0,
             "NULL model dim/orbitals = 0");

  tb_model_free(NULL);

  tb_model_t *m = tb_model_alloc(1, 2);
  check_true(m != NULL, "1D 2-orbital model");
  if (m) {
    check_true(tb_set_onsite(m, 2, 1.0) == TB_ERR_INVALID,
               "onsite orbital out of range");
    check_true(tb_set_onsite(m, 0, NAN) == TB_ERR_INVALID, "non-finite onsite");
    check_true(tb_add_hopping(m, 0, 0, 0, 0, c_new(1, 0)) == TB_ERR_INVALID,
               "on-cell self hopping rejected");
    check_true(tb_add_hopping(m, 0, 1, 0, 1, c_new(1, 0)) == TB_ERR_INVALID,
               "n2 != 0 rejected in 1D");
    check_true(tb_add_hopping(m, 0, 2, 0, 0, c_new(1, 0)) == TB_ERR_INVALID,
               "hop orbital out of range");
    check_true(tb_add_hopping(m, 0, 1, TB_MAX_HOP_RANGE + 1, 0, c_new(1, 0)) ==
                   TB_ERR_INVALID,
               "hop range exceeded");
    check_true(tb_add_hopping(m, 0, 1, 0, 0, c_new(INFINITY, 0)) ==
                   TB_ERR_INVALID,
               "non-finite amplitude");
    check_true(tb_add_hopping(m, 0, 1, 0, 0, c_new(1, 0)) == TB_OK,
               "valid hopping accepted");

    double e[2];
    check_true(tb_bands(m, 0.0, 0.0, NULL) == TB_ERR_INVALID, "NULL energies");
    check_true(tb_bands(NULL, 0.0, 0.0, e) == TB_ERR_INVALID, "NULL model");
    check_true(tb_bands(m, NAN, 0.0, e) == TB_ERR_INVALID, "NaN k");

    double c = 0.0;
    check_true(tb_chern_number(m, 0, 1, 10, &c) == TB_ERR_INVALID,
               "Chern needs a 2D model");
    check_true(tb_hamiltonian(NULL, 0, 0) == NULL, "hamiltonian NULL model");

    tb_model_free(m);
  }

  tb_model_t *g = tb_model_graphene(1.0);
  if (g) {
    double c = 0.0;
    double ph = 0.0;

    check_true(tb_chern_number(g, 0, 0, 10, &c) == TB_ERR_INVALID,
               "zero bands rejected");
    check_true(tb_chern_number(g, 1, 2, 10, &c) == TB_ERR_INVALID,
               "band group past top rejected");
    check_true(tb_chern_number(g, 0, 1, 1, &c) == TB_ERR_INVALID,
               "n_k too small");
    check_true(tb_chern_number(g, 0, 1, TB_MAX_GRID + 1, &c) == TB_ERR_INVALID,
               "n_k too large");
    check_true(tb_chern_number(g, 0, 1, 12, &c) == TB_ERR_DEGENERATE,
               "gapless graphene: DEGENERATE");
    check_true(tb_chern_number(g, 0, 1, 12, NULL) == TB_ERR_INVALID,
               "NULL output");
    check_true(tb_wilson_phase(g, 0, 1, 0.3, 1, &ph) == TB_ERR_INVALID,
               "wilson n_k too small");
    check_true(tb_wilson_winding(g, 0, 1, 20, 3, &ph) == TB_ERR_INVALID,
               "winding n_k2 < 4");

    double e[16];
    check_true(tb_dos(g, 8, 0.0, -1.0, 4, 0.1, e) == TB_ERR_INVALID,
               "DOS e_min >= e_max");
    check_true(tb_dos(g, 8, -1.0, 1.0, 4, 0.0, e) == TB_ERR_INVALID,
               "DOS sigma <= 0");
    check_true(tb_bands_along_line(g, 0, 0, 1, 1, 1, e) == TB_ERR_INVALID,
               "band line needs >= 2 points");

    tb_model_free(g);
  }

  check_true(tb_model_chain(NAN) == NULL, "chain(NaN) rejected");
  check_true(tb_model_hofstadter(3, 3, 1.0) == NULL, "flux p >= q rejected");
  check_true(tb_model_hofstadter(0, 3, 1.0) == NULL, "flux p < 1 rejected");
  check_true(tb_model_haldane(1.0, NAN, 0.0, 0.0) == NULL,
             "haldane(NaN) rejected");
  check_true(tb_haldane_chern(1.0, 0.1, M_PI / 2.0, 3.0 * sqrt(3.0) * 0.1) ==
                 TB_CHERN_UNDEFINED,
             "phase boundary: UNDEFINED");
  check_true(tb_haldane_chern(1.0, 0.9, M_PI / 2.0, 0.0) == TB_CHERN_UNDEFINED,
             "outside |t2| <= |t1|/3: UNDEFINED");
}

int main(void) {
  printf(" > Tight-binding tests:\n");

  test_dispersions();
  test_haldane();
  test_ssh();
  test_kagome();
  test_hofstadter();
  test_wilson_winding();
  test_invalid();

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll tight-binding tests passed\n");
  return 0;
}
