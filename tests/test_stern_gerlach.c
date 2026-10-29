/*
 * Test: Stern-Gerlach analyzers
 *
 * 1. One analyzer: P(+) = \cos^2(\theta / 2) for a z-polarised beam, outputs
 *    sum to inpu
 * 2. Sequential analyzers: z+ -> x+ -> z+ passes 1/4 while z+ -> z- passes
 *    nothing; three-analyzer chain agrees with a spinor calculation
 * 3. Rotating analyzers \by pi / N each: transmission \cos^(2N)(\pi / 2N) tends
 *    to 1
 * 4. Passing both outputs: perfect recombination leaves beam alone (Feynman
 *    "Humpty Dumpty" box), no recombination decoheres it
 * 5. Screen pattern: normalised, two spots, packets overlap \exp(-d^2 / 8
 *    \sigma^2) (numerical quadrature), classical smear normalised
 * 6. Spreading width and invalid input
 */

#include "../physics/stern_gerlach.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);

  printf("  %s: got=%.8g expected=%.8g err=%.2e\n", label, got, expected, err);
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

/* Independent spinor calculation: |<n, +|\psi>|^2 with \psi = (a, b)
 * real-pair amplitudes, projecting onto +n eigenspinor
 * (\cos(\theta/2), \exp(i \phi) \sin(\theta/2))
 * Returns probability and new spinor */
typedef struct {
  double ar, ai, br, bi;
} spinor_t;

static double project(spinor_t *s, double theta, double phi, int down) {
  double c = cos(theta / 2);
  double sn = sin(theta / 2);
  // eigenspinors: up = (c, e^{i phi} sn), down = (-\exp(-i \phi) sn, c)
  double ur = down ? -sn * cos(phi) : c, ui = down ? sn * sin(phi) : 0.0;
  double vr = down ? c : sn * cos(phi), vi = down ? 0.0 : sn * sin(phi);
  // amplitude <u|\psi> = conj(u) a + conj(v) b
  double amr = ur * s->ar + ui * s->ai + vr * s->br + vi * s->bi;
  double ami = ur * s->ai - ui * s->ar + vr * s->bi - vi * s->br;
  double p = amr * amr + ami * ami;

  if (p > 1e-15) {
    double nr = amr / sqrt(p), ni = ami / sqrt(p);

    s->ar = nr * ur - ni * ui;
    s->ai = nr * ui + ni * ur;
    s->br = nr * vr - ni * vi;
    s->bi = nr * vi + ni * vr;
  }

  return p;
}

static void test_single_analyzer(void) {
  printf("  === Test: Single analyzer ===\n");

  sg_vec_t zup = {0, 0, 1};
  for (int k = 0; k <= 4; k++) {
    double th = k * M_PI / 4;
    sg_analyzer_t up = {th, 0.3, SG_PASS_UP, 1.0};
    sg_analyzer_t dn = {th, 0.3, SG_PASS_DOWN, 1.0};
    double iu[2];
    double id[2];
    sg_vec_t ru[2];
    sg_vec_t rd[2];

    sg_chain(zup, 1.0, &up, 1, iu, ru);
    sg_chain(zup, 1.0, &dn, 1, id, rd);
    check_close(iu[1], cos(th / 2) * cos(th / 2), 1e-14,
                "P(+) = \\cos^2(\\theta/2)");
    check_close(iu[1] + id[1], 1.0, 1e-14, "outputs sum to input");
  }
}

static void test_sequential_analyzers(void) {
  printf("  === Test: Sequential analyzers ===\n");

  sg_vec_t zup = {0, 0, 1};
  sg_analyzer_t a[4];
  double inten[5];
  sg_vec_t r[5];

  sg_analyzer_t zxz[3] = {{0, 0, SG_PASS_UP, 1},
                          {M_PI / 2, 0, SG_PASS_UP, 1},
                          {0, 0, SG_PASS_UP, 1}};

  sg_chain(zup, 1.0, zxz, 3, inten, r);
  check_close(inten[3], 0.25, 1e-14, "z+ -> x+ -> z+ passes 1/4");

  sg_analyzer_t zz[2] = {{0, 0, SG_PASS_UP, 1}, {0, 0, SG_PASS_DOWN, 1}};

  sg_chain(zup, 1.0, zz, 2, inten, r);
  check_close(inten[2], 0.0, 1e-14, "z+ -> z- passes nothing");
  check_close(fabs(r[1].z - 1.0), 0.0, 1e-14,
              "beam after first analyzer is z+");

  // Random chain against a spinor calculation
  double th[4] = {0.7, 2.1, 1.3, 2.8};
  double ph[4] = {0.4, 5.0, 1.1, 3.3};
  int mode[4] = {SG_PASS_UP, SG_PASS_DOWN, SG_PASS_DOWN, SG_PASS_UP};
  sg_vec_t r0 = sg_axis(1.1, 0.6);
  spinor_t s = {cos(0.55), 0.0, sin(0.55) * cos(0.6), sin(0.55) * sin(0.6)};
  double pref = 1.0;

  for (int k = 0; k < 4; k++) {
    a[k].theta = th[k];
    a[k].phi = ph[k];
    a[k].mode = mode[k];
    a[k].coherence = 1.0;

    pref *= project(&s, th[k], ph[k], mode[k] == SG_PASS_DOWN);
  }

  sg_chain(r0, 1.0, a, 4, inten, r);
  check_close(inten[4], pref, 1e-12, "4 analyzers versus spinor calculation");
}

static void test_zeno_staircase(void) {
  printf("  === Test: Zeno staircase ===\n");

  sg_vec_t zup = {0, 0, 1};
  for (int n = 2; n <= 8; n += 3) {
    sg_analyzer_t st[8];
    double inten[9];
    sg_vec_t r[9];

    for (int k = 0; k < n; k++) {
      st[k].theta = (k + 1) * M_PI / n;
      st[k].phi = 0.0;
      st[k].mode = SG_PASS_UP;
      st[k].coherence = 1.0;
    }

    // last analyzer points along -z, so beam is turned by pi through n small
    // steps
    sg_chain(zup, 1.0, st, n, inten, r);
    check_close(inten[n], pow(cos(M_PI / (2.0 * n)), 2.0 * n), 1e-13,
                "transmission \\cos^(2N)(\\pi/2N)");
  }
  check_true(pow(cos(M_PI / 400.0), 200.0) > 0.99,
             "N = 100 passes more than 99 %");
}

static void test_recombination(void) {
  printf("  === Test: Recombination ===\n");

  sg_vec_t r0 = sg_axis(1.1, 0.6);
  sg_analyzer_t both = {1.2, 0.5, SG_PASS_BOTH, 1.0};
  double inten[9];
  sg_vec_t r[9];

  sg_chain(r0, 1.0, &both, 1, inten, r);
  check_close(inten[1], 1.0, 1e-15, "intensity kept");

  double dx = r[1].x - r0.x;
  double dy = r[1].y - r0.y;
  double dz = r[1].z - r0.z;
  check_close(sqrt(dx * dx + dy * dy + dz * dz), 0.0, 1e-14,
              "perfect recombination restores the state");

  both.coherence = 0.0;
  sg_chain(r0, 1.0, &both, 1, inten, r);
  sg_vec_t n1 = sg_axis(1.2, 0.5);
  double along = r0.x * n1.x + r0.y * n1.y + r0.z * n1.z;

  check_close(r[1].x * n1.x + r[1].y * n1.y + r[1].z * n1.z, along, 1e-14,
              "decohered beam keeps r.n");
  check_close(sqrt(r[1].x * r[1].x + r[1].y * r[1].y + r[1].z * r[1].z),
              fabs(along), 1e-14, "and loses transverse part");

  // decohered beam into a perpendicular analyzer: 1/2, coherent recombination
  // into original axis: full memory
  sg_analyzer_t chain2[2] = {{M_PI / 2, 0, SG_PASS_BOTH, 0.0},
                             {0, 0, SG_PASS_UP, 1.0}};
  sg_vec_t xup = {1, 0, 0};

  sg_chain(xup, 1.0, chain2, 2, inten, r);
  check_close(inten[2], 0.5, 1e-14,
              "x+ beam through a decohering x box then z+: 1/2");

  chain2[0].theta = 0.0;
  chain2[0].coherence = 0.0;
  sg_chain(xup, 1.0, chain2, 2, inten, r);
  check_close(inten[2], 0.5, 1e-14,
              "x+ beam through a decohering z box then z+: 1/2");

  chain2[0].coherence = 1.0;
  chain2[1].theta = M_PI / 2;
  sg_chain(xup, 1.0, chain2, 2, inten, r);
  check_close(inten[2], 1.0, 1e-14,
              "x+ beam through a perfect z box then x+: 1");

  chain2[0].coherence = 0.0;
  sg_chain(xup, 1.0, chain2, 2, inten, r);
  check_close(inten[2], 0.5, 1e-14,
              "x+ beam through a decohering z box then x+: 1/2");
}

static void test_screen_pattern(void) {
  printf("  === Test: Screen pattern ===\n");

  const double d = 3.0;
  const double sg = 0.6;
  const double pup = 0.3;
  double sum = 0;
  double sumc = 0;
  double ov = 0;
  double dz = 0.002;

  for (double z = -12; z < 12; z += dz) {
    sum += sg_screen_density(z, pup, d, sg) * dz;
    sumc += sg_classical_density(z, d, sg) * dz;
    ov += sqrt(sg_gauss(z, 0.5 * d, sg) * sg_gauss(z, -0.5 * d, sg)) * dz;
  }

  check_close(sum, 1.0, 1e-6, "pattern normalised");
  check_close(sumc, 1.0, 1e-6, "classical smear normalised");
  check_close(ov, sg_coherence(d, sg), 1e-6,
              "Bhattacharyya overlap of packets = \\exp(-d^2/8 \\sigma^2)");

  check_true(sg_screen_density(0.5 * d, pup, d, sg) >
                 3 * sg_screen_density(0.0, pup, d, sg),
             "2 spots for d = 5 \\sigma");

  check_true(sg_classical_density(0.0, d, sg) >
                 0.9 * sg_classical_density(0.4 * d, d, sg),
             "classical smear has no gap");

  check_close(sg_screen_density(-1.1, pup, d, sg),
              sg_screen_density(1.1, 1.0 - pup, d, sg), 1e-15,
              "mirror symmetry");
}

static void test_spreading_and_invalid_input(void) {
  printf("  === Test: Spreading and invalid input ===\n");

  sg_vec_t zup = {0, 0, 1};
  sg_analyzer_t both = {1.2, 0.5, SG_PASS_BOTH, 1.0};
  double inten[9];
  sg_vec_t r[9];

  check_close(sg_width(0.5, 0.0), 0.5, 1e-15, "width at t = 0");
  check_close(sg_width(0.5, 0.5), 0.5 * sqrt(2.0), 1e-14,
              "width at t = 2 \\sigma0^2");

  sg_vec_t big = {2, 0, 0};
  sg_analyzer_t bad = {0, 0, SG_PASS_UP, 2.0};

  check_true(sg_chain(big, 1.0, &both, 1, inten, r) == -1, "|r| > 1 rejected");
  check_true(sg_chain(zup, 1.0, &bad, 1, inten, r) == -1,
             "coherence out of range rejected");
  check_true(sg_chain(zup, -1.0, &both, 1, inten, r) == -1,
             "negative intensity rejected");
  check_true(sg_chain(zup, 1.0, &both, -1, inten, r) == -1,
             "negative count rejected");
  check_true(sg_chain(zup, NAN, &both, 1, inten, r) == -1, "NaN rejected");
}

int main(void) {
  printf(" > Stern-Gerlach analyzers tests:\n");

  test_single_analyzer();
  test_sequential_analyzers();
  test_zeno_staircase();
  test_recombination();
  test_screen_pattern();
  test_spreading_and_invalid_input();

  if (failures == 0) {
    printf("\nAll test_stern_gerlach checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
