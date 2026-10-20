/* Native reference for wasm/qmc_hubbard.c; diffed against hubbard_check.mjs. */
#include <math.h>
#include <stdio.h>

int qmc_hb_curves(double umax);
double qmc_hb_levels(double u);
int qmc_hb_quench(double u, double tmax);
int qmc_hb_probs(double u, int mode, double time);
int qmc_hb_thermal(double u, double tmax);
double *qmc_hb_spec_buffer(void);
double *qmc_hb_gs_d_buffer(void);
double *qmc_hb_gs_ss_buffer(void);
double *qmc_hb_gs_s_buffer(void);
double *qmc_hb_q_d_buffer(void);
double *qmc_hb_q_m_buffer(void);
double *qmc_hb_q_s_buffer(void);
double *qmc_hb_q_ss_buffer(void);
double *qmc_hb_th_d_buffer(void);
double *qmc_hb_th_ss_buffer(void);
double *qmc_hb_th_s_buffer(void);
double *qmc_hb_th_e_buffer(void);
double *qmc_hb_probs_buffer(void);
double *qmc_hb_levels_buffer(void);

/* Printed to 5 digits so libm last-bit differences cannot show up. */
int main(void) {
  qmc_hb_curves(12.0);
  const double *sp = qmc_hb_spec_buffer(), *d = qmc_hb_gs_d_buffer(), *ss = qmc_hb_gs_ss_buffer(), *se = qmc_hb_gs_s_buffer();
  printf("spec l0 i0=%.5f i100=%.5f i200=%.5f l5 i100=%.5f\n", sp[0], sp[100], sp[200], sp[5 * 201 + 100]);
  printf("gs d i0=%.5f i50=%.5f i200=%.5f\n", d[0], d[50], d[200]);
  printf("gs ss i0=%.5f i50=%.5f i200=%.5f\n", ss[0], ss[50], ss[200]);
  printf("gs s i0=%.5f i50=%.5f i200=%.5f\n", se[0], se[50], se[200]);
  double e0 = qmc_hb_levels(4.0);
  const double *lv = qmc_hb_levels_buffer();
  printf("levels e0=%.5f %.5f %.5f %.5f %.5f %.5f\n", e0, lv[1], lv[2], lv[3], lv[4], lv[5]);
  qmc_hb_quench(8.0, 20.0);
  const double *qd = qmc_hb_q_d_buffer(), *qm = qmc_hb_q_m_buffer(), *qs = qmc_hb_q_s_buffer(), *qss = qmc_hb_q_ss_buffer();
  printf("quench d q0=%.5f q100=%.5f q400=%.5f\n", qd[0], qd[100], qd[400]);
  printf("quench m q0=%.5f q100=%.5f q400=%.5f\n", qm[0], qm[100], qm[400]);
  printf("quench s q100=%.5f q300=%.5f ss q100=%.5f\n", qs[100], qs[300], qss[100]);
  qmc_hb_probs(4.0, 0, 0.0);
  const double *p = qmc_hb_probs_buffer();
  printf("probs gs %.5f %.5f %.5f %.5f | %.5f %.5f %.5f %.5f\n", p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
  qmc_hb_probs(4.0, 1, 1.3);
  printf("probs neel %.5f %.5f %.5f %.5f | %.5f %.5f %.5f %.5f\n", p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
  qmc_hb_thermal(4.0, 5.0);
  const double *td = qmc_hb_th_d_buffer(), *tss = qmc_hb_th_ss_buffer(), *ts = qmc_hb_th_s_buffer(), *te = qmc_hb_th_e_buffer();
  printf("thermal d %.5f %.5f %.5f ss %.5f %.5f %.5f\n", td[0], td[100], td[200], tss[0], tss[100], tss[200]);
  printf("thermal s %.5f %.5f %.5f e %.5f %.5f %.5f\n", ts[0], ts[100], ts[200], te[0], te[100], te[200]);
  if (qmc_hb_curves(0.0) != -1 || qmc_hb_quench(-1.0, 5.0) != -1 || qmc_hb_quench(2.0, 0.0) != -1 || qmc_hb_probs(2.0, 2, 0.0) != -1 ||
      qmc_hb_thermal(2.0, 0.0) != -1 || !isnan(qmc_hb_levels(99.0))) {
    printf("bad args accepted\n");

    return 1;
  }
  printf("bad args rejected\n");

  return 0;
}
