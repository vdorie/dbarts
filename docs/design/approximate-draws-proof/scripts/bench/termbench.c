/* Time per term evaluation of the certified evaluator's inner loops in plain double, without the
 * rigorous radius bookkeeping of Appendix A (which adds arithmetic but no further library calls):
 *   fnode: one PG node by eq. (7): complex sqrt, expm1, exp, log1p, exp, a few products;
 *   qterm: one term a_j log1p(c_j Delta) of the finite family's node (log1p of a complex number);
 *   qcell: one term of a cell majorant of q: log of |z^2 + tau'^2|^2 at the clamped point.
 * Build: cc -O2 -o termbench termbench.c -lm.  Run: ./termbench */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

static double complex cexpm1_(double complex z) {
  double a = creal(z), b = cimag(z), s = sin(0.5 * b);
  return (expm1(a) * cos(b) - 2.0 * s * s) + I * (exp(a) * sin(b));
}
static double complex clog1p_(double complex z) {
  double a = creal(z), b = cimag(z);
  return 0.5 * log1p(2.0 * a + a * a + b * b) + I * atan2(b, 1.0 + a);
}
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

int main(void) {
  const int n = 20000000;
  volatile double sink = 0.0;
  double rho = 1.3, x = 2.0, bb = 64.0, W0r = 3.0;
  double complex w0 = csqrt(W0r), e2w0 = 1.0 + cexp(-2.0 * w0);
  double t = now(), acc = 0.0;
  for (int k = 1; k <= n; ++k) {
    double y = 1e-7 * k;
    double complex D = -y * y + I * (2.0 * rho * y), w = csqrt(W0r + D), d = D / (w + w0);
    double complex X = cexpm1_(d) * (1.0 - cexp(-(w + w0))) / e2w0;
    double complex E = 2.0 * x * D - bb * clog1p_(X);
    acc += creal(cexp(E) * (rho + I * y));
  }
  double tf = (now() - t) / n; sink += acc;
  t = now(); acc = 0.0;
  for (int k = 1; k <= n; ++k) {
    double c = 1.0 / (10.0 + 1e-6 * k);
    double complex D = -0.01 + I * (1e-8 * k);
    acc += 0.7 * creal(clog1p_(c * D));
  }
  double tq = (now() - t) / n; sink += acc;
  t = now(); acc = 0.0;
  for (int k = 1; k <= n; ++k) {
    double tau = 10.0 + 1e-6 * k, Y = 0.3, u = 1.2, A = u * u + tau - Y;
    acc += 0.7 * log(A * A + 4.0 * u * u * Y);
  }
  double tc = (now() - t) / n; sink += acc;
  printf("ns per term: fnode %.1f  qterm %.1f  qcell %.1f  (sink %g)\n", 1e9 * tf, 1e9 * tq, 1e9 * tc, (double) sink);
  return 0;
}
