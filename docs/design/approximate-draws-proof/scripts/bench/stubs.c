#include <stdio.h>
#include <stdlib.h>
#include <Rmath.h>
void ext_throwError(const char* format, ...) { fprintf(stderr, "error\n"); exit(1); }
void ext_issueWarning(const char* format, ...) { fprintf(stderr, "warning\n"); }
int rc_getRuntimeVersion(int* major, int* minor, int* revision) { *major = 4; *minor = 5; *revision = 0; return 0; }
double dbarts_nmath_pnorm5(double x, double m, double s, int l, int g) { return pnorm5(x, m, s, l, g); }
double dbarts_nmath_qnorm5(double p, double m, double s, int l, int g) { return qnorm5(p, m, s, l, g); }
