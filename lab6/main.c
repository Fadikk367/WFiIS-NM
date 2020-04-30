#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 5
#define IT_MAX 30

double licz_r(double*, double*, int, double);

int main(void) {

  double* a = (double*) calloc(N + 1, sizeof(double));
  double* b = (double*) calloc(N + 1, sizeof(double));
  double* c = (double*) calloc(N, sizeof(double));
  
  a[0] = 240.0;
  a[1] = -196.0;
  a[2] = -92.0;
  a[3] = 33.0;
  a[4] = 14.0;
  a[5] = 1.0;

  int n;
  double x0, x1, Rj, Rj_prim;
  for (int L = 1; L <= N; ++L) {
    n = N - L + 1;
    x0 = 0.0;
    for (int it = 1; it <= IT_MAX; ++it) {

      Rj = licz_r(a, b, n, x0);
      Rj_prim = licz_r(b, c, n-1, x0);
      x1 = x0 - Rj/Rj_prim;

      printf("%d, %d, %.8lf, %.8lf, %.8lf \n", L, it, x1, Rj, Rj_prim);

      if (fabs(x1 - x0) < 1.0e-7) break;
      x0 = x1;
    }
    for (int i = 0; i <= (n - 1); ++i) a[i] = b[i];
    printf("\n=========================\n\n");
  }
  return 0;
}


double licz_r(double* a, double* b, int n, double x_0) {
  b[n] = 0.0;
  for (int k = n - 1; k >= 0; --k) {
    b[k] = a[k + 1] + x_0*b[k + 1];
  }

  return a[0] + x_0*b[0];
}

