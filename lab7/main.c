#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define x_max 5.0
#define x_min -5.0

double fun(double x) {
  return 1.0 / (1.0 + x*x);
}

double interpolateNewton(double x, int n, double* xm, double** fm) {
  double result = 0;
  for (int j = 0; j <= n; ++j) {
    double product = fm[j][j];
    for (int i = 0; i <= j - 1; i++) 
      product *= x - xm[i];

    result += product;
  }
  return result;
}

double calcCzebyszewZero(int i, int n) {
  return 0.5*((x_min - x_max)*cos(M_PI*(2*i + 1)/(2*n + 2)) + (x_min + x_max));
}

int main(void) {
  
  int steps[4] = {5, 10, 15, 20};
  for (int k = 0; k < 4; ++k) {
    int n = steps[k];
    double* xm = (double*) calloc(n + 1, sizeof(double));
    double* ym = (double*) calloc(n + 1, sizeof(double));

    double** fm = (double**) calloc(n + 1, sizeof(double*));
    for (int i = 0; i <= n; ++i)
      fm[i] = (double*) calloc(n + 1, sizeof(double));


    double h = (x_max - x_min) / n;
    for (int i = 0; i <= n; ++i) {
	  /// xm[i] = x_min + i*h;
      xm[i] = calcCzebyszewZero(i, n);
      ym[i] = fun(xm[i]);
    }

    for (int i = 0; i <= n; ++i) {
      fm[i][0] = ym[i];
    }

    for (int i = 1; i <= n; ++i)
      for (int j = 1; j <= i; ++j)
        fm[i][j] = (fm[i][j-1] - fm[i-1][j-1]) / (xm[i] - xm[i-j]);

    for (float x = x_min; x <= x_max; x += 0.01) 
      printf("%f  %f\n", x, interpolateNewton(x, n, xm, fm));
    
    printf("\n\n");

    free(ym);
    free(xm);
    for (int i = 0; i <= n; ++i)
      free(fm[i]);
    free(fm);
  }

  return 0;
}
