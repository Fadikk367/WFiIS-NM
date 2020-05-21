#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/tred2.c"
#include "/opt/NR/numerical_recipes.c/tqli.c"
#include "/opt/NR/numerical_recipes.c/pythag.c"

#define nx 20
#define ny 20
#define n nx*ny
#define m 10
#define t -0.021

int main(void) {

  float** H = matrix(1, n, 1, n);
  float** Y = matrix(1, n, 1, n);
  float** X = matrix(1, n, 1, n);

  float* d = vector(1, n);
  float* e = vector(1, n);
  int *indx = ivector(1, n);

  for (int i = 1; i <= nx; i++) {
    for (int j = 1; j <= ny; j++) {
      int l = j + (i - 1)*ny;
      for (int k = 1; k <= n; k++) H[l][k] = 0.;
      if (i >1 ) H[l][l - ny] = t; 
      if (i < nx) H[l][l + ny] = t; 
      H[l][l] = -4*t;
      if (j > 1) H[l][l-1] = t;
      if (j < ny) H[l][l+1] = t;
    }
  }

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      if (i == j) Y[i][j] = 1;
      else Y[i][j] = 0;
    }
  }

  tred2(H, n, d, e);

  float** P = H;

  tqli(d, e, n, Y);

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      float row_sum = 0;
      for (int k = 1; k <= n; k++)
        row_sum += P[i][k] * Y[k][j];
      X[i][j] = row_sum;
    }
  }


  for (int l = 1; l <= n; l++) indx[l] = l; 
  for (int l = 1; l <= n - 1; l++) {
    for (int k = n; k >= l + 1; k--) {
      float e1 = d[k - 1];
      float e2 = d[k];
      int l1 = indx[k - 1];
      int l2 = indx[k];
      if(e2 < e1) { 
        d[k] = e1;
        d[k - 1] = e2;
        indx[k] = l1;
        indx[k - 1] = l2;
      }
    }
  }


  FILE *fp;
  fp = fopen("dane.dat","w");
  for (int i = 1; i <= nx; i++) {
    for (int j = 1; j <= ny; j++) {
      int l = j + (i - 1)*ny;
      fprintf(fp, "%6d %6d ", i, j);
      for (int k = 1; k <= m; k++) fprintf(fp, " %12.6f ", X[l][indx[k]]);
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
  fclose(fp);


  free_matrix(H, 1, n, 1, n);
  free_matrix(Y, 1, n, 1, n);
  free_matrix(X, 1, n, 1, n);

  free_vector(d, 1, n);
  free_vector(e, 1, n);
  free_ivector(indx, 1, n);
  return 0;
}

