#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/tred2.c"
#include "/opt/NR/numerical_recipes.c/tqli.c"
#include "/opt/NR/numerical_recipes.c/pythag.c"


void print_matrix(float** M, int dim_n, int dim_m) {
  printf("\n");
  for (int i = 1; i <= dim_n; ++i) {
		for (int j = 1; j <= dim_m; ++j)
			printf("%5.3f ", M[i][j]);
    printf("\n");
	} 
  printf("\n");
}

float dot(float* a, float* b, int n) {
  float res = 0;
  for (int i = 1; i <= n; i++)
    res += a[i] * b[i];
  return res;
}

float norm_vec(float* vec, int n) {
  float norm = sqrt(fabs(dot(vec, vec, n)));
  for (int i = 1; i <= n; i++)
    vec[i] /= norm;
}

void aprox_x1(float* x1, float** Wk, float* x, int n) {
  for (int i = 1; i <= n; i++)
    x1[i] = dot(Wk[i], x, n);
} 

void hot_red(float** W, float eigen_val, float* x, int n) {
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      W[i][j] -= eigen_val * x[i] * x[j];
    }
  } 
}

float calc_eigen_val(float* x1, float* x0, int n) {
  return dot(x1, x0, n) / dot(x0, x0, n);
}

void copy_vec_to(float* src_vec, float* dest_vec, int n) {
  for (int i = 1; i <= n; i++) 
    dest_vec[i] = src_vec[i];
}

int main(void) {
  int n = 7;

  float** A = matrix(1, n, 1, n);
  float** W = matrix(1, n, 1, n);

  float* d = vector(1, n);
  float* e = vector(1, n);
  
  float* x0 = vector(1, n);
  float* x1 = vector(1, n);

  float eigen_val;

  // initializing arrays

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      A[i][j] = W[i][j] = sqrt(i + j);
    }
  }

  printf("=== INPUT MATRIX A: === \n");
  print_matrix(A, n, n);

  tred2(A, n, d, e);

  tqli(d, e, n, A);

  printf("=== A containing eigen vectors === \n");
  print_matrix(A, n, n);

  printf("=== vector d with eigen valeus: === \n\n");
  for (int i = 1; i <= n; i++)
	  printf("%.10f \n", d[i]);
  printf("\n\n");


  printf("=== ITERATION METHOD === \n\n");
  for (int k = 1; k <= n; k++) {
    for (int j = 0; j <= n; j++)
	    x0[j] = 1;
	
	for (int i = 1; i<= 8; i++) {
      aprox_x1(x1, W, x0, n);
      eigen_val = calc_eigen_val(x1, x0, n);
      norm_vec(x1, n);
      copy_vec_to(x1, x0, n);
    }
      hot_red(W, eigen_val, x0, n);
      printf("%.10f\n", eigen_val);
  }


  free_matrix(A, 1, n, 1, n);
  free_matrix(W, 1, n, 1, n);

  free_vector(d, 1, n);
  free_vector(e, 1, n);
  free_vector(x0, 1, n);
  free_vector(x1, 1, n);
  return 0;
}

