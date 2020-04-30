#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/gaussj.c"

#define x_min -5.f
#define x_max 5.f

#define alpha 0.f
#define betha 0.f

float fun1(float x) { return 1 / (1 + x*x); }
float fun2(float x) { return cos(2*x); }

void printMatrix(float** M, int n, int m) {
  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= m; ++j) 
      printf("%.4f ", M[i][j]);
    printf("\n");
  }
  printf("\n");
}

void printVector(float* v, int n) {
  for (int i = 1; i <= n; ++i) 
	printf("%f ", v[i]);
  printf("\n\n");
}

void fillZerosVector(float* v, int n) {
  for (int i = 1; i <= n; ++i)
    v[i] = 0.f;
}

void fillZerosMatrix(float** M, int n, int m) {
  for (int i = 1; i <= n; ++i)
    fillZerosVector(M[i], m);
}

float derivative2nd(float (*fun)(float), float x, float dx) {
  return (fun(x - dx) - 2*fun(x) + fun(x + dx)) / (dx*dx);
}

void wyzM(float* xw, float* yw, float* m, int n, float alph, float beth) {
  float** A = matrix(1, n, 1, n);
  float** d = matrix(1, n, 1, 1);
  fillZerosMatrix(A, n, n);
  fillZerosMatrix(d, n, 1);

  // Fill M matrix and d (m) vector
  A[1][1] = A[n][n] = 1.f;
  A[1][2] = A[n][n-1] = 0.f;

  d[1][1] = m[1] = alph;
  d[n][1] = m[n] = beth;

  float lambda, my, di, h = (x_max - x_min) / (n - 1);
  for (int i = 2; i < n; ++i) {
    lambda = h / (h + h);
    my = 1 - lambda;
    di = 6/(h + h) * ((yw[i+1] - yw[i])/h - (yw[i] - yw[i-1])/h);
    A[i][i] = 2;
    A[i][i-1] = my;
    A[i][i+1] = lambda;
    d[i][1] = di;
  }

  gaussj(A, n, d, 1);

  // copy result vector to m
  for (int i = 2; i < n; ++i)
    m[i] = d[i][1];

  free_matrix(A, 1, n, 1, n);
  free_matrix(d, 1, n, 1, 1);
}

float wyzSx(float* xw, float* yw, float* m, int n, float x) {
  int i= 1;
  for (int k = 2; k <= n; ++k)
    if (xw[k-1] <= x && x <= xw[k]) i = k;
  float h = (x_max - x_min) / (n - 1); 
  float A = (yw[i] - yw[i-1])/h - h/6 * (m[i] - m[i-1]);
  float B = yw[i-1] - m[i-1]*h*h/6;
  float s = m[i-1] * pow(xw[i] - x, 3)/(6*h) + m[i] * pow(x - xw[i-1], 3)/(6*h) + A*(x - xw[i-1]) + B;
  return s;
}

void printResults(float* xw, float* yw, float* m, float a, float b, int n, float step) {
  for (float x = a; x <= b; x += step)
    printf("%f  %f\n", x, wyzSx(xw, yw, m, n, x));
  printf("\n\n");

  //for (int i = 1; i <= n; ++i)
  //  printf("%f  %f  %f\n", xw[i], m[i], derivative2nd(fun1, xw[i], 0.01));
  //printf("\n\n");
}


void interpolate(float (*fun)(float), float a, float b, int n) {
  float* xw = vector(1, n);
  float* yw = vector(1, n);
  float* m = vector(1, n);

  float h = (b - a) / (n - 1);

  // Fill nodes and nodes value vectors
  xw[1] = x_min;
  xw[n] = x_max;
  yw[1] = fun(x_min);
  yw[n] = fun(x_max);
  for (int i = 2; i < n; ++i) {
    xw[i] = x_min + (i-1)*h;
    yw[i] = fun(xw[i]);
  }

  wyzM(xw, yw, m, n, alpha, betha);

  printResults(xw, yw, m, a, b, n, 0.01);

  free_vector(xw, 1, n);
  free_vector(yw, 1, n);
  free_vector(m, 1, n);
}


int main(void) {
  int nodeCount[3] = {5, 8, 21};
  float (*fun)(float) = fun1;
  //float (*fun)(float) = fun2;

  for (int k = 0; k < 3; ++k) 
    interpolate(fun, x_min, x_max, nodeCount[k]);

  //interpolate(fun1, x_min, x_max, 10);

  //return 0;
}
