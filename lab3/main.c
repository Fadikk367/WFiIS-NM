#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// some constans
#define N 2000
#define F_0 0.1
#define Beta 0.4
#define Omega 0.8
#define h 0.02

// initial conditions
#define V_0  0
#define x_0  1.0
#define w  1.0

// helper function
float sum_squares(const float* const v, int len) {
  float result = 0;
  for (int i = 0; i < len; i++) {
    result += v[i]*v[i];
  }
  return result;
}

int main(void) {

  // allocate resources
  float* b = (float*) calloc(sizeof(float), (N + 1));
  float* d_0 = (float*) calloc(sizeof(float), (N + 1));
  float* d_1 = (float*) calloc(sizeof(float), (N + 1));
  float* d_2 = (float*) calloc(sizeof(float), (N + 1));
  float* x_s = (float*) malloc(sizeof(float)*(N + 1));
  float* x_n = (float*) calloc(sizeof(float), (N + 1));

  // calculate coefficients
  float a_1 = 1.0;
  float a_2 = w*w * h*h - 2 - Beta*h;
  float a_3 = 1 + Beta*h;

  // fill b and d-s vectors
  b[0] = 1;
  b[1] = 0;

  d_0[0] = 1;
  d_0[1] = 1;

  d_1[0] = 0;
  d_1[1] = -1;

  d_2[0] = 0;
  d_2[1] = 0;
  
  for (int i = 2; i < N + 1; i++) {
    b[i] = F_0 * sin(Omega*(h * i)) * h*h;
    d_0[i] = a_3;
    d_1[i] = a_2;
    d_2[i] = a_1;
  }

  // main loop - calculating next x_n vector
  int it = 0;
  while (it < 1e5) {
    ++it;
    // avoiding negative indexes
    x_n[0] = 1/d_0[0] * b[0];
    x_n[1] = (1/d_0[1]) * (b[1] - d_1[1]*x_n[0]);

    for (int i = 2; i < N + 1; i++) 
      x_n[i] = (1/d_0[i]) * (b[i] - d_1[i]*x_n[i - 1] - d_2[i]*x_n[i - 2]);

    float s_n = sum_squares(x_n, N + 1);
    float s_s = sum_squares(x_s, N + 1);

    // convergence test
    if (fabs(s_n - s_s) < 1e-6) break;

    for (int i = 0; i < N + 1; i++) 
      x_s[i] = x_n[i];
  }

  // results printing
  for (int i = 0; i < N + 1; i++) 
    printf("%.2f, %f \n", i*h, x_n[i]);
  printf("\n\n");

  printf("iterations to cevergence: %d", it);

  // memory clear
  free(b);
  free(d_0);
  free(d_1);
  free(d_2);
  free(x_s);
  free(x_n);

  return 0;
}