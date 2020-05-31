#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_errno.h>

#define T 1.0
#define T_MAX 3.0*T
#define omega 2.0*M_PI/T
#define sigma T/20.0

double noise() {
	double result = (double) rand()/(RAND_MAX + 1.0);
	return result - 0.5;
}

double f0(double t) {
	return sin(omega*t) + sin(2*omega*t) + sin(3*omega*t);
}

double g0(double t) {
	double result = 1.0/(sigma*sqrt(2*M_PI));
	result *= exp(-t*t/(2*sigma*sigma));
	return result;
}


int main(void) {
	srand(time(0));
	int kNums[] = {8, 10, 12};		
	FILE* file8 = fopen("k8.dat", "w");
	FILE* file10 = fopen("k10.dat", "w");
	FILE* file12 = fopen("k12.dat", "w");

	for (int p = 0; p < 3; ++p) {
		int k = kNums[p];
		int N = pow(2, k);
		double dt = T_MAX / N;

		FILE* file = NULL;
		switch(k) {
			case 8:
				file = file8;
				break;
			case 10:
				file = file10;
				break;
			case 12:
				file = file12;
				break;
			default:
			break;
		}

		double* f = (double*) calloc(2*N, sizeof(double));
		double* g = (double*) calloc(2*N, sizeof(double));
		double* g1 = (double*) calloc(2*N, sizeof(double));
		double* g2 = (double*) calloc(2*N, sizeof(double));

		for (int i = 0; i < N - 1; ++i) {
			double t = i * dt;
			f[2*i] = f0(t) + noise();
			f[2*i + 1] = 0.0; 

			g1[2*i] = g0(t);
			g1[2*i + 1] = 0.0; 

			g2[2*i] = g0(t);
			g2[2*i + 1] = 0.0; 

			fprintf(file, "%g %g\n", t, f[2*i]);
		}


		gsl_fft_complex_radix2_forward(f, 1, N);
		gsl_fft_complex_radix2_forward(g1, 1, N);
		gsl_fft_complex_radix2_backward(g2, 1, N);

		double a1, a2, b1, b2;
		for (int i = 0; i < N - 1; ++i) {
			g[2*i] = g1[2*i] + g2[2*i];
			g[2*i + 1] = g1[2*i + 1] + g2[2*i + 1];

			a1 = f[2*i];
			b1 = f[2*i+1];
			a2 = g[2*i];
			b2 = g[2*i+1];

			f[2*i] = a1 * a2 - b1 * b2;
			f[2*i+1] = a1 * b2 + a2 * b1;
		}

		gsl_fft_complex_radix2_backward(f, 1, N);

		double fmax = -999999999.;
		for (int i = 0; i < N - 1; ++i)
			if (abs(f[2*i]) > fmax) fmax = abs(f[2*i]);
	
        printf("f_max(k = %d) = %g\n", k, fmax);
		fprintf(file, "\n\n");
		for (int i = 0; i < N - 1; ++i)
			fprintf(file, "%g  %g\n", dt*i, f[2*i]*2.5/fmax);

		free(f);
		free(g);
		free(g1);
		free(g2);
	}

	return 0;
}