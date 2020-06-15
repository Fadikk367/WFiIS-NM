#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/gauher.c"
#include "/opt/NR/numerical_recipes.c/gammln.c"
#include "/opt/NR/numerical_recipes.c/gaulag.c"
#include "/opt/NR/numerical_recipes.c/gauleg.c"

#define c1a M_PI/3.f
#define c2a -0.8700577f
#define c3a 2.f/13.f

float fun1(float x) {
	return 1/(x*sqrt(x*x - 1));
}

float fun2Hermit(float x) {
	return log(fabsf(x))/2.f;
}

float fun2Legendre(float x) {
	return log(x)*exp(-x*x);
}

float fun3(float x) {
	return sin(2*x)*exp(-2*x);
}

int main(void) {
	for (int n = 2; n <= 100; ++n) {
		float* x = vector(1, n);
		float* w = vector(1, n);

		float x1 = 1.f;
		float x2 = 2.f;

		gauleg(x1, x2, x, w, n);

		float c1 = 0.0;
		for (int i = 1; i <= n; ++i) {
			c1 += w[i]*fun1(x[i]);
		}

		float err = fabsf(c1 - c1a);

		printf("%d	%f\n", n, err);
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}

	printf("\n\n");

	for (int n = 2; n <= 100; n+=2) {
		float* x = vector(1, n);
		float* w = vector(1, n);
		
		gauher(x, w, n);

		float c2 = 0.0;
		for (int i = 1; i <= n; ++i) {
			c2 += w[i] * fun2Hermit(x[i]);
		}

		float err = fabsf(c2 - c2a);

		printf("%d	%f\n", n, err);
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}

    printf("\n\n");
    printf("\n\n");

    for (int n = 1; n <= 100; ++n) {
		float* x = vector(1, n);
		float* w = vector(1, n);
		
		float x1 = 0.f;
		float x2 = 5.f;

		gauleg(x1, x2, x, w, n);


		float c2 = 0.0;
		for (int i = 1; i <= n; ++i) {
			c2 += w[i]*fun2Legendre(x[i]);
		}

		float err = fabsf(c2 - c2a);

		printf("%d	%f\n", n, err);
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}

	printf("\n\n");
	printf("\n\n");

	for (int n = 2; n <= 20; ++n) {
		float* x = vector(1, n);
		float* w = vector(1, n);

		gaulag(x, w, n, 0.f);

		float c3 = 0.0;
		for (int i = 1; i <= n; ++i) {
			c3 += w[i]*fun3(x[i]);
		}

		float err = fabsf(c3 - c3a);

		printf("%d	%f\n", n, err);
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}

	return 0;
}