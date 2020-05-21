#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 201
#define x_min -4.0
#define x_max 4.0
#define x_0 2.0
#define sigma ((x_max - x_min) / 16)

double generateNoise() {
    double Y = ((double) rand()) / (RAND_MAX + 1.0);
    return (Y - 0.5) / 5.0;
}


double f(double x) {
    double result = sin(14*M_PI*x / (x_max - x_min));
    result *= (exp(-pow(x - x_0, 2) / (2*sigma*sigma)) + exp(-pow(x + x_0, 2) / (2*sigma*sigma)));
    return result;
}

double calculateAlpha(double* phi, double* x) {
    double numerator = 0.0;
    double denominator = 0.0;
    for (int i = 0; i < N; ++i) {
        numerator += x[i] * phi[i] * phi[i];
        denominator += phi[i] * phi[i];
    }
    return numerator / denominator;
}

double calculateBetha(double* prevPhi, double* phi, double* x) {
    double numerator = 0.0;
    double denominator = 0.0;
    for (int i = 0; i < N; ++i) {
        numerator += x[i] * prevPhi[i] * phi[i];
        denominator += prevPhi[i] * prevPhi[i];
    }
    return numerator / denominator;
}

double calculateC(double* phi, double* y) {
    double result = 0.0;
    for (int i = 0; i < N; ++i)
        result += y[i]*phi[i];

    return result;
}

double calculateS(double* phi) {
    double result = 0.0;
    for (int i = 0; i < N; ++i)
        result += phi[i]*phi[i];

    return result;
}

int main(void) {
    srand(time(0));

    // declare and initialize output files
    FILE *approxFile, *nodesFile, *gramFile;
    approxFile = fopen("approx.dat", "w");
    nodesFile = fopen("pkt.dat", "w");
    gramFile = fopen("Gram.dat", "w");

    // allocate memory
    double* x = (double*) calloc(N, sizeof(double));
    double* y = (double*) calloc(N, sizeof(double));

    double** phi = (double**) calloc(51, sizeof(double*));
    double** F = (double**) calloc(3, sizeof(double*));

    double* c = (double*) calloc(51, sizeof(double));
    double* s = (double*) calloc(51, sizeof(double));

    for (int i = 0; i < 51; ++i)
        phi[i] = (double*) calloc(N, sizeof(double));

    for (int i = 0; i < 3; ++i)
        F[i] = (double*) calloc(N, sizeof(double));

    int mNumber[3] = {10, 30, 50};


    // fill nodes vectors and apply noise
    double h = (x_max - x_min) / (N - 1);
    for (int i = 0; i < N - 1; ++i) {
        x[i] = x_min + i*h;
        //y[i] = f(x[i]) + generateNoise();
        y[i] = f(x[i]);
        fprintf(nodesFile, "%g  %g\n", x[i], y[i]);
    }
    x[N - 1] = x_max;
    //y[N - 1] = f(x[N]) + generateNoise();
    y[N - 1] = f(x[N]);
    fprintf(nodesFile, "%g  %g\n", x[N - 1], y[N - 1]);

    // fill phi matrix
    for (int k = 0; k < N; ++k)
        phi[0][k] = 1.0;

    double alpha = calculateAlpha(phi[0], x);
    for (int k = 0; k < N; ++k) 
        phi[1][k] = (x[k] - alpha) * phi[0][k];

    for (int j = 2; j < 51; ++j) {
        double alpha = calculateAlpha(phi[j - 1], x);
        double betha = calculateBetha(phi[j-  2], phi[j - 1], x);
        for (int k = 0; k < N; ++k) {
            phi[j][k] = (x[k] - alpha)*phi[j - 1][k] - betha*phi[j - 2][k];
        }
    }    

	for(int i = 0; i < N; ++i) {
		fprintf(gramFile, "%g  ", x[i]);
		for(int j = 0; j < 7; ++j)
			fprintf(gramFile, "%g  ", phi[j][i] / phi[j][0]);
		fprintf(gramFile, "\n");
    }

    // calculate c and s coefficients
    for (int j = 0; j < 51; ++j) {
        c[j] = calculateC(phi[j], y);
        s[j] = calculateS(phi[j]);
    }

    // approximation
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < N; ++k) {
            double result = 0.0;
            int m = mNumber[i];
            for (int j = 0; j < m; ++j)
                result += phi[j][k] * c[j]/s[j];

            F[i][k] = result;
            fprintf(approxFile, "%g  %g\n", x[k], result);
        }
        fprintf(approxFile, "\n\n");
    }

    // clear memory and resources
	fclose(nodesFile);	
	fclose(gramFile);	
	fclose(approxFile);	

    free(x);
    free(y);
    free(c);
    free(s);

    for (int i = 0; i < 51; ++i)
        free(phi[i]);
    free(phi);

    for (int i = 0; i < 3; ++i)
        free(F[i]);
    free(F);

    return 0;
}