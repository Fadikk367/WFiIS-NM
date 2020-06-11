#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define a 0.0
#define b 1.0

#define n 8

void printMatrix(double m[][n+1]);
double fun(double x);
double calculateSimpson(double h, int N);
double calculateMilne(double h, int N);
double calculateRichardsonElement(double D[][n+1], int w, int k);

int main(void) {
    double D[n + 1][n + 1] = {{0}};
    FILE* file = fopen("out.dat", "w");
    for (int w = 0; w <= 8; ++w) {
        int N = pow(2, w+1);
        double h_w = (b - a)/N;
        D[w][0] = calculateSimpson(h_w, N/2-1);
        if (w >= 1) {
            for (int k = 1; k <= w; ++k) {
                D[w][k] = calculateRichardsonElement(D, w, k);
            }
        }
        fprintf(file, "%.10lf  %.10lf\n", D[w][0], D[w][w]);
    }
    printMatrix(D);
    fprintf(file, "\n\n");
    double M[n + 1][n + 1] = {{0}};
    for (int w = 0; w <= 8; ++w) {
        int N = pow(2, w+2);
        double h_w = (b - a)/N;
        M[w][0] = calculateMilne(h_w, N/4-1);
        if (w >= 1) {
            for (int k = 1; k <= w; ++k) {
                M[w][k] = calculateRichardsonElement(M, w, k);
            }
        }
        fprintf(file, "%.10lf  %.10lf\n", M[w][0], M[w][w]);
    }
    printMatrix(M);
    return 0;
}


void printMatrix(double m[][n+1]) {
    for (int i = 0; i < n+1; ++i) {
        for (int j = 0; j < n+1; ++j) {
            printf("%.10lf  ", m[i][j]);
        }
        printf("\n");
    }
}

double fun(double x) {
    double result =  log(x*x*x + 3*x*x + x + 0.1);
    result *= sin(18*x);
    return result;
}

double calculateSimpson(double h, int N) {
    double result = 0.0;
    for (int i = 0; i <= N; ++i) {
        double x = a + 2*i*h;
        result += h/3.0 * (fun(x) + 4*fun(x + h) + fun(x + 2*h));
    }
    return result;
}

double calculateMilne(double h, int N) {
    double result = 0.0;
    for (int i = 0; i <= N; ++i) {
        double x = a + 4*i*h;
        result += 4.0*h/90.0 * (7*fun(x) + 32*fun(x + h) + 12*fun(x + 2*h) + 32*fun(x + 3*h) + 7*fun(x + 4*h));
    }
    return result;
}

double calculateRichardsonElement(double D[][n+1], int w, int k) {
    return (pow(4.0, k)*D[w][k-1] - D[w-1][k-1]) / (pow(4.0, k) - 1.0);
}