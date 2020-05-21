#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define N 200
#define steps 100
#define it 20
#define x_max 10.0
#define x_min -10.0
#define y_max 10.0
#define y_min -10.0

struct Point2D {
    double x;
    double y;
};

typedef struct Point2D Point;

double d_rand (const double min , const double max) {
    double r = (double) rand()/RAND_MAX ; // Przedzial [0 , 1]
    r = r*(max - min) + min; // Przeskalowanie do [min , max]
    return r;
}

double fun(double x, double y) {
    double result = sin(x)*sin(y) - exp(-pow(x + M_PI/2.0, 2.0) - pow(y - M_PI/2.0, 2.0));
    return result;
}

int main(void) {
    srand(time(0));
    Point* voyagers = (Point*) calloc(N, sizeof(Point));
    FILE* voyagerFile = fopen("w0.dat", "w");
    FILE* snapshotFile = fopen("T.dat", "w");

    for (int i = 0; i < N; ++i) 
        voyagers[i].x = voyagers[i].y = 5.0;

    for (int i = 0; i <= it; ++i) {
        double T = 10.0 / pow(2, i);
        for (int k = 0; k < steps; ++k) {
            for (int j = 0; j < N; ++j) {
                double dx = d_rand(-1.0, 1.0);
                double dy = d_rand(-1.0, 1.0);

                double newPosX = voyagers[j].x + dx;
                double newPosY = voyagers[j].y + dy;

                if (newPosX > x_max) newPosX = x_max;
                else if (newPosX < x_min) newPosX = x_min;

                if (newPosY > y_max) newPosY = y_max;
                else if (newPosY < y_min) newPosY = y_min;

                double newValue = fun(newPosX, newPosY);
                double currentValue = fun(voyagers[j].x, voyagers[j].y);

                if (newValue < currentValue || d_rand(0.0, 1.0) < exp(-(newValue - currentValue)/T)) {
                    voyagers[j].x = newPosX;
                    voyagers[j].y = newPosY;
                }
            }
            fprintf(voyagerFile, "%g\n", fun(voyagers[0].x, voyagers[0].y));
        }
        if (i == 0 || i == 7 || i == 20) {
            for (int j = 0; j < N; ++j)
                fprintf(snapshotFile, "%g, %g\n", voyagers[j].x, voyagers[j].y);
            fprintf(snapshotFile, "\n\n");
        }
    }

    double minValue = fun(voyagers[0].x, voyagers[0].y);
    int minValueIndex = 0;
    for (int i = 1; i < N; ++i) {
        double value = fun(voyagers[i].x, voyagers[i].y);
        if (value < minValue) {
            minValue = value;
            minValueIndex = i;
        }
    }

    printf("\nGlobal minimum: f(%g, %g) = %g\n", voyagers[minValueIndex].x, voyagers[minValueIndex].y, minValue);

    free(voyagers);
    return 0;
}
