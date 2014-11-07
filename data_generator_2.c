#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>

int i, j;

void stampaMatriceSuFile(int righe, int col, double mat[righe][col], FILE *punt_file) {
    int i, j;
    for (i = 0; i < righe; i++) {
        for (j = 0; j < col; j++) {
            fprintf(punt_file, "%lf", mat[i][j]);
            fprintf(punt_file, " ");
        }
        fprintf(punt_file, "\n");
    }
    fflush(punt_file);
}

double uni_rand() /* distribuzione uniforme, (0..1] */ {
    return (lrand48() + 1.0) / (RAND_MAX + 1.0);
}

double random_normal() {
    /* distribuzione normale, centrata su 0, std dev 1 */
    return sqrt(-2 * log(uni_rand())) * cos(2 * M_PI * uni_rand());
}

int main() {
    FILE *out_X;
    out_X = fopen("x.dat", "w");
    int n, d;
    double sigma_gauss;
    printf("numero di punti: ");
    scanf("%d", &n); //numero di punti totale in input
    printf("numero di dimensioni: ");
    scanf("%d", &d);
    /*printf("mi gauss: ");
    scanf("%lf", &mi_gauss);
    printf("sigma gauss: ");
    scanf("%lf", &sigma_gauss);*/
    sigma_gauss = 2;

    double X[n][d];
    //INIT X
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            //gaussiana con media mi_gauss e devstd sigma_gauss
            srand48(rand());
            if (j == 0 && i < 100)//x1
                X[i][j] = 0 + (sigma_gauss * random_normal());
            else if (j == 0 && i >= 100 && i < 200)//x2
                X[i][j] = 10 + (sigma_gauss * random_normal());
            else if (j == 0 && i >= 200)//x3
                X[i][j] = 20 + (sigma_gauss * random_normal());

            if (j == 1 && i < 100)//y1
                X[i][j] = 10 + (4 * sigma_gauss * random_normal());
            else if (j == 1 && i >= 100 && i < 200)//y2
                X[i][j] = 30 + (4 * sigma_gauss * random_normal());
            else if (j == 1 && i >= 200)//y3
                X[i][j] = 20 + (4 * sigma_gauss * random_normal());
        }
    }
    stampaMatriceSuFile(n, d, X, out_X);
    puts("stampata su x.dat");
    return (0);
}
