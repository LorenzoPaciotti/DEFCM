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
    double mi_gauss, sigma_gauss;
    printf("numero di punti: ");
    scanf("%d", &n); //numero di punti totale in input
    printf("numero di dimensioni: ");
    scanf("%d", &d);
    printf("sigma gauss: ");
    scanf("%lf", &sigma_gauss);

    double X[n][d];
    //INIT X
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            //gaussiana con media mi_gauss e devstd sigma_gauss
            srand48(rand());
            //if(j%2 == 0)
                X[i][j] = mi_gauss + (sigma_gauss * random_normal());
            //else
                //X[i][j] = -mi_gauss + (sigma_gauss * random_normal());
        }
        if (i == 50) {
            mi_gauss = 8;
        }
        if (i == 100) {
            mi_gauss = -8;
        }
        if (i == 150) {
            mi_gauss = 2;
        }
        if (i == 200) {
            mi_gauss = 4;
        }


    }
    stampaMatriceSuFile(n, d, X, out_X);
    puts("stampata su x.dat");
    return (0);
}
