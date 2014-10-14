#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#define c 4 //numero di centri di cluster
#define n 2000 //numero di punti totale in input
double m = 2.0; //fuzzification
#define d 2 //dimensioni spaziali
double epsilon = 0.001; //minima distanza per arrestare
double distanze[c]; //vettore con le dist fra centroidi dopo l'aggiornamento

int mi_gauss = 2;
double sigma_gauss = 2.0;


double X[n][d]; //dati input
double U[c][n]; //partition matrix
double V[c][d]; //matr centroidi
double max;
double coordXCentroidiAttese[c];
void stampaMatrice(int righe, int col, double mat[][col]) {
    int i, j;
    printf("\n\n");
    for (i = 0; i < righe; i++) {
        for (j = 0; j < col; j++) {
            printf("%lf", mat[i][j]);
            printf(" ");
        }
        puts("");
    }
}

void stampaMatriceSuFile(int righe, int col, double mat[righe][col], FILE *punt_file) {
    int i, j;
    for (i = 0; i < righe; i++) {
        for (j = 0; j < col; j++) {
            fprintf(punt_file, "%lf", mat[i][j]);
            fprintf(punt_file, " ");
        }
        fprintf(punt_file, "\n");
    }
}

void stampaVett(double x[c]) {
    puts("____");
    printf("|%lf|\n", x[0]);
    printf("|%lf|\n", x[1]);
    puts("____");
}

double calcDistanza(double a[d], double b[d]) {
    double ris = 0;
    int i;
    for (i = 0; i < d; i++)
        ris += pow(a[i] - b[i], 2);
    return sqrt(ris);
}

void prodottoScalareVettore(double scal, double vett_in[], double vett_out[]) {

    int i;
    for (i = 0; i < d; i++) {
        vett_out[i] = vett_in[i] * scal;
    }
}

double maxDistCentroidi() {
    double max = 0.0;
    int i;
    for (i = 0; i < c; i++) {
        if (distanze[i] > max)
            max = distanze[i];
    }
    return max;
}

double drand() /* distribuzione uniforme, (0..1] */ {
    return (rand() + 1.0) / (RAND_MAX + 1.0);
}

double random_normal() {
    /* distribuzione normale, centrata su 0, std dev 1 */
    return sqrt(-2 * log(drand())) * cos(2 * M_PI * drand());
}

int main(int argc, char** argv) {
    //PUNTATORI A FILE DI OUTPUT
    FILE *out_V, *out_X, *out_U;
    out_V = fopen("v.dat", "w");
    out_X = fopen("x.dat", "w");
    out_U = fopen("u.dat", "w");
    /*
     * INPUT: X, c, m
     * OUTPUT: U,V
     */  
    int i,j;
    //INIT X
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++)
            //gaussiana con media mi_gauss e devstd sigma_gauss
            X[i][j] = mi_gauss + (sigma_gauss * random_normal());
        if (i == 0)
            coordXCentroidiAttese[0] = mi_gauss;
        if (i == 50) {
            mi_gauss *= 4;
            coordXCentroidiAttese[1] = mi_gauss;
        }
        if (i == 100) {
            mi_gauss *= 2;
            coordXCentroidiAttese[2] = mi_gauss;
        }
        if (i == 150) {
            mi_gauss *= 2;
            coordXCentroidiAttese[3] = mi_gauss;
        }
    }
    puts("matrice X:");
        stampaMatrice(n, d, X);
        puts("");
    
    //INIT V
    srand48(3);
    for (i = 0; i < c; i++)
        for (j = 0; j < d; j++)
            V[i][j] = 10 * drand48() - 5;

    puts("\ninizializzazione matrice V:");
    stampaMatrice(c, d, V);
    printf("#######################\n\n");

    int contPassi = 0;
    max = 0.0;
    do {
        //CALCOLO PARTITION MATRIX
        printf("PASSO: %d\n\n", ++contPassi);
        for (i = 0; i < c; i++) {
            for (j = 0; j < n; j++) {
                double esponente = 2.0 / (m - 1.0);
                double denom = 0.0;
                double dist_x_j__v_i = calcDistanza(X[j], V[i]);

                int k; //SOMMATORIA 1
                for (k = 0; k < c; k++) {
                    double dist_xj_vk = calcDistanza(X[j], V[k]);
                    denom += pow((dist_x_j__v_i / dist_xj_vk), esponente);
                }
                U[i][j] = 1.0 / denom;
            }
        }


        //RICALCOLO POSIZIONE CENTROIDI
        int i, j, z, k;
        double old[d];
        double denom;
        for (i = 0; i < c; i++) {
            for (z = 0; z < d; z++)
                old[z] = V[i][z]; //per confronto diff
            denom = 0.0;
            for(j=0; j<n;j++)//sommatoria denom (fatta una sola volta a centr.)
                denom += pow(U[i][j], m);
            for (k = 0; k < d; k++) {
                double num = 0.0;
                for (j = 0; j < n; j++) {//SOMMATORIA numeratore
                    num += X[j][k] * pow(U[i][j], m);
                }
                V[i][k] = num / denom;
            }
            distanze[i] = pow(calcDistanza(V[i], old), 2.0);
        }



        
        puts("matrice U:");
        stampaMatrice(c, n, U);
        puts("");
        puts("matrice V:");
        stampaMatrice(c, d, V);
        puts("");
        stampaMatriceSuFile(c,d,V,out_V);
        stampaMatriceSuFile(c,n,U,out_U);
        stampaMatriceSuFile(n,d,X,out_X);
        printf("MAXDISTCENTROIDI: %lf", maxDistCentroidi());
        printf("\n\n\n\n\n");
    } while (maxDistCentroidi() > epsilon);
    return (0);
}
