#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>

#define c 4 //numero di centri di cluster
#define n 200 //numero di punti totale in input
#define d 4 //dimensioni spaziali

double m = 2.0; //fuzzification
double epsilon = 0.001; //minima distanza per arrestare
double distanze[c]; //vettore con le dist fra centroidi dopo l'aggiornamento

double X[n][d]; //dati input
double U[c][n]; //partition matrix
double V[c][d]; //matr centroidi
double max;

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

double calcDistanza(double a[d], double b[d]) {
    double ris = 0;
    int i;
    for (i = 0; i < d; i++)
        ris += pow(a[i] - b[i], 2);
    return sqrt(ris);
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

int main(int argc, char** argv) {
    //PUNTATORI A FILE DI OUTPUT
    FILE *out_V, *out_X, *out_U;
    out_V = fopen("v_fcm.dat", "w");
    out_U = fopen("u_fcm.dat", "w");
    int i, j;
    
    //lettura X
    out_X = fopen("x.dat", "r");
    i = 0; j = 0;
    //lettura X da file
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            if (!fscanf(out_X, "%lf", &X[i][j]))
                break;
        }
    }
    fclose(out_X);

    //INIT V
    srand48(3);
    for (i = 0; i < c; i++)
        for (j = 0; j < d; j++)
            V[i][j] = 10 * drand48() - 5;

    puts("\ninizializzazione matrice V:");
    stampaMatrice(c, d, V);
    printf("#######################\n\n");
    max = 0.0;
    do {
        //CALCOLO PARTITION MATRIX
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
            for (j = 0; j < n; j++)//sommatoria denom (fatta una sola volta a centr.)
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
        
    } while (maxDistCentroidi() > epsilon);

    //puts("matrice U:");
    //stampaMatrice(c, n, U);
    puts("");
    puts("matrice V:");
    stampaMatrice(c, d, V);
    puts("");
    stampaMatriceSuFile(c, d, V, out_V);
    stampaMatriceSuFile(c, n, U, out_U);
    printf("MAXDISTCENTROIDI: %lf", maxDistCentroidi());
    printf("\n\n\n\n\n");
    return (0);
}
