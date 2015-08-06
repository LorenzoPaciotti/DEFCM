#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>

int c, n, d;

double m = 2.0; //fuzzification
double epsilon = 0.000000000001; //minima distanza per arrestare


int i, j;
double **X; //dati input
double **U; //partition matrix
double **V; //matr centroidi
double max;

int attivaGnuPlot = 0;

void stampaMatrice(int righe, int col, double **mat) {
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

void stampaMatriceSuFile(int righe, int col, double **mat, FILE *punt_file) {
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

double calcDistanza(double a[d], double b[d]) {
    double ris = 0;
    int i;
    for (i = 0; i < d; i++)
        ris += pow(a[i] - b[i], 2.0);
    return sqrt(ris);
}

double maxDistCentroidi(double distanze[c]) {
    double max = 0.0;
    int i;
    for (i = 0; i < c; i++) {
        if (distanze[i] > max)
            max = distanze[i];
    }
    return max;
}

double calcolaXB(double **V, double **U, int debug) {
    /*
     XB funzione del rapporto fra la variazione totale sigma
     e la separazione minima fra i centroidi
     */
    //if(debug == 1)
    //    puts("debug XB");
    //CALCOLO MIN_SEP
    double min_sep = DBL_MAX;
    int i, j;
    double dist_tmp = 0;
    j = 0;
    for (i = 0; i < c; i++) {
        if (j == i)
            j++;
        while (j < c) {
            if (j == i)
                j++;
            if (j < c) {
                dist_tmp = pow(calcDistanza(V[i], V[j]), 2.0);
                if (dist_tmp < min_sep)
                    min_sep = dist_tmp;
                j++;
            }
        }
        j = 0;
    }

    //CALCOLO SIGMA
    double sigma = 0.0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < c; j++) {
            sigma += pow(U[j][i], m) * pow(calcDistanza(V[j], X[i]), 2.0);
        }
    }
    return sigma / (n * min_sep);
}

void plot() {
    if ((d == 2 || d == 3) && attivaGnuPlot) {
        char *commandsForGnuplot[] = {"set key off", "set term x11 1", "set title \"matrice X\"", "", "set term x11 2", "set key off", "set title \"FCM - matrice V\"", ""};
        if (d == 2) {
            commandsForGnuplot[3] = "plot 'x.dat' pointtype 3";
            commandsForGnuplot[7] = "plot 'v_fcm.dat' pointtype 3";
        } else if (d == 3) {
            commandsForGnuplot[3] = "splot 'x.dat' pointtype 3";
            commandsForGnuplot[7] = "splot 'v_fcm.dat' pointtype 3";
        }

        FILE * gnuplotPipe = popen("gnuplot -persistent", "w");
        int i;
        for (i = 0; i < 8; i++) {
            fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);
            fflush(gnuplotPipe);
        }
    }
}

int main(int argc, char** argv) {
    //PUNTATORI A FILE DI OUTPUT
    FILE *out_V, *out_X, *out_U;
    out_V = fopen("v_fcm.out", "w");
    out_U = fopen("u_fcm.out", "w");
    out_X = fopen("dataset/gauss4.data", "r");

    //letture da utente
    puts("numero di punti in input");
    scanf("%d", &n);
    puts("numero di dimensioni");
    scanf("%d", &d);
    puts("numero di centroidi");
    scanf("%d", &c);

    double distanze[c]; //vettore con le dist fra centroidi dopo l'aggiornamento

    //allocazione X
    int row;
    X = malloc(n * sizeof (double*));
    for (row = 0; row < n; row++) {
        X[row] = malloc(d * sizeof (double));
    }

    //lettura X da file
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            if (!fscanf(out_X, "%lf", &X[i][j]))
                break;
        }
    }
    fclose(out_X);

    //alloc V
    V = malloc(c * sizeof (double*));
    for (row = 0; row < c; row++) {
        V[row] = malloc(d * sizeof (double));
    }

    //INIT V
    srand48(rand());
    for (i = 0; i < c; i++)
        for (j = 0; j < d; j++)
            V[i][j] = 10 * drand48() - 5;

    //alloc U
    U = malloc(c * sizeof (double*));
    for (row = 0; row < c; row++) {
        U[row] = malloc(n * sizeof (double));
    }


    printf("########END INIT##########\n\n");
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

    } while (maxDistCentroidi(distanze) > epsilon);

    puts("matrice V:");
    stampaMatrice(c, d, V);
    stampaMatriceSuFile(c, d, V, out_V);
    stampaMatriceSuFile(c, n, U, out_U);
    printf("MAXDISTCENTROIDI: %lf\n", maxDistCentroidi(distanze));
    puts("indice XB:");
    double xb = calcolaXB(V, U, 0);
    printf("%lf\n", xb);
    plot();
    return (0);
}
