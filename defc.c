#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#define c 8 //numero di centri di cluster
#define n 32 //numero di punti totale in input
double m = 2.0; //fuzzification
#define d 2 //dimensioni spaziali
double epsilon = 0.001; //minima distanza per arrestare
double distanze[c]; //vettore con le dist fra centroidi dopo l'aggiornamento

double CR = 0.4; //crossover rate [0,1]
double XB = DBL_MAX; //XB index

double X[n][d]; //dati input
double U[c][n]; //partition matrix
double V[c][d]; //matr centroidi

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

void stampaVett(int dim, double x[dim]) {
    int i;
    for (i = 0; i < dim; i++) {
        printf("|%lf|\n", x[i]);
    }
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

void copiaVettore(int dim, double input[dim], double output[dim]) {
    int i;
    for (i = 0; i < dim; i++) {
        output[i] = input[i];
    }
}

int calcolaXB() {
    /*
     XB è la funzione del rapporto fra la variazione totale sigma e la separazione minima
     fra i centroidi
     */

    double sep = DBL_MAX;

    int i, j;
    for (i = 0; i < c; i++) {
        double dist = calcDistanza(V[i], V[i + 1]);
        if (dist < sep)
            sep = dist;
    }

    double sigma = 0.0;
    for (i = 0; i < c; i++) {
        for (j = 0; j < n; j++) {
            sigma += U[i][j] * pow(calcDistanza(V[j], X[i]), 2.0);
        }
    }


    return sigma / (n * sep);
}

long random_at_most(long max) {
    unsigned long
    // max <= RAND_MAX < ULONG_MAX, so this is okay.
    num_bins = (unsigned long) max + 1,
            num_rand = (unsigned long) RAND_MAX + 1,
            bin_size = num_rand / num_bins,
            defect = num_rand % bin_size;

    long x;
    // This is carefully written not to overflow
    while (num_rand - defect <= (unsigned long) (x = random()));

    // Truncated division is intentional
    return x / bin_size;
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

/*la fitness di un vettore è migliore di quella di un altro se la sua introduzione
 * porta ad avere un indice Xie-Beni di bontà del clustering migliore
 */

int main(int argc, char** argv) {


    int i, j;
    //INIT X
    srand48(time(0));
    for (i = 0; i < n; i++)
        for (j = 0; j < d; j++)
            X[i][j] = 10 * drand48() + 1;
    //INIT U
    for (i = 0; i < c; i++)
        for (j = 0; j < n; j++)
            U[i][j] = X[i][j];
    puts("INIT matrice X:");
    stampaMatrice(n, d, X);
    puts("");
    puts("INIT matrice U:");
    stampaMatrice(c, n, U);
    puts("");

    sleep(1);
    //INIT V
    srand48(time(0));
    for (i = 0; i < c; i++)
        for (j = 0; j < d; j++)
            V[i][j] = 10 * drand48() - 5;
    
    sleep(1);
    


    double g1[d], g2[d], g3[d], trial[d];
    int indiceTarget = 0;
    while (indiceTarget < c) {
        //test, i tre vettori devono essere scelti a caso nella popolazione
        copiaVettore(d, V[random_at_most((long)c)], g1);
        copiaVettore(d, V[random_at_most((long)c)], g2);
        copiaVettore(d, V[random_at_most((long)c)], g3);

        //MUTAZIONE
        //f è un numero fra 0 e 1
        double f = fRand(0.0,1.0);
        int i;
        for (i = 0; i < d; i++) {
            trial[i] = g3[i] + f * (g1[i] - g2[i]);
        }

        //INCROCIO

        //CALCOLO XB
        double new_XB = calcolaXB();
        if (calcolaXB() < XB) {
            puts("##################");
            printf("NUOVO XB: %lf",new_XB);
            puts("##################");
            XB = new_XB;
            copiaVettore(d, trial, V[indiceTarget]); //il nuovo vett sostituisce
        }


        printf("\n\nINDICE TARGET:%d\n", indiceTarget);
        puts("matrice U");
        stampaMatrice(c, n, U);
        puts("matrice V");
        stampaMatrice(c, d, V);

        indiceTarget++;
    }

    return (EXIT_SUCCESS);
}

