#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>

#define c 4 //numero di centri di cluster
#define n 200 //numero di punti totale in input
double m = 2.0; //fuzzification
#define d 2 //dimensioni spaziali
#define moltiplicatore_popolazione = 20;

double CR = 0.8; //crossover rate [0,1]
double XB = DBL_MAX; //XB index

int numero_generazioni = 200;
int conteggio_crossover = 0;

double X[n][d]; //dati input
double U[c][n]; //partition matrix
double V[c][d]; //matr centroidi

double POP[c * 20][d];

void stampaMatrice(int righe, int col, double mat[righe][col]) {
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

double calcolaXB() {
    /*
     XB funzione del rapporto fra la variazione totale sigma e la separazione minima
     fra i centroidi
     */
    
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
                dist_tmp = pow(calcDistanza(V[i], V[j]),2.0);
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

double fRand(double fMin, double fMax) {
    double f = (double) rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double drand() /* distribuzione uniforme, (0..1] */ {
    return (rand() + 1.0) / (RAND_MAX + 1.0);
}

double random_normal() {
    /* distribuzione normale, centrata su 0, std dev 1 */
    return sqrt(-2 * log(drand())) * cos(2 * M_PI * drand());
}

int main(int argc, char** argv) {
    FILE *out_V, *out_X, *out_U;
    out_V = fopen("v.dat", "w");
    out_X = fopen("x.dat", "w");
    out_U = fopen("u.dat", "w");

    double coordXCentroidiAttese[c];

    int i, j;
    int conteggio_selezioni = 0;
    //INIT X
    int mi_gauss = 1;
    double sigma_gauss = 1.0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++)
            X[i][j] = mi_gauss + (sigma_gauss * random_normal()); //gaussiana con media mi_gauss e devstd sigma_gauss
        if (i == 0)
            coordXCentroidiAttese[0] = mi_gauss;
        if (i == 50) {
            mi_gauss *= 8;
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
    stampaMatriceSuFile(n, d, X, out_X);

    puts("ATTESA");
    sleep(3);
    //INIT V
    //V viene inizializzata con alcuni dei punti di input (spostati)
    srand48(time(0));
    for (i = 0; i < c; i++)
        for (j = 0; j < d; j++)
            V[i][j] = X[random_at_most(n)][random_at_most(d)] + 1;
    //V[i][j] = 10 * drand48() - 5;

    puts("\ninizializzazione matrice V:");
    stampaMatrice(c, d, V);

    //INIT POP
    for (i = 0; i < c * 20; i++) {
        for (j = 0; j < d; j++)
            POP[i][j] = X[random_at_most(n)][random_at_most(d)];
    }
    puts("\ninizializzazione matrice POP:");
    stampaMatrice(c * 20, d, POP);
    printf("\n########## FINE INIT #############\n");

    //CALCOLO PARTITION MATRIX (INIT)
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

    //f è un numero fra 0 e 2
    double f = fRand(0.0, 2.0);
    double prob_crossover;

    do {
        printf("\n\nCOUNTDOWN GENERAZIONE: %d", numero_generazioni);
        //DIFFERENTIAL EVOLUTION (SINGOLO PASSO)
        double trial[d], temp[d];
        int indiceTarget = 0;
        //crea un nuovo mutante per ogni vettore della popolazione
        while (indiceTarget < c) {
            //i tre vettori devono essere scelti a caso nella popolazione
            //diversi dal target e mutualmente
            int indice_1, indice_2, indice_3;
            do {
                indice_1 = random_at_most((long) c * 20);
            } while (indice_1 == indiceTarget);

            do {
                indice_2 = random_at_most((long) c * 20);
            } while (indice_2 == indiceTarget || indice_2 == indice_1);

            do {
                indice_3 = random_at_most((long) c * 20);
            } while (indice_3 == indiceTarget || indice_3 == indice_1 || indice_3 == indice_2);

            if (indice_1 == indiceTarget || indice_2 == indiceTarget || indice_3 == indiceTarget)
                puts("INDICE NON VALIDO!!!!");

            //MUTAZIONE E INCROCIO
            for (i = 0; i < d; i++) {
                prob_crossover = fRand(0.0, 1.0);
                if (prob_crossover < CR) {
                    trial[i] = POP[indice_3][i] + f * (POP[indice_1][i] - POP[indice_2][i]);
                    conteggio_crossover++;
                } else
                    trial[i] = V[indiceTarget][i];
            }


            //copia vettore originale in temp
            copiaVettore(d, V[indiceTarget], temp);
            //inserimento temporaneo del mutato
            copiaVettore(d, trial, V[indiceTarget]);

            //SELEZIONE
            //CALCOLO XB
            double new_XB = calcolaXB();
            if (new_XB < XB) {//XB è migliorato, tengo il mutante e ricalcolo tutto U
                printf("\nNUOVO XB: %lf\n", new_XB);
                XB = new_XB;
                conteggio_selezioni++;

                //RICALCOLO PARTITION MATRIX
                for (i = 0; i < c; i++) {
                    for (j = 0; j < n; j++) {
                        double esponente = 2.0 / (m - 1.0);
                        double denom = 0.0;
                        double dist_x_j__v_i = calcDistanza(X[j], V[i]);

                        int k;
                        for (k = 0; k < c; k++) {
                            double dist_xj_vk = calcDistanza(X[j], V[k]);
                            denom += pow((dist_x_j__v_i / dist_xj_vk), esponente);
                        }
                        U[i][j] = 1.0 / denom;
                    }
                }
            } else {// XB NON è migliorato
                //ripristino del vettore target
                copiaVettore(d, temp, V[indiceTarget]);
            }
            //END CALCOLO XB
            indiceTarget++;
        }
        //END DE
        numero_generazioni--;
    } while (numero_generazioni > 0);

    puts("***********************************************");
    puts("matrice X:");
    //stampaMatrice(n, d, X);
    puts("matrice U FINALE");
    //stampaMatrice(c, n, U);
    stampaMatriceSuFile(c, n, U, out_U);
    puts("matrice V FINALE");
    stampaMatrice(c, d, V);
    stampaMatriceSuFile(c, d, V, out_V);
    puts("\ninizializzazione matrice POP:");
    stampaMatrice(c * 20, d, POP);
    printf("\nultimo XB:%lf", XB);
    printf("\nconteggio_crossover:%d\n", conteggio_crossover);
    printf("\nCR: %lf\n", CR);
    printf("\nconteggio selezioni:%d\n", conteggio_selezioni);
    stampaVett(c, coordXCentroidiAttese);
    puts("***********************************************");


    return (EXIT_SUCCESS);
}

