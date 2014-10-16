#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#include <string.h>

FILE *out_V, *out_X, *out_U, *out_LOG, *out_LOG_RIS; //puntatori a file di output

double **X;
double m = 2.0; //fuzzification
double CR = 0.9; //crossover rate [0,1]

#define num_pop 100

typedef struct el_pop {//individuo della popolazione
    double **V_p;
    double **U_p;
    double indice_xb;
} el_pop;

el_pop *POP_NEW[num_pop];
el_pop *POP_NOW[num_pop];

//attiva e disattiva GnuPlot
int attivaGnuPlot = 1;



int num_pop_iniziale, numero_generazioni, i, j, k, pop_index;
double esponente_U;
double xb_selezionato;
int n, c, d;

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
        ris += pow(a[i] - b[i], 2);

    if (ris == 0) {
        puts("!!!CALCOLATA UNA DISTANZA NULLA!!!");
        exit(-1);
    }
    return sqrt(ris);
}

void copiaVettore(int dim, double input[dim], double output[dim]) {
    int i;
    for (i = 0; i < dim; i++) {
        output[i] = input[i];
    }
}

///////////////////////////////////////////////////////////////////////////////

long random_at_most(long max) {
    unsigned long

    num_bins = (unsigned long) max + 1,
            num_rand = (unsigned long) RAND_MAX + 1,
            bin_size = num_rand / num_bins,
            defect = num_rand % bin_size;

    long x;

    while (num_rand - defect <= (unsigned long) (x = random()));


    return x / bin_size;
}

double fRand(double fMin, double fMax) {//random double in un range
    double f = (double) lrand48() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

///////////////////////////////////////////////////////////////////////////////

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



    /*if (min_sep == 0) {
        puts("calcolaXB: DISTANZA NULLA MIN SEP");
        exit(-1);
    }*/
    if (sigma <= 0) {
        puts("calcolaXB: SIGMA NULLO");
        exit(-1);
    }

    /*double ris = (sigma) / (n * min_sep);

    if (ris <= 0) {
        puts("calcolato XB nullo");
        exit(-1);
    }*/


    return sigma/(c*min_sep);
}

void init(int n, int c, int d) {
    puts("###start init###");
    int row;
    //###INIT POPOLAZIONE 0
    //init V e calcolo U

    for (pop_index = 0; pop_index < num_pop_iniziale; pop_index++) {
        //alloc struttura
        POP_NEW[pop_index] = malloc(sizeof (el_pop));

        //alloc V_p
        POP_NEW[pop_index] -> V_p = malloc(c * sizeof (double*));
        for (row = 0; row < c; row++) {
            POP_NEW[pop_index] -> V_p[row] = malloc(d * sizeof (double));
        }

        //init V_p
        for (i = 0; i < c; i++) {
            for (j = 0; j < d; j++) {
                POP_NEW[pop_index] -> V_p[i][j] = X[random_at_most(n - 1)][random_at_most(d - 1)] * drand48();
            }
        }
        ////////////SINGOLO PASSO FCM

        //alloc U_p
        POP_NEW[pop_index] -> U_p = malloc(c * sizeof (double*));
        for (row = 0; row < c; row++) {
            POP_NEW[pop_index] -> U_p[row] = malloc(n * sizeof (double));
        }

        //init U
        for (i = 0; i < c; i++) {
            for (j = 0; j < n; j++) {
                double denom = 0.0;
                double dist_x_j__v_i = calcDistanza(X[j], POP_NEW[pop_index] -> V_p[i]);

                int k;
                for (k = 0; k < c; k++) {
                    double dist_xj_vk = calcDistanza(X[j], POP_NEW[pop_index] -> V_p[k]);
                    denom += pow((dist_x_j__v_i / dist_xj_vk), esponente_U);
                }
                POP_NEW[pop_index] -> U_p[i][j] = 1.0 / denom;
                if (POP_NEW[pop_index] -> U_p[i][j] <= 0) {
                    printf("ERR: INIT POP:inizializzato U di un membro con negativo::%lf", POP_NEW[pop_index] -> U_p[i][j]);
                    exit(-1);
                }
            }
        }
        //RICALCOLO POSIZIONE CENTROIDI
        for (i = 0; i < c; i++) {
            double denom = 0.0;
            for (j = 0; j < n; j++)//sommatoria denom (fatta una sola volta a centr.)
                denom += pow(POP_NEW[pop_index] ->U_p[i][j], m);
            for (k = 0; k < d; k++) {
                double num = 0.0;
                for (j = 0; j < n; j++) {//SOMMATORIA numeratore
                    num += X[j][k] * pow(POP_NEW[pop_index] ->U_p[i][j], m);
                }
                POP_NEW[pop_index] ->V_p[i][k] = num / denom;
            }
        }
        //////END FCM

        //calcolo fitness
        double xb = calcolaXB(POP_NEW[pop_index]->V_p, POP_NEW[pop_index]->U_p, 0);
        POP_NEW[pop_index]->indice_xb = xb;
    }
    //###END INIT POPOLAZIONE 0
    printf("\n########## FINE INIT #############\n");
}

void lavora(int n, int c, int d) {
    int numero_generazioni_utente = numero_generazioni;
    int row;
    do {//NUOVA GENERAZIONE
        //SCAMBIO VETTORI POPOLAZIONE
        int i_CPop;
        for (i_CPop = 0; i_CPop < num_pop_iniziale; i_CPop++) {
            if (POP_NEW[i_CPop] != POP_NOW[i_CPop]) {
                if (POP_NOW[i_CPop] != 0) {
                    for (row = 0; row < c; row++) {
                        free(POP_NOW[i_CPop] -> V_p[row]);
                        free(POP_NOW[i_CPop] -> U_p[row]);
                    }
                    free(POP_NOW[i_CPop]->V_p);
                    free(POP_NOW[i_CPop]->U_p);
                    free(POP_NOW[i_CPop]);
                }
                POP_NOW[i_CPop] = POP_NEW[i_CPop];
            }
        }

        /////////////////DE////////////////////
        for (i_CPop = 0; i_CPop < num_pop_iniziale; i_CPop++) {//PER OGNI COMPONENTE DELLA POP
            //SCELTA CANDIDATI
            //tre vettori devono essere scelti a caso nella popolazione
            //diversi dal target (indice i) e mutualmente
            int indice_1, indice_2, indice_3;
            do {
                indice_1 = random_at_most(((long) num_pop_iniziale) - 1);
            } while (indice_1 == i_CPop);

            do {
                indice_2 = random_at_most(((long) num_pop_iniziale) - 1);
            } while (indice_2 == i_CPop || indice_2 == indice_1);

            do {
                indice_3 = random_at_most(((long) num_pop_iniziale) - 1);
            } while (indice_3 == i_CPop || indice_3 == indice_1 || indice_3 == indice_2);

            if (indice_1 == i_CPop || indice_2 == i_CPop || indice_3 == i_CPop) {
                puts("INDICE NON VALIDO!!!!");
                exit(-1);
            }

            //l'elemento mutante
            el_pop *mutant = malloc(sizeof (el_pop));
            //alloc V_p del mutante
            mutant -> V_p = malloc(c * sizeof (double*));
            for (row = 0; row < c; row++) {
                mutant -> V_p[row] = malloc(d * sizeof (double));
            }
            //alloc U_p mutante
            mutant -> U_p = malloc(c * sizeof (double*));
            for (row = 0; row < c; row++) {
                mutant -> U_p[row] = malloc(n * sizeof (double));
            }

            //MUTATION
            int i1, j1;
            //scambio di geni (tipo 1)
            for (i1 = 0; i1 < c; i1++) {
                double f = fRand(0.1, 1.0); //differential weight
                for (j1 = 0; j1 < d; j1++) {
                    mutant->V_p[i1][j1] = POP_NOW[indice_3]->V_p[i1][j1] + f * (POP_NOW[indice_1]->V_p[i1][j1] - POP_NOW[indice_2]->V_p[i1][j1]);
                }
            }

            //CROSSOVER CON IL VETTORE ATTUALE (TIPO 1)
            for (i1 = 0; i1 < c; i1++) {
                double prob_crossover = fRand(0.0, 1.0);
                for (j1 = 0; j1 < d; j1++) {
                    if (prob_crossover < CR) {
                        //prendo il gene del vettore attuale
                        mutant->V_p[i1][j1] = POP_NOW[i_CPop]->V_p[i1][j1];
                    }
                }
            }

            //CROSSOVER CON IL VETTORE ATTUALE (TIPO 2)
            /*for (i1 = 0; i1 < c; i1++) {
                for (j1 = 0; j1 < d; j1++) {
                    double prob_crossover = fRand(0.0, 1.0);
                    if (prob_crossover < CR) {
                        //prendo il gene del vettore attuale
                        trial->V_p[i1][j1] = POP_NOW[i_CPop]->V_p[i1][j1];
                    }
                }
            }*/

            //CROSSOVER CON IL VETTORE ATTUALE (TIPO 3)
            /*for (i1 = 0; i1 < c; i1++) {
                    double prob_crossover = fRand(0.0, 1.0);
                    if (prob_crossover < CR) {
                        //prendo il gene del vettore attuale
                        copiaVettore(d,POP_NOW[i_CPop]->V_p[i1],trial->V_p[i1]);
                    }
            }*/

            ///CALCOLO FITNESS DEL MUTANTE
            //calcolo U mutante
            for (i = 0; i < c; i++) {
                for (j = 0; j < n; j++) {
                    double denom = 0.0;
                    double dist_x_j__v_i = calcDistanza(X[j], mutant -> V_p[i]);
                    if (dist_x_j__v_i == 0) {
                        puts("calcolo U mutante, distanza nulla");
                        exit(-1);
                    }
                    //printf("\ndist_x_j__v_i: %lf,%d,%d\n",dist_x_j__v_i,c,n);
                    int k;
                    for (k = 0; k < c; k++) {
                        double dist_xj_vk = calcDistanza(X[j], mutant -> V_p[k]);
                        if (dist_xj_vk == 0) {
                            puts("calcolo U mutante, distanza nulla");
                            exit(-1);
                        }
                        denom += pow((dist_x_j__v_i / dist_xj_vk), esponente_U);
                        if (denom == 0) {
                            puts("calcolo U mutante, denom nullo");
                            exit(-1);
                        }

                    }
                    mutant -> U_p[i][j] = 1.0 / denom;
                }
            }

            //calcolo XB mutante            
            mutant->indice_xb = calcolaXB(mutant->V_p, mutant->U_p, 1);
            /*if (trial -> indice_xb <= guardia_xb_1) {
                puts("*ERROR: indice_xb trial inferiore a guardia, scartato");
                trial->indice_xb = DBL_MAX;
            }*/
            //////END CALCOLO FITNESS MUTANTE

            //SELECTION
            if (mutant->indice_xb < POP_NOW[i_CPop]->indice_xb) {
                xb_selezionato = mutant->indice_xb;
                POP_NEW[i_CPop] = mutant;
            } else {
                for (row = 0; row < c; row++) {
                    free(mutant -> V_p[row]);
                    free(mutant -> U_p[row]);
                }
                free(mutant->V_p);
                free(mutant->U_p);
                free(mutant);
                POP_NEW[i_CPop] = POP_NOW[i_CPop];
                xb_selezionato = POP_NOW[i_CPop]->indice_xb;
            }
        }//END DE
        numero_generazioni--;
    } while (numero_generazioni > 0);

    puts("calcolo xb finale");
    //computazione XB della popolazione finale
    double best_xb = DBL_MAX;
    int indice_best;
    for (pop_index = 0; pop_index < num_pop_iniziale; pop_index++) {
        POP_NEW[pop_index]->indice_xb = calcolaXB(POP_NEW[pop_index]->V_p, POP_NEW[pop_index]->U_p, 0);
        if (POP_NEW[pop_index]->indice_xb < best_xb) {
            best_xb = POP_NEW[pop_index]->indice_xb;
            indice_best = pop_index;
        }
    }


    puts("***********************************************");
    /*if (xb_selezionato < guardia_xb_2)
        puts("***SUPERATA GUARDIA XB_2***");*/
    int numero_gen_fatte = numero_generazioni_utente - numero_generazioni;
    printf("\nnumero di elementi delle popolazioni:%d\n", num_pop_iniziale);
    fprintf(out_LOG_RIS, "\nnumero di elementi delle popolazioni:%d\n", num_pop_iniziale);
    printf("\nnumero di generazioni max:%d\n", numero_generazioni_utente);
    fprintf(out_LOG_RIS, "\nnumero di generazioni max:%d\n", numero_generazioni_utente);
    printf("\nnumero generazioni fatte:%d\n", numero_gen_fatte);
    fprintf(out_LOG_RIS, "\nnumero generazioni fatte:%d\n", numero_gen_fatte);
    printf("\nmiglior XB:%lf\n", best_xb);
    fprintf(out_LOG_RIS, "\nmiglior XB:%lf\n", best_xb);


    puts("matrice V:");
    stampaMatrice(c, d, POP_NEW[indice_best]->V_p);
    puts("");
    stampaMatriceSuFile(c, d, POP_NEW[indice_best]->V_p, out_V);
    stampaMatriceSuFile(c, n, POP_NEW[indice_best]->U_p, out_U);
    puts("***********************************************");
}

void plot() {
    if ((d == 2 || d == 3) && attivaGnuPlot) {
        char *commandsForGnuplot[] = {"set key off", "set term wxt 1", "set title \"matrice X\"", "", "set term wxt 2", "set key off", "set title \"matrice V\"", ""};
        if (d == 2) {
            commandsForGnuplot[3] = "plot 'x.dat'";
            commandsForGnuplot[7] = "plot 'v.dat'";
        } else if (d == 3) {
            commandsForGnuplot[3] = "splot 'x.dat'";
            commandsForGnuplot[7] = "splot 'v.dat'";
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
    //stream file
    out_X = fopen("x.dat", "r");
    out_V = fopen("v.dat", "w");
    out_U = fopen("u.dat", "w");
    out_LOG = fopen("log", "w");
    out_LOG_RIS = fopen("log_ris", "a");

    esponente_U = 2.0 / (m - 1.0);

    //letture da utente
    puts("numero di punti in input");
    scanf("%d", &n);
    puts("numero di dimensioni");
    scanf("%d", &d);
    puts("numero di centroidi");
    scanf("%d", &c);
    puts("numero di generazioni:");
    scanf("%d", &numero_generazioni);

    num_pop_iniziale = num_pop;


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

    fputs("\n######\n", out_LOG_RIS);
    fprintf(out_LOG_RIS, "\nnumero dimensioni:%d\n", d);
    fprintf(out_LOG_RIS, "\nnumero centroidi:%d\n", c);

    init(n, c, d);
    lavora(n, c, d);
    plot();

    return (0);
}

