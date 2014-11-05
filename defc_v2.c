#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#include <string.h>

FILE *out_V, *out_X, *out_U, *out_LOG, *out_LOG_RIS, *out_csv; //puntatori a file di output

double **X;
double m = 2.0; //fuzzification
double CR; //crossover rate
//bound del differential weight
double dw_lowerbound;
double dw_upperbound;
int dw_adattivo;
int tipo_dw;
int tipo_dataset;

#define num_pop 20

typedef struct el_pop {//individuo della popolazione
    double **V_p;
    double **U_p;
    double fitness;
} el_pop;

el_pop *POP_NEW[num_pop];
el_pop *POP_NOW[num_pop];

//attiva e disattiva GnuPlot
int attivaGnuPlot = 0;



int num_pop_iniziale, numero_generazioni, i, j, k, pop_index, numero_generazione_attuale;
double esponente_U;
double xb_selezionato;
double best_xb;
int n, c, d;
int range_init_min, range_init_max;

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

double dbl_rnd_inRange(double fMin, double fMax) {//random double in un range
    double f = (double) lrand48() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

///////////////////////////////////////////////////////////////////////////////

double calcolaFitness(double **V, double **U, int debug) {
    //CALCOLO SIGMA
    double sigma = 0.0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < c; j++) {
            sigma += pow(U[j][i], m) * pow(calcDistanza(V[j], X[i]), 2.0);
        }
    }
    return sigma;
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
        srand(rand());
        //init V_p
        for (i = 0; i < c; i++) {
            for (j = 0; j < d; j++) {
                //POP_NEW[pop_index] -> V_p[i][j] = X[random_at_most(n-1)][random_at_most(d-1)]*drand48();//;
                //POP_NEW[pop_index] -> V_p[i][j] = drand48();
                //POP_NEW[pop_index] -> V_p[i][j] = dbl_rnd_inRange(0,range_init_max);
                POP_NEW[pop_index] -> V_p[i][j] = random_at_most(range_init_max);
            }
        }

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
                if (POP_NEW[pop_index] -> U_p[i][j] < 0) {
                    printf("ERR: INIT POP:inizializzato U di un membro con negativo::");
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
        double xb = calcolaFitness(POP_NEW[pop_index]->V_p, POP_NEW[pop_index]->U_p, 0);
        POP_NEW[pop_index]->fitness = xb;
    }
    //###END INIT POPOLAZIONE 0
    printf("\n########## FINE INIT #############\n");
}

void lavora(int n, int c, int d) {
    int row;
    numero_generazione_attuale = 0;
    do {//NUOVA GENERAZIONE
        //SCAMBIO VETTORI POPOLAZIONE
        numero_generazione_attuale++;
        printf(".");
        fflush(stdout);
        int i_target;
        for (i_target = 0; i_target < num_pop_iniziale; i_target++) {
            if (POP_NEW[i_target] != POP_NOW[i_target]) {
                if (POP_NOW[i_target] != 0) {
                    for (row = 0; row < c; row++) {
                        free(POP_NOW[i_target] -> V_p[row]);
                        free(POP_NOW[i_target] -> U_p[row]);
                    }
                    free(POP_NOW[i_target]->V_p);
                    free(POP_NOW[i_target]->U_p);
                    free(POP_NOW[i_target]);
                }
                POP_NOW[i_target] = POP_NEW[i_target];
            }
        }

        //SPERIMENTALE: cambio dinamico di dw_upperbound
        //if (dw_adattivo && numero_generazione_attuale > 0 && dw_upperbound > dw_lowerbound)
        //    dw_upperbound = dw_upperbound / 1.1;

        int indice_base;
        //indice_base = random_at_most(((long) num_pop_iniziale) - 1);
        
        /////////////////DE////////////////////
        for (i_target = 0; i_target < num_pop_iniziale; i_target++) {//PER OGNI COMPONENTE DELLA POP
            /*if(indice_base == i_target){
                printf("!");
                continue;
            }*/
                
            //SCELTA CANDIDATI
            //tre vettori devono essere scelti a caso nella popolazione
            //diversi dal target (indice i) e mutualmente
            
            int indice_1, indice_2;
            do {
                indice_1 = random_at_most(((long) num_pop_iniziale) - 1);
            } while (indice_1 == i_target); // || indice_1 == indice_base

            do {
                indice_2 = random_at_most(((long) num_pop_iniziale) - 1);
            } while (indice_2 == i_target || indice_2 == indice_1); // || indice_2 == indice_base

            do {
                indice_base = random_at_most(((long) num_pop_iniziale) - 1);
            } while (indice_base == i_target || indice_base == indice_1 || indice_base == indice_2);

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
            for (i1 = 0; i1 < c; i1++) {
                double f;
                if (tipo_dw == 2)//dithering
                    f = dbl_rnd_inRange(dw_lowerbound, dw_upperbound);
                else if (tipo_dw == 3) //dither + jitter
                    f = dbl_rnd_inRange(dw_lowerbound, dw_upperbound) + 0.001 * (dbl_rnd_inRange(0, 1) - 0.5);
                else
                    f = dw_upperbound;

                for (j1 = 0; j1 < d; j1++) {
                    mutant->V_p[i1][j1] = POP_NOW[indice_base]->V_p[i1][j1] + f * (POP_NOW[indice_1]->V_p[i1][j1] - POP_NOW[indice_2]->V_p[i1][j1]);
                }
            }
            //CROSSOVER CON IL VETTORE TARGET (TIPO 1)
            for (i1 = 0; i1 < c; i1++) {
                double prob_crossover = dbl_rnd_inRange(0.0, 1.0);
                for (j1 = 0; j1 < d; j1++) {
                    if (prob_crossover < CR) {
                        //prendo il cromosoma del target
                        mutant->V_p[i1][j1] = POP_NOW[i_target]->V_p[i1][j1];
                    }
                }
            }

            //CROSSOVER CON IL VETTORE ATTUALE (TIPO 2) (BAD)
            /*for (i1 = 0; i1 < c; i1++) {
                for (j1 = 0; j1 < d; j1++) {
                    double prob_crossover = fRand(0.0, 1.0);
                    if (prob_crossover < CR) {
                        mutant->V_p[i1][j1] = POP_NOW[i_CPop]->V_p[i1][j1];
                    }
                }
            }*/

            //CROSSOVER CON IL VETTORE ATTUALE (TIPO 3)
            /*for (i1 = 0; i1 < c; i1++) {
                    double prob_crossover = fRand(0.0, 1.0);
                    if (prob_crossover < CR) {
                        //prendo tutto il vettore del compagno
                        copiaVettore(d,POP_NOW[i_CPop]->V_p[i1],mutant->V_p[i1]);
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
            mutant->fitness = calcolaFitness(mutant->V_p, mutant->U_p, 1);

            //SELECTION
            //TIPO 1
            /*int indice_sostituito = 0;
            if (mutant->fitness < POP_NOW[i_CPop]->fitness) {
                indice_sostituito = i_CPop;
                xb_selezionato = mutant->fitness;
                POP_NEW[i_CPop] = mutant;
            } else if (mutant->fitness < POP_NOW[indice_1]->fitness) {
                indice_sostituito = indice_1;
                xb_selezionato = mutant->fitness;
                POP_NEW[indice_1] = mutant;
            } else if (mutant->fitness < POP_NOW[indice_2]->fitness) {
                indice_sostituito = indice_2;
                xb_selezionato = mutant->fitness;
                POP_NEW[indice_2] = mutant;
            } else if (mutant->fitness < POP_NOW[indice_3]->fitness) {
                indice_sostituito = indice_3;
                xb_selezionato = mutant->fitness;
                POP_NEW[indice_3] = mutant;
            } else {
                for (row = 0; row < c; row++) {
                    free(mutant -> V_p[row]);
                    free(mutant -> U_p[row]);
                }
                free(mutant->V_p);
                free(mutant->U_p);
                free(mutant);
                POP_NEW[i_CPop] = POP_NOW[i_CPop];
                xb_selezionato = POP_NOW[i_CPop]->fitness;
            }*/
            //TIPO 2 (STANDARD)
            if (mutant->fitness < POP_NOW[i_target]->fitness) {
                xb_selezionato = mutant->fitness;
                POP_NEW[i_target] = mutant;
                
            } else {
                for (row = 0; row < c; row++) {
                    free(mutant -> V_p[row]);
                    free(mutant -> U_p[row]);
                }
                free(mutant->V_p);
                free(mutant->U_p);
                free(mutant);
                POP_NEW[i_target] = POP_NOW[i_target];
                xb_selezionato = POP_NOW[i_target]->fitness;
            }
        }//END DE
        numero_generazioni--;
    } while (numero_generazioni > 0);


    //computazione fitness della popolazione finale
    double best_fitness = DBL_MAX;
    int indice_best;
    for (pop_index = 0; pop_index < num_pop_iniziale; pop_index++) {
        POP_NEW[pop_index]->fitness = calcolaFitness(POP_NEW[pop_index]->V_p, POP_NEW[pop_index]->U_p, 0);
        if (POP_NEW[pop_index]->fitness < best_fitness) {
            best_fitness = POP_NEW[pop_index]->fitness;
            indice_best = pop_index;
        }
    }
    //calcolo xb del migliore
    best_xb = calcolaXB(POP_NEW[indice_best]->V_p, POP_NEW[indice_best]->U_p, 0);

    puts("\n\n***********************************************");
    printf("miglior XB:%lf\n", best_xb);
    fprintf(out_LOG_RIS, "\nmiglior XB:%lf\n\n", best_xb);
    stampaMatriceSuFile(c, d, POP_NEW[indice_best]->V_p, out_LOG_RIS);
    puts("matrice V:");
    stampaMatrice(c, d, POP_NEW[indice_best]->V_p);
    stampaMatriceSuFile(c, d, POP_NEW[indice_best]->V_p, out_V);
    stampaMatriceSuFile(c, n, POP_NEW[indice_best]->U_p, out_U);
    puts("***********************************************");
}

void plot() {
    if ((d == 2 || d == 3) && attivaGnuPlot) {
        char *commandsForGnuplot[] = {"set key off", "set term x11 1", "set title \"matrice X\"", "", "set term x11 2", "set key off", "set title \"DEFC - matrice V\"", ""};
        if (d == 2) {
            commandsForGnuplot[3] = "plot 'x.dat' pointtype 3";
            commandsForGnuplot[7] = "plot 'v_defc.dat' pointtype 3";
        } else if (d == 3) {
            commandsForGnuplot[3] = "splot 'x.dat' pointtype 3";
            commandsForGnuplot[7] = "splot 'v_defc.dat' pointtype 3";
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
    out_V = fopen("v_defc.dat", "w");
    out_U = fopen("u_defc.dat", "w");
    out_LOG_RIS = fopen("log_ris", "a");
    out_csv = fopen("output.csv", "a");

    esponente_U = 2.0 / (m - 1.0);
    num_pop_iniziale = num_pop;


    //letture da utente
    printf("tipo dataset: ");
    scanf("%d",&tipo_dataset);
    printf("numero di punti in input: ");
    scanf("%d", &n);
    printf("numero di dimensioni: ");
    scanf("%d", &d);
    printf("numero di centroidi: ");
    scanf("%d", &c);
    printf("numero di generazioni: ");
    scanf("%d", &numero_generazioni);
    printf("crossover rate (reale tra 0 e 1): ");
    scanf("%lf", &CR);
    printf("differential weight upperbound: ");
    scanf("%lf", &dw_upperbound);
    printf("range init max: ");
    scanf("%d", &range_init_max);
    printf("tipo diff. weight: 1-costante 2-dithered 3-dithered/jittered: ");
    scanf("%d", &tipo_dw);

    int ngenIniziali = numero_generazioni;
    dw_lowerbound = 0.001;
    fputs("\n######\n", out_LOG_RIS); //nuova log entry


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

    init(n, c, d);
    lavora(n, c, d);
    plot();

    //aggiornamento log
    fprintf(out_LOG_RIS, "punti in input:%d\n", n);
    fprintf(out_LOG_RIS, "dimensioni:%d\n", d);
    fprintf(out_LOG_RIS, "centroidi:%d\n", c);
    fprintf(out_LOG_RIS, "popolazione:%d\n", num_pop_iniziale);
    fprintf(out_LOG_RIS, "numero generazioni:%d\n", numero_generazione_attuale);
    fprintf(out_LOG_RIS, "crossover rate:%lf\n", CR);
    fprintf(out_LOG_RIS, "dw lowerbound:%lf\n", dw_lowerbound);
    fprintf(out_LOG_RIS, "dw upperbound:%lf\n", dw_upperbound);

    //scrittura csv
    fprintf(out_csv,"%d,",tipo_dataset);
    fprintf(out_csv, "%d,", n);
    fprintf(out_csv, "%d,", c);
    fprintf(out_csv, "%d,", d);
    fprintf(out_csv, "%lf,", CR);
    fprintf(out_csv, "%d,", num_pop);
    fprintf(out_csv, "%d,", ngenIniziali);
    fprintf(out_csv, "%d,", tipo_dw);
    fprintf(out_csv, "%lf,", dw_lowerbound);
    fprintf(out_csv, "%lf,", dw_upperbound);
    fprintf(out_csv, "%d,", range_init_min);
    fprintf(out_csv, "%d,", range_init_max);
    fprintf(out_csv, "%lf\n", best_xb);

    return (0);
}

