#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#include <string.h>

FILE *in_V, *out_V, *out_X, *out_U, *out_LOG, *out_LOG_RIS, *out_csv, *debug_log, *dataset; //puntatori a file di output

double **X;
double m; //fuzzification factor
double CR; //crossover rate
int tipo_dataset;
int starting_age;
int ultimo_conteggio_successi;
int conteggio_successi_generazione_attuale;
int conteggio_reset;
int reset_threshold;
int bestFitIndex;
int bestXBIndex;
double best_xb;
double best_fit;
int numero_generazioni, i, j, k, pop_index, numero_generazione_attuale, numero_generazioni_iniziale;
double esponente_U;
int n, c, d; //numero di punti in input, numero di centroidi, numero di dimensioni
int abilita_partitioning; //riordina gli array V per far accoppiare solo i centroidi della stessa zona dello spazio
int abilita_invecchiamento; //invecchiamento e reinizializzazione
int usa_xb_per_fitness; //usa XB invece di sigma
int abilita_shuffle; //shuffle posizione centroidi
int abilita_reset; //reset generale popolazione
int attivaGnuPlot; //attiva e disattiva GnuPlot
int testLoadVIdeale; //test solo per dataset gauss 4
int random_init; //inizializza V popolazione iniziale completamente in modo casuale
int aggiungi_peso_sigma;
int usa_sumsep;
int attiva_partitioned_init;
int attiva_sigma_separate;
double peso_sigma;
double soglia_conteggio;
double soglia_peso_sigma;

//numero di elementi della popolazione - fare parametrico
#define num_pop 100 // 30, 50, 100

//struttura elemento popolazione

typedef struct el_pop {
    double **V_p;
    double **U_p;
    double fitness;
    double XB;
    int age;
    //jDE
    double f;
    double CR;
} el_pop;

//vettore nuova popolazione
el_pop *POP_NEW[num_pop];

//vettore popolazione attuale
el_pop *POP_NOW[num_pop];

//per confronto veloce delle fitness conosciute
double fitness_vector[num_pop];

//per confronti XB
double xb_vector[num_pop];

void stampaMatrice(int righe, int col, double **mat) {
    int i, j;
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

unsigned int rand_int_in_interval(unsigned int min, unsigned int max) {
    int r;
    const unsigned int range = 1 + max - min;
    const unsigned int buckets = RAND_MAX / range;
    const unsigned int limit = buckets * range;

    do {
        r = rand();
    } while (r >= limit);

    return min + (r / buckets);
}

int rand_int(int n) {
    int limit = RAND_MAX - RAND_MAX % n;
    int rnd;

    do {
        rnd = rand();
    } while (rnd >= limit);
    return rnd % n;
}

///////////////////////////////////////////////////////////////////////////////

double calcolaSigma(double **V, double **U) {
    double sigma = 0;
    //CALCOLO SIGMA
    if (!attiva_sigma_separate) {
        for (i = 0; i < c; i++) {

            V[i][d + 1] = 0; //reset contatore elementi

            for (j = 0; j < n; j++) {
                sigma += pow(U[i][j], m) * pow(calcDistanza(V[i], X[j]), 2.0);

                //conteggio elementi del cluster
                if (U[i][j] > soglia_conteggio)
                    V[i][d + 1]++;
            }



            if (aggiungi_peso_sigma && V[i][d + 1] > soglia_peso_sigma) {
                peso_sigma = 1 - (V[i][d + 1] / n); //sperimentale
                //aggiunta di un peso per via del numero di elementi vicini al centroide
                sigma = sigma * peso_sigma;
            }
        }
    } else {
        //calcolo le sigma separatamente
        for (i = 0; i < c; i++) {

            V[i][d + 1] = 0; //reset conteggio elementi cluster

            for (j = 0; j < n; j++) {
                V[i][d] += pow(U[i][j], m) * pow(calcDistanza(V[i], X[j]), 2.0);

                //conteggio elementi del cluster
                if (U[i][j] > soglia_conteggio)
                    V[i][d + 1]++;
            }

            if (aggiungi_peso_sigma && V[i][d + 1] > soglia_peso_sigma) {
                peso_sigma = 1 - (V[i][d + 1] / n); //1 - (V[i][d + 1] / n); //sperimentale
                //aggiunta di un peso per via del numero di elementi vicini al centroide
                V[i][d] = V[i][d] * peso_sigma;
            }
        }

        //somma delle sigma separate
        for (i = 0; i < c; i++) {
            sigma += V[i][d];
        }
    }
    return sigma;
}

double calcolaMinSep(double **V) {
    //CALCOLO MIN_SEP
    double den;
    if (!usa_sumsep) {
        double min_sep = DBL_MAX;
        double dist_tmp = 0;
        j = 0;
        for (i = 0; i < c; i++) {
            if (j == i)
                j++;
            while (j < c) {
                if (j == i)
                    j++;
                if (j < c) {
                    dist_tmp = (pow(calcDistanza(V[i], V[j]), 2.0));
                    if (dist_tmp < min_sep)
                        min_sep = dist_tmp;
                    j++;
                }
            }
            j = 0;
        }

        den = min_sep;
    } else {
        double distsum = 0;
        //CALCOLO SUM_SEP
        //invece di calcolare la minima separazione si prende la somma delle distanze fra i centroidi
        for (i = 0; i < c; i++) {
            if (j == i)
                j++;
            while (j < c) {
                if (j == i)
                    j++;
                if (j < c) {
                    distsum += (pow(calcDistanza(V[i], V[j]), 2.0));
                    j++;
                }
            }
            j = 0;
        }
        den = distsum;
    }

    return den;
}

double calcolaXB(double **V, double **U) {
    //double sigma;
    double sigma = calcolaSigma(V, U);
    /*for (i = 0; i < n; i++) {
        for (j = 0; j < c; j++) {
            sigma += pow(U[j][i], m) * pow(calcDistanza(V[j], X[i]), 2.0);
        }
    }
     */
    double den = calcolaMinSep(V);

    return sigma / (n * den);
}

double calcolaFitness(double **V, double **U) {
    if (usa_xb_per_fitness)
        return calcolaXB(V, U);
    else
        return calcolaSigma(V, U);
}

//usata per partizionare in sottopopolazioni

void sortMatrice(double **V) {
    if (abilita_partitioning) {
        double temp;
        int scambiare;
        int alto, riga, colonna;

        for (alto = c - 1; alto > 0; alto--) {
            for (riga = 0; riga < alto; riga++) {
                scambiare = 0;
                for (colonna = 0; colonna < d; colonna++) {
                    if (colonna == 0) {//ordinamento rispetto solo alla prima dimensione
                        if (V[riga][colonna] > V[riga + 1][colonna]) {
                            scambiare = 1;
                            temp = V[riga][colonna];
                            V[riga][colonna] = V[riga + 1][colonna];
                            V[riga + 1][colonna] = temp;
                        }
                    } else {
                        if (scambiare == 1) {//spostamento degli altri elementi del vettore da spostare
                            temp = V[riga][colonna];
                            V[riga][colonna] = V[riga + 1][colonna];
                            V[riga + 1][colonna] = temp;
                        }
                    }
                }
            }
        }
    }
}

//usato per lo shake della matrice V

void shuffleMatrice(double **V) {
    if (abilita_shuffle) {
        int riga, colonna, destinazione;
        double temp;

        for (riga = c - 1; riga > 0; riga--) {//riga centroide
            for (colonna = 0; colonna < d; colonna++) {//colonna dimensione
                if (colonna == 0) {
                    destinazione = rand_int(riga + 1);
                    temp = V[riga][colonna];
                    V[riga][colonna] = V[destinazione][colonna];
                    V[destinazione][colonna] = temp;
                } else {
                    temp = V[riga][colonna];
                    V[riga][colonna] = V[destinazione][colonna];
                    V[destinazione][colonna] = temp;
                }
            }
        }
    }
}

void aggiornaBestFitIndex() {
    double bestFit = DBL_MAX;
    for (i = 0; i < num_pop; i++) {
        if (fitness_vector[i] < bestFit) {
            bestFit = fitness_vector[i];
            bestFitIndex = i;
        }
    }
}

void aggiornaBestXBIndex() {
    double temp = DBL_MAX;
    for (i = 0; i < num_pop; i++) {
        if (xb_vector[i] < temp) {
            temp = xb_vector[i];
            bestXBIndex = i;
        }
    }
}

//INIZIALIZZAZIONE GENERALE DI TUTTA LA POPOLAZIONE INIZIALE

void init(int n, int c, int d) {
    int row;
    //###INIT POPOLAZIONE 0
    //init V e calcolo U

    for (pop_index = 0; pop_index < num_pop; pop_index++) {
        //alloc struttura
        POP_NEW[pop_index] = malloc(sizeof (el_pop));

        //alloc V_p
        //in posizione V_p[d] va messo numero di elementi del cluster
        POP_NEW[pop_index] -> V_p = malloc((c * sizeof (double*)));
        for (row = 0; row < c; row++) {
            POP_NEW[pop_index] -> V_p[row] = malloc((d * sizeof (double)) + (2 * sizeof (double)));
        }

        int rigaX;
        if (testLoadVIdeale && pop_index == 0) {//carica una V da file per testare la funzione obiettivo
            //SOLO PER TEST, leggo V da file esterno
            for (i = 0; i < c; i++) {
                for (j = 0; j < d; j++) {
                    if (!fscanf(in_V, "%lf", &POP_NEW[pop_index] -> V_p[i][j]))
                        break;
                }
                //servono per peso sigma e sigma separate
                POP_NEW[pop_index]->V_p[i][d] = 0;
                POP_NEW[pop_index]->V_p[i][d + 1] = 0;
            }
        } else {
            if (attiva_partitioned_init) {
                int conteggio_bin = 0;
                int grandezza_bin = floor(n / c);
                //inizializzo c-esimo centroide nella n/c - esima partizione dell'input
                for (i = 0; i < c; i++) {
                    if (i == c - 1)
                        rigaX = rand_int_in_interval((grandezza_bin * conteggio_bin) - grandezza_bin, (grandezza_bin * conteggio_bin));
                    else
                        rigaX = rand_int_in_interval(grandezza_bin*conteggio_bin, (grandezza_bin * conteggio_bin) + grandezza_bin);
                    for (j = 0; j < d; j++) {
                        POP_NEW[pop_index] -> V_p[i][j] = X[rigaX][j] + drand48();
                    }
                    conteggio_bin++;
                    //servono per peso sigma e sigma separate
                    POP_NEW[pop_index]->V_p[i][d] = 0;
                    POP_NEW[pop_index]->V_p[i][d + 1] = 0;
                }
            } else {
                //init V_p STANDARD
                for (i = 0; i < c; i++) {
                    if (!random_init)
                        rigaX = random_at_most(n - 1);
                    for (j = 0; j < d; j++) {
                        if (random_init)
                            POP_NEW[pop_index] -> V_p[i][j] = drand48();
                        else
                            POP_NEW[pop_index] -> V_p[i][j] = X[rigaX][j] + drand48();
                    }
                    //servono per peso sigma e sigma separate
                    POP_NEW[pop_index]->V_p[i][d] = 0;
                    POP_NEW[pop_index]->V_p[i][d + 1] = 0;
                }
            }
        }

        //riordino matrice
        sortMatrice(POP_NEW[pop_index]->V_p);

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

        //calcolo fitness
        double fit = calcolaFitness(POP_NEW[pop_index]->V_p, POP_NEW[pop_index]->U_p);
        POP_NEW[pop_index]->fitness = fit;
        fitness_vector[pop_index] = fit;
        //calcolo XB
        double xb = calcolaXB(POP_NEW[pop_index]->V_p, POP_NEW[pop_index]->U_p);
        POP_NEW[pop_index]->XB = xb;
        xb_vector[pop_index] = xb;


        //testLoadVIdeale
        if (testLoadVIdeale && pop_index == 0) {
            printf("fitness ideale:\t%lf\n", POP_NEW[pop_index]->fitness);
            printf("XB ideale\t%lf\n", POP_NEW[pop_index]->XB);
        }


        //impostazione età
        POP_NEW[pop_index] -> age = starting_age;

        //impostazione f e CR
        POP_NEW[pop_index]->f = dbl_rnd_inRange(0.1, 1.0);
        POP_NEW[pop_index]->CR = dbl_rnd_inRange(0.0, 1.0);
    }
    //###END INIT POPOLAZIONE 0
    printf("#end init#\n");

}

void lavora(int n, int c, int d) {
    int row;
    numero_generazione_attuale = 0;
    do {//NUOVA GENERAZIONE
        //SCAMBIO VETTORI POPOLAZIONE
        conteggio_successi_generazione_attuale = 0;
        numero_generazione_attuale++;
        int i_target;
        for (i = 0; i < num_pop; i++) {
            if (POP_NEW[i] != POP_NOW[i]) {//l'individuo era stato sostituito da un trial
                if (POP_NOW[i] != 0) {//free del vecchio individuo
                    for (row = 0; row < c; row++) {
                        free(POP_NOW[i] -> V_p[row]);
                        free(POP_NOW[i] -> U_p[row]);
                    }
                    free(POP_NOW[i]->V_p);
                    free(POP_NOW[i]->U_p);
                    free(POP_NOW[i]);
                }
                POP_NOW[i] = POP_NEW[i]; //sostituzione puntatore
            }
        }

        aggiornaBestXBIndex();
        aggiornaBestFitIndex();

        //DEBUG/////////////////////////
        if (numero_generazione_attuale % 50 == 0) {
            printf("\n#####DEBUG#####\n");
            double best_xb = POP_NOW[bestXBIndex]->XB;
            double best_fit = POP_NOW[bestFitIndex]->fitness;
            double best_xb_and_fit = POP_NOW[bestFitIndex]->XB;
            printf("Generazione:%d      best_fit_index:%d       best_xb_index:%d\n", numero_generazione_attuale, bestFitIndex, bestXBIndex);
            printf("best_XB_assoluto:%lf    best_fit:%lf    xb_of_best_fit:%lf\n", best_xb, best_fit, best_xb_and_fit);
            fprintf(debug_log, "%d\n", numero_generazione_attuale);
            stampaMatriceSuFile(c, d, POP_NOW[bestXBIndex]->V_p, debug_log);
            fprintf(debug_log, "\n");
            stampaMatriceSuFile(c, d, POP_NOW[bestFitIndex]->V_p, debug_log);
            fprintf(debug_log, "\n");
            fprintf(debug_log, "\n");
        }
        /////////////////////////////

        /////////////////DE////////////////////
        for (i_target = 0; i_target < num_pop; i_target++) {//PER OGNI COMPONENTE DELLA POP
            //SCELTA CANDIDATI
            //tre vettori devono essere scelti a caso nella popolazione
            //diversi dal target (indice i) e mutualmente

            int indice_1, indice_2, indice_base;
            do {
                indice_1 = random_at_most(((long) num_pop) - 1);
            } while (indice_1 == i_target); // || indice_1 == indice_base

            do {
                indice_2 = random_at_most(((long) num_pop) - 1);
            } while (indice_2 == i_target || indice_2 == indice_1); // || indice_2 == indice_base

            do {
                indice_base = random_at_most(((long) num_pop) - 1);
            } while (indice_base == i_target || indice_base == indice_1 || indice_base == indice_2);

            //l'elemento mutante
            el_pop * mutant = malloc(sizeof (el_pop));
            //alloc V_p del mutante
            mutant -> V_p = malloc(c * sizeof (double*));
            for (row = 0; row < c; row++) {
                mutant -> V_p[row] = malloc((d * sizeof (double)) + (2 * sizeof (double))); //usato per conservare la dimensione del cluster
            }

            //alloc U_p mutante
            mutant -> U_p = malloc(c * sizeof (double*));
            for (row = 0; row < c; row++) {
                mutant -> U_p[row] = malloc(n * sizeof (double));
            }


            int i1, j1;

            //jDE - generazione di nuovi parametri f e CR che possibilmente diventeranno i nuovi del target o del trial
            //a seconda di quale dei due sopravvive
            double f, CR;

            double prob_newCR, prob_newF;
            prob_newCR = dbl_rnd_inRange(0.0, 1.0);
            if (prob_newCR < 0.1)
                f = POP_NOW[i_target]->f;
            else
                f = dbl_rnd_inRange(0.1, 1.0);

            prob_newF = dbl_rnd_inRange(0.0, 1.0);
            if (prob_newF < 0.1)
                CR = POP_NOW[i_target]->CR;
            else
                CR = dbl_rnd_inRange(0.0, 1.0);
            //////


            //MUTATION (TIPO 1 rand-rand-rand) singolo punto
            for (i1 = 0; i1 < c; i1++) {//centroide
                for (j1 = 0; j1 < d; j1++) {//dimensione
                    mutant->V_p[i1][j1] = POP_NOW[indice_base]->V_p[i1][j1] + f * (POP_NOW[indice_1]->V_p[i1][j1] - POP_NOW[indice_2]->V_p[i1][j1]);
                }
                //servono per peso sigma e sigma separate
                mutant->V_p[i1][d] = 0;
                mutant->V_p[i1][d + 1] = 0;
            }

            //CROSSOVER CON IL VETTORE TARGET UNIFORME
            double prob_crossover;
            for (i1 = 0; i1 < c; i1++) {
                prob_crossover = dbl_rnd_inRange(0.0, 1.0);
                for (j1 = 0; j1 < d; j1++) {
                    if (prob_crossover < CR) {
                        mutant->V_p[i1][j1] = POP_NOW[i_target]->V_p[i1][j1];
                    }
                }
            }

            //SORT MATRICE MUTANTE (ORA DETTO TRIAL)
            sortMatrice(mutant->V_p);

            //calcolo U trial
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

            //calcolo fitness trial   
            mutant->fitness = calcolaFitness(mutant->V_p, mutant->U_p);

            //calcolo XB trial
            mutant->XB = calcolaXB(mutant->V_p, mutant->U_p);

            //impostazione timer età trial
            mutant -> age = starting_age;

            //SELECTION
            //selezione può essere fatta su fitness o su XB
            if (mutant->fitness < POP_NOW[i_target]->fitness) {// && mutant->XB < POP_NOW[i_target]->XB) {// || (mutant->XB < POP_NOW[i_target]->XB && mutant->fitness <= POP_NOW[i_target]->fitness)) {
                //IL TRIAL RIMPIAZZA IL TARGET
                POP_NEW[i_target] = mutant;
                fitness_vector[i_target] = mutant->fitness;
                xb_vector[i_target] = mutant->XB;
                //jDE - aggiornameto f e CR che hanno avuto successo
                POP_NEW[i_target]->CR = CR;
                POP_NEW[i_target]->f = f;
                conteggio_successi_generazione_attuale++;
            } else {//IL TRIAL è SCARTATO
                for (row = 0; row < c; row++) {
                    free(mutant -> V_p[row]);
                    free(mutant -> U_p[row]);
                }
                free(mutant->V_p);
                free(mutant->U_p);
                free(mutant);

                if (abilita_invecchiamento) {
                    //INVECCHIAMENTO del sopravvissuto
                    aggiornaBestFitIndex();
                    if (i_target != bestFitIndex) { //se non è il migliore invecchia
                        if (POP_NOW[i_target] -> age > 0)//non sono in modalità reset
                            POP_NOW[i_target] -> age--;
                        if (POP_NOW[i_target] -> age <= 0) {//morte
                            printf("*");
                            //RINASCITA
                            //REINIT V_p
                            int rigaX;
                            for (i = 0; i < c; i++) {
                                if (!random_init)
                                    rigaX = random_at_most(n - 1);
                                for (j = 0; j < d; j++) {
                                    if (!random_init)
                                        POP_NOW[i_target] -> V_p[i][j] = X[rigaX][j] + drand48();
                                    else
                                        POP_NOW[i_target] -> V_p[i][j] = drand48();
                                    //in alternativa basato sul best fit?
                                }
                                //reset contatore membri del cluster
                                POP_NOW[i_target] -> V_p[i][d] = 0;
                                POP_NOW[i_target]->V_p[i][d + 1] = 0;
                            }
                            //SORT MATRICE V
                            sortMatrice(POP_NOW[i_target]->V_p);

                            //calcola U
                            for (i = 0; i < c; i++) {
                                for (j = 0; j < n; j++) {
                                    double denom = 0.0;
                                    double dist_x_j__v_i = calcDistanza(X[j], POP_NOW[i_target] -> V_p[i]);

                                    int k;
                                    for (k = 0; k < c; k++) {
                                        double dist_xj_vk = calcDistanza(X[j], POP_NOW[i_target] -> V_p[k]);
                                        denom += pow((dist_x_j__v_i / dist_xj_vk), esponente_U);
                                    }
                                    POP_NOW[i_target] -> U_p[i][j] = 1.0 / denom;
                                    if (POP_NOW[i_target] -> U_p[i][j] < 0) {
                                        printf("ERR: INIT POP:inizializzato U di un membro con negativo::");
                                        exit(-1);
                                    }
                                }
                            }

                            //calcolo fitness
                            double fit = calcolaFitness(POP_NOW[i_target]->V_p, POP_NOW[i_target]->U_p);
                            POP_NOW[i_target]->fitness = fit;
                            fitness_vector[i_target] = fit;

                            //calcolo XB
                            double xb = calcolaXB(POP_NOW[i_target]->V_p, POP_NOW[i_target]->U_p);
                            POP_NOW[i_target]->XB = xb;
                            xb_vector[i_target] = xb;

                            //impostazione età
                            POP_NOW[i_target] -> age = starting_age;
                        }
                    }
                }//END INVECCHIAMENTO

                //shake matrice target sopravvissuto alla selezione
                //e ricalcolo della sua matrice U

                if (abilita_shuffle) {
                    double probShuffle = dbl_rnd_inRange(0.0, 1.0);
                    if (probShuffle < 0.1) {
                        shuffleMatrice(POP_NOW[i_target]->V_p);
                        //ricalcola U
                        for (i = 0; i < c; i++) {
                            for (j = 0; j < n; j++) {
                                double denom = 0.0;
                                double dist_x_j__v_i = calcDistanza(X[j], POP_NOW[i_target] -> V_p[i]);

                                int k;
                                for (k = 0; k < c; k++) {
                                    double dist_xj_vk = calcDistanza(X[j], POP_NOW[i_target] -> V_p[k]);
                                    denom += pow((dist_x_j__v_i / dist_xj_vk), esponente_U);
                                }
                                POP_NOW[i_target] -> U_p[i][j] = 1.0 / denom;
                                if (POP_NOW[i_target] -> U_p[i][j] < 0) {
                                    printf("ERR: INIT POP:inizializzato U di un membro con negativo::");
                                    exit(-1);
                                }
                            }
                        }
                    }
                }

                //il target passa alla nuova generazione, il TRIAL era stato scartato
                POP_NEW[i_target] = POP_NOW[i_target];
            }//END SELECTION
        }//END DE//END GENERATION


        numero_generazioni--;




        if (conteggio_successi_generazione_attuale == 0)
            conteggio_reset++;

        if (abilita_reset) {
            if (conteggio_reset == reset_threshold) {
                puts("#RESET GLOBALE#");
                for (i = 0; i < num_pop; i++) {
                    aggiornaBestFitIndex();
                    if (i != bestFitIndex)
                        POP_NOW[i] -> age = 0;
                }
                conteggio_reset = 0; //reset contatore
            }
        }
        ultimo_conteggio_successi = conteggio_successi_generazione_attuale;
    } while (numero_generazioni > 0);
    //END GENERAZIONI

    puts("\n*END MAIN LOOP*");

    //ULTIMO SCAMBIO POPOLAZIONI
    for (i = 0; i < num_pop; i++) {
        if (POP_NEW[i] != POP_NOW[i]) {
            if (POP_NOW[i] != 0) {
                for (row = 0; row < c; row++) {

                    free(POP_NOW[i] -> V_p[row]);
                    free(POP_NOW[i] -> U_p[row]);
                }
                free(POP_NOW[i]->V_p);
                free(POP_NOW[i]->U_p);
                free(POP_NOW[i]);
            }
            POP_NOW[i] = POP_NEW[i];
        }
    }

    //computazione fitness della popolazione finale
    aggiornaBestFitIndex();
    aggiornaBestXBIndex();

    best_xb = POP_NOW[bestFitIndex]->XB;
    best_fit = POP_NOW[bestFitIndex]->fitness;
    printf("\n**MIGLIOR XB FINALE: \t%lf\n", best_xb);
    printf("\n**MIGLIOR FITNESS FINALE: \t%lf\n", best_fit);
    printf("\n**f: \t%lf\n", POP_NOW[bestFitIndex]->f);
    printf("\n**CR: \t%lf\n", POP_NOW[bestFitIndex]->CR);
    puts("matrice del best XB:");
    stampaMatrice(c, d, POP_NOW[bestXBIndex]->V_p);
    puts("matrice del best fitness:");
    stampaMatrice(c, d, POP_NOW[bestFitIndex]->V_p);
    stampaMatriceSuFile(c, d, POP_NOW[bestFitIndex]->V_p, out_V);
    stampaMatriceSuFile(c, n, POP_NOW[bestFitIndex]->U_p, out_U);
    puts("###################################################################");
    puts("###################################################################");
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

    int leggi_parametri_esterni = 1; //leggere parametri da CL
    if (leggi_parametri_esterni) {
        if (argc < 6) {
            //letture da utente
            printf("tipo dataset: ");
            scanf("%d", &tipo_dataset);
            printf("numero di punti in input: ");
            scanf("%d", &n);
            printf("numero di dimensioni: ");
            scanf("%d", &d);
            printf("numero di centroidi: ");
            scanf("%d", &c);
            printf("numero di generazioni: ");
            scanf("%d", &numero_generazioni);
        } else if (argc == 6) {
            tipo_dataset = atoi(argv[1]);
            n = atoi(argv[2]);
            d = atoi(argv[3]);
            c = atoi(argv[4]);
            numero_generazioni = atoi(argv[5]);
        } else {
            puts("SINTASSI: ./defc9b tipo_ds n d c num_gen");
            exit(-1);
        }
    } else {
        //test
        puts("!override parametri input attivo");
        n = 600;
        tipo_dataset = 0; //gauss = 0, s = 1
        d = 2;
        c = 15;
        numero_generazioni = 500;
    }

    //PARAMETRI INIZIALI
    m = 2.0; //fuzzification factor
    esponente_U = 2.0 / (m - 1.0);
    starting_age = numero_generazioni / 10; //timer iniziale
    abilita_invecchiamento = 0;
    abilita_reset = 0; //richiede invecchiamento
    reset_threshold = numero_generazioni / 2;
    abilita_partitioning = 0; //riodina vettori delle V secondo la prima coordinata
    abilita_shuffle = 0; //mescola i centroidi di V, con bassa probabilità (10%)
    attivaGnuPlot = 0;
    int output_csv = 1; //accende output su csv
    testLoadVIdeale = 0; //SOLO TEST, carica da file una matrice V predeterminata e la assegna al primo della popolazione
    random_init = 0; //se a 0 utilizza punti dell'input (con disturbo) per inizializzare
    usa_xb_per_fitness = 0; //DIVERGE
    usa_sumsep = 0; //richiede usa xb per fitness, usa somma delle distanza al denominatore di XB, DIVERGE
    attiva_partitioned_init = 0; //divide equamente in bins la posizione iniziale dei centroidi all'inizializzazione

    //per cluster sbilanciati
    aggiungi_peso_sigma = 1; //aumenta la sigma di una soluzione con il numero di punti appartenenti ai centroidi oltre la soglia indicata da soglia_conteggio
    attiva_sigma_separate = 1; //(richiede peso sigma) calcola sigma come somma delle diverse sigma dei cluster, il peso sarà dato in modo separato

    //peso_sigma = 0.5; //"compattazione" dei cluster molto grandi, sopra alla soglia_peso_sigma come numero di elementi conteggiati secondo soglia_conteggio
    soglia_conteggio = 0.90; //soglia di appartenenza di un punto per conteggiarlo come facente strettamente parte del cluster
    soglia_peso_sigma = n / c; //soglia (sul numero di punti nel cluster) oltre la quale si applica il peso sigma

    puts("v9b");
    numero_generazioni_iniziale = numero_generazioni;
    conteggio_reset = 0;
    //stream file
    //matrice di input
    out_X = fopen("dataset/aggregation.data", "r");
    //matrice di output centroidi
    out_V = fopen("v_defc9b.out", "w");
    //matrice output appartenenze
    out_U = fopen("u_defc9b.out", "w");
    if (testLoadVIdeale) {
        //matrice V per test funzione obiettivo
        in_V = fopen("v_test.dat", "r");
        puts("TEST LOAD V IDEALE ATTIVO!");
    }
    debug_log = fopen("debug/debug_defcv9b.debug", "w");

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

    if (random_init)
        srand48(time(NULL));
    init(n, c, d);
    lavora(n, c, d);
    plot();

    fclose(out_U);
    fclose(out_V);
    if (testLoadVIdeale)
        fclose(in_V);
    fclose(debug_log);

    //scrittura csv
    if (output_csv) {
        out_csv = fopen("csv/output_defcv9b.csv", "a");
        fprintf(out_csv, "dataset:%d,", tipo_dataset);
        fprintf(out_csv, "best_XB:%lf,", best_xb);
        fprintf(out_csv, "best_fit:%lf,", best_fit);
        fprintf(out_csv, "usa peso sigma:%d,", aggiungi_peso_sigma);
        fprintf(out_csv, "val peso sigma:%lf,", peso_sigma);
        fprintf(out_csv, "sog peso sigma:%lf,", soglia_peso_sigma);
        fprintf(out_csv, "sog count:%lf,", soglia_conteggio);
        fprintf(out_csv, "separa sigma:%d,", attiva_sigma_separate);
        fprintf(out_csv, "n:%d,", n);
        fprintf(out_csv, "c:%d,", c);
        fprintf(out_csv, "d:%d,", d);
        fprintf(out_csv, "num_fitness_eval:%d,", num_pop * numero_generazioni_iniziale);
        fprintf(out_csv, "start_age:%d,", starting_age);
        fprintf(out_csv, "aging:%d,", abilita_invecchiamento);
        fprintf(out_csv, "reset:%d,", abilita_reset);
        fprintf(out_csv, "random_init:%d,", random_init);
        fprintf(out_csv, "shuffle:%d\n", abilita_shuffle);

        fclose(out_csv);
    }



    return (0);
}

