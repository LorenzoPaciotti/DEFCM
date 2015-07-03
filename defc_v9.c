#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#include <string.h>

FILE *out_V, *out_X, *out_U, *out_LOG, *out_LOG_RIS, *out_csv, *debug_log, *dataset; //puntatori a file di output

double **X;
double m = 2.0; //fuzzification factor
double CR; //crossover rate
int tipo_dataset;
int starting_age;
int ultimo_conteggio_successi;
int conteggio_successi_generazione_attuale;
int conteggio_adattamenti;
int reset_threshold;
int bestFitIndex;
double best_xb;
double best_fit;
int num_pop_iniziale, numero_generazioni, i, j, k, pop_index, numero_generazione_attuale, numero_generazioni_iniziale;
double esponente_U;
int n, c, d; //numero di punti in input, numero di centroidi, numero di dimensioni
int abilita_partitioning; //riordina gli array V per far accoppiare solo i centroidi della stessa zona dello spazio
int abilita_invecchiamento;
int usa_xb_per_fitness;
int abilita_shuffle;
int attivaGnuPlot; //attiva e disattiva GnuPlot

//numero di elementi della popolazione - fare parametrico
#define num_pop 100 // 50, 100

//struttura elemento popolazione

typedef struct el_pop {
    double **V_p;
    double **U_p;
    double fitness;
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

int rand_int(int n) {
    int limit = RAND_MAX - RAND_MAX % n;
    int rnd;

    do {
        rnd = rand();
    } while (rnd >= limit);
    return rnd % n;
}

void shuffleArray(int *array, int n) {
    int i, j, tmp;

    for (i = n - 1; i > 0; i--) {
        j = rand_int(i + 1);
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
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
    return sigma / (n * min_sep);
}

double calcolaFitness(double **V, double **U, int debug) {
    if (usa_xb_per_fitness) {
        return calcolaXB(V, U, 0);
    } else {
        //CALCOLO SIGMA
        double sigma = 0.0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < c; j++) {
                sigma += pow(U[j][i], m) * pow(calcDistanza(V[j], X[i]), 2.0);
            }
        }
        return sigma;
    }

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

void aggiornaBestFitIndex() {
    double bestFit = DBL_MAX;
    for (i = 0; i < num_pop; i++) {
        if (fitness_vector[i] < bestFit) {
            bestFit = fitness_vector[i];
            bestFitIndex = i;
        }
    }
}



//INIZIALIZZAZIONE GENERALE DI TUTTA LA POPOLAZIONE INIZIALE

void init(int n, int c, int d) {
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
            //int rigaX = random_at_most(n-1);
            for (j = 0; j < d; j++) {

                //POP_NEW[pop_index] -> V_p[i][j] = X[rigaX][j] + drand48();
                //POP_NEW[pop_index] -> V_p[i][j] = X[random_at_most(n - 1)][random_at_most(d - 1)] + drand48();
                //POP_NEW[pop_index] -> V_p[i][j] = X[random_at_most(n - 1)][j] + drand48();
                POP_NEW[pop_index] -> V_p[i][j] = drand48();
                //POP_NEW[pop_index] -> V_p[i][j] = dbl_rnd_inRange(0,range_init_max);
                //POP_NEW[pop_index] -> V_p[i][j] = random_at_most(range_init_max);
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
        double fit = calcolaFitness(POP_NEW[pop_index]->V_p, POP_NEW[pop_index]->U_p, 0);
        POP_NEW[pop_index]->fitness = fit;
        fitness_vector[pop_index] = fit;


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
        for (i = 0; i < num_pop_iniziale; i++) {
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

        //DEBUG/////////////////////////
        if (numero_generazione_attuale % 50 == 0) {
            aggiornaBestFitIndex();
            double best_xb = calcolaXB(POP_NOW[bestFitIndex]->V_p, POP_NOW[bestFitIndex]->U_p, 0);
            double best_fit = POP_NOW[bestFitIndex]->fitness;
            printf("generazione:%d  best XB:%lf  best fit:%lf\n", numero_generazioni, best_xb, best_fit);
        }
        /////////////////////////////

        /////////////////DE////////////////////
        for (i_target = 0; i_target < num_pop_iniziale; i_target++) {//PER OGNI COMPONENTE DELLA POP
            //SCELTA CANDIDATI
            //tre vettori devono essere scelti a caso nella popolazione
            //diversi dal target (indice i) e mutualmente

            int indice_1, indice_2, indice_base;
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
            //la f cambia da vettore centroide a vettore centroide
            //fattore scala mutazione
            for (i1 = 0; i1 < c; i1++) {//centroide
                for (j1 = 0; j1 < d; j1++) {//dimensione
                    mutant->V_p[i1][j1] = POP_NOW[indice_base]->V_p[i1][j1] + f * (POP_NOW[indice_1]->V_p[i1][j1] - POP_NOW[indice_2]->V_p[i1][j1]);
                }
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

            //CROSSOVER SPLIT
            /*for (i1 = 0; i1 < c; i1++) {//centroide
                for (j1 = 0; j1 < d; j1++) {//coordinata centroide
                    if (i1 >= floor((double) c / 2)) {//prendo i cromosomi dal target so oltre il punto di CO
                        mutant->V_p[i1][j1] = POP_NOW[i_target]->V_p[i1][j1];
                    }
                }
            }*/

            //SORT MATRICE MUTANTE (ORA DETTO TRIAL)
            sortMatrice(mutant->V_p);
            
            //shake matrice mutante
            //shuffle(mutant->V_p,c);

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
            mutant->fitness = calcolaFitness(mutant->V_p, mutant->U_p, 1);

            //impostazione timer età trial
            mutant -> age = starting_age;

            //SELECTION
            if (mutant->fitness < POP_NOW[i_target]->fitness) {//IL TRIAL RIMPIAZZA IL TGT
                POP_NEW[i_target] = mutant;
                fitness_vector[i_target] = mutant->fitness;
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
                    aggiornaBestFitIndex(); //MIGLIORE
                    if (i_target != bestFitIndex) { //se non è il migliore invecchia
                        if (POP_NOW[i_target] -> age > 0)//non sono in modalità reset
                            POP_NOW[i_target] -> age--;
                        if (POP_NOW[i_target] -> age <= 0) {//morte
                            //RINASCITA
                            //REINIT V_p
                            for (i = 0; i < c; i++) {
                                //int rigaX = random_at_most(n-1);
                                //la mutazione è usata solo con rinascita basata su best fit
                                //double mutazione = drand48() * intens_mutazione_neonato;
                                for (j = 0; j < d; j++) {
                                    //POP_NOW[i_target] -> V_p[i][j] = POP_NOW[bestFitIndex]->V_p[i][j] + drand48();
                                    //POP_NOW[i_target] -> V_p[i][j] = X[random_at_most(n-1)][j]+drand48();
                                    //POP_NOW[i_target] -> V_p[i][j] = X[random_at_most(n - 1)][random_at_most(d - 1)] + drand48();
                                    //POP_NOW[i_target] -> V_p[i][j] = lrand48();
                                    POP_NOW[i_target] -> V_p[i][j] = drand48();
                                }
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
                            double fit = calcolaFitness(POP_NOW[i_target]->V_p, POP_NOW[i_target]->U_p, 0);
                            POP_NOW[i_target]->fitness = fit;
                            fitness_vector[i_target] = fit;

                            //impostazione età
                            POP_NOW[i_target] -> age = starting_age;
                        }
                    }
                }//END INVECCHIAMENTO

                //il target passa alla nuova generazione, il TRIAL era stato scartato
                POP_NEW[i_target] = POP_NOW[i_target];
            }//END SELECTION
        }//END DE//END GENERATION

        numero_generazioni--;
        ultimo_conteggio_successi = conteggio_successi_generazione_attuale;
    } while (numero_generazioni > 0);
    //END GENERAZIONI

    puts("\n*END MAIN LOOP*");

    //ULTIMO SCAMBIO POPOLAZIONI
    for (i = 0; i < num_pop_iniziale; i++) {
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

    //calcolo xb del migliore
    best_xb = calcolaXB(POP_NOW[bestFitIndex]->V_p, POP_NOW[bestFitIndex]->U_p, 0);

    best_fit = POP_NOW[bestFitIndex]->fitness;

    printf("\n**MIGLIOR XB FINALE: \t%lf\n", best_xb);
    printf("\n**MIGLIOR FITNESS FINALE: \t%lf\n", best_fit);
    stampaMatrice(c, d, POP_NOW[bestFitIndex]->V_p);
    stampaMatriceSuFile(c, d, POP_NOW[bestFitIndex]->V_p, out_V);
    stampaMatriceSuFile(c, n, POP_NOW[bestFitIndex]->U_p, out_U);
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
    esponente_U = 2.0 / (m - 1.0);
    num_pop_iniziale = num_pop;

    //PARAMETRI INIZIALI
    starting_age = 10;
    abilita_partitioning = 1;
    abilita_invecchiamento = 1;
    abilita_shuffle = 1;
    usa_xb_per_fitness = 0; //mai funzionato
    attivaGnuPlot = 0;
    int output_csv = 0; //accende output su csv
    int leggi_parametri_esterni = 0; //leggere parametri da CL

    if (leggi_parametri_esterni) {
        if (argc == 0) {
            //letture da utente
            //printf("tipo dataset: ");
            //scanf("%d", &tipo_dataset);
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
            puts("defc9 tipo_ds n d c num_gen");
            exit(0);
        }
    } else {
        //test
        puts("!override parametri input attivo");
        n = 10000;
        tipo_dataset = 0;
        d = 2;
        c = 4;
        numero_generazioni = 100;
    }

    puts("v9: jDE");
    numero_generazioni_iniziale = numero_generazioni;

    //stream file
    //matrice di input
    out_X = fopen("dataset/gauss4.data", "r");
    //matrice di output centroidi
    out_V = fopen("v_defc9.dat", "w");
    //matrice output appartenenze
    out_U = fopen("u_defc9.dat", "w");
    //out_LOG_RIS = fopen("log_ris", "a");
    //debug_log = fopen("debug", "w");
    //fputs("\n######\n", out_LOG_RIS); //nuova log entry

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

    //scrittura csv
    if (output_csv) {
        out_csv = fopen("csv/output_defcv9.csv", "a");
        fprintf(out_csv, "%d,", tipo_dataset);
        fprintf(out_csv, "%d,", n);
        fprintf(out_csv, "%d,", c);
        fprintf(out_csv, "%d,", d);
        fprintf(out_csv, "%d,", num_pop);
        fprintf(out_csv, "%d,", numero_generazioni_iniziale);
        fprintf(out_csv, "%d,", starting_age);
        fprintf(out_csv, "%lf\n", best_xb);
        fclose(out_csv);
    }



    return (0);
}

