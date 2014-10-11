#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>

double m = 2.0; //fuzzification
#define c 9 //numero di centri di cluster
#define n 200 //numero di punti totale in input

#define d 2 //dimensioni spaziali
double CR = 0.9; //crossover rate [0,1]
int numero_generazioni = 1000;//numero massimo di generazioni
double X[n][d]; //dati input

//parametro per l'arresto automatico
double guardia_xb = 0.001;

//parametri dell'inizializzazione dell'input
int mi_gauss = 4;
double sigma_gauss = 2.0;

//attiva e disattiva GnuPlot
int attivaGnuPlot = 1;

//individuo della popolazione
typedef struct el_pop {
    double V_p[c][d];
    double U_p[c][n];
    double indice_xb;
} el_pop;

/*
molt_pop moltiplicato per c numero di cluster
regola la grandezza della
popolazione
 */
int molt_pop = 10;
el_pop *POP_NEW[c * 10]; //VETTORE POPOLAZIONE NUOVA
el_pop *POP_NOW[c * 10]; //VETTORE POPOLAZIONE ATTUALE

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

///////////////////////////////////////////////////////////////////////////////

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

double calcolaXB(double V[c][d], double U[c][n], int debug) {
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
                if (dist_tmp == 0)
                    puts("dist_tmp nullo");
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
    
    
    
    if (min_sep == 0){
        puts("calcolaXB: DISTANZA NULLA MIN SEP");
        exit(-1);
    }
    if (sigma == 0){
        puts("calcolaXB: SIGMA NULLO");
        exit(-1);
    }
    
    double ris = (sigma)/(n*min_sep);
    
    if(ris <= 0){
        puts("calcolato XB nullo");
        exit(-1);
    }
        

    return ris;
}

int main(int argc, char** argv) {
    //PUNTATORI A FILE DI OUTPUT
    FILE *out_V, *out_X, *out_U;
    out_V = fopen("v.dat", "w");
    out_X = fopen("x.dat", "w");
    out_U = fopen("u.dat", "w");

    int i, j;
    //INIT X
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++)
            //gaussiana con media mi_gauss e devstd sigma_gauss
            X[i][j] = mi_gauss + (sigma_gauss * random_normal());
        if (i == 50) {
            mi_gauss *= 8;
        }
        if (i == 100) {
            mi_gauss *= 2;
        }
        if (i == 150) {
            mi_gauss *= 2;
        }
        
    }
    stampaMatriceSuFile(n, d, X, out_X);

    double esponente_U = 2.0 / (m - 1.0);

    //###INIT POPOLAZIONE 0
    //init V e calcolo U
    int pop_index;
    for (pop_index = 0; pop_index < c * molt_pop; pop_index++) {
        POP_NEW[pop_index] = malloc(sizeof (el_pop));
        //init V
        for (i = 0; i < c; i++)
            for (j = 0; j < d; j++) {
                POP_NEW[pop_index] -> V_p[i][j] = X[random_at_most(n - 1)][random_at_most(d - 1)] + 1;
                if (POP_NEW[pop_index] -> V_p[i][j] <= 0) {
                    printf("WARN: INIT POP:inizializzato V di un membro con negativo::%lf\n", POP_NEW[pop_index] -> V_p[i][j]);
                }
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

        //computazione XB della popolazione iniziale (POPOLAZIONE 0)
        double xb = calcolaXB(POP_NEW[pop_index]->V_p, POP_NEW[pop_index]->U_p, 0);
        if (xb <= guardia_xb) {
            puts("INIT POP: xb di un elemento inferiore a guardia_xb");
            stampaMatrice(c,d,POP_NEW[pop_index]->V_p);
            exit(-1);
        }
        POP_NEW[pop_index]->indice_xb = xb;
    }
    //###END INIT POPOLAZIONE 0
    printf("\n########## FINE INIT #############\n");

    double xb_selezionato;
    int numero_generazioni_utente = numero_generazioni;
    do {//NUOVA GENERAZIONE
        //printf("\n\nCOUNTDOWN GENERAZIONE: %d\n", numero_generazioni);
        //SCAMBIO VETTORI POPOLAZIONE
        //REINIT POP_NEW
        int i_CPop;
        for (i_CPop = 0; i_CPop < c * molt_pop; i_CPop++) {
            if (POP_NEW[i_CPop] != POP_NOW[i_CPop]) {
                if(POP_NOW[i_CPop != 0])
                    free(POP_NOW[i_CPop]);
                POP_NOW[i_CPop] = POP_NEW[i_CPop];
            }
        }

        /////////////////DE////////////////////
        for (i_CPop = 0; i_CPop < c * molt_pop; i_CPop++) {//PER OGNI COMPONENTE DELLA POP
            //SCELTA CANDIDATI
            //tre vettori devono essere scelti a caso nella popolazione
            //diversi dal target (indice i) e mutualmente
            int indice_1, indice_2, indice_3;
            do {
                indice_1 = random_at_most(((long) c * molt_pop) - 1);
            } while (indice_1 == i_CPop);

            do {
                indice_2 = random_at_most(((long) c * molt_pop) - 1);
            } while (indice_2 == i_CPop || indice_2 == indice_1);

            do {
                indice_3 = random_at_most(((long) c * molt_pop) - 1);
            } while (indice_3 == i_CPop || indice_3 == indice_1 || indice_3 == indice_2);

            if (indice_1 == i_CPop || indice_2 == i_CPop || indice_3 == i_CPop)
                puts("INDICE NON VALIDO!!!!");

            //l'elemento mutante
            el_pop *trial = malloc(sizeof (el_pop));
            //MUTATION
            //scambio di geni
            int i1, j1;
            for (i1 = 0; i1 < c; i1++) {
                for (j1 = 0; j1 < d; j1++) {
                    double f = fRand(0.0, 1.0);
                    trial->V_p[i1][j1] = POP_NOW[indice_3]->V_p[i1][j1] + f * (POP_NOW[indice_1]->V_p[i1][j1] - POP_NOW[indice_2]->V_p[i1][j1]);
                }
            }

            //CROSSOVER CON IL VETTORE ATTUALE (TIPO 1)
            for (i1 = 0; i1 < c; i1++) {
                for (j1 = 0; j1 < d; j1++) {
                    double prob_crossover = fRand(0.0, 1.0);
                    if (prob_crossover < CR) {
                        //prendo il gene del vettore attuale
                        trial->V_p[i1][j1] = POP_NOW[i_CPop]->V_p[i1][j1];
                    }
                }
            }

            //CROSSOVER CON IL VETTORE ATTUALE (TIPO 2)
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
                    double dist_x_j__v_i = calcDistanza(X[j], trial -> V_p[i]);
                    if (dist_x_j__v_i == 0) {
                        puts("calcolo U mutante, distanza nulla");
                        exit(-1);
                    }
                    //printf("\ndist_x_j__v_i: %lf,%d,%d\n",dist_x_j__v_i,c,n);
                    int k;
                    for (k = 0; k < c; k++) {
                        double dist_xj_vk = calcDistanza(X[j], trial -> V_p[k]);
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
                    trial -> U_p[i][j] = 1.0 / denom;
                }
            }

            //calcolo XB mutante            
            trial->indice_xb = calcolaXB(trial->V_p, trial->U_p, 1);
            if (trial -> indice_xb <= guardia_xb) {
                puts("*ERROR: indice_xb trial inferiore a guardia, scartato");
                trial->indice_xb = DBL_MAX;
            }
            //////END CALCOLO FITNESS MUTANTE

            //SELECTION
            if (trial->indice_xb < POP_NOW[i_CPop]->indice_xb) {
                xb_selezionato = trial->indice_xb;
                POP_NEW[i_CPop] = trial;
            } else {
                free(trial);
                POP_NEW[i_CPop] = POP_NOW[i_CPop];
                xb_selezionato = POP_NOW[i_CPop]->indice_xb;
            }
        }//END DE
        numero_generazioni--;
    } while (numero_generazioni > 0);
    
    
    //computazione XB della popolazione finale
    double best_xb = DBL_MAX;
    int indice_best;
    for (pop_index = 0; pop_index < c * molt_pop; pop_index++) {
        POP_NEW[pop_index]->indice_xb = calcolaXB(POP_NEW[pop_index]->V_p, POP_NEW[pop_index]->U_p, 0);
        if (POP_NEW[pop_index]->indice_xb < best_xb) {
            best_xb = POP_NEW[pop_index]->indice_xb;
            indice_best = pop_index;
        }
    }


    puts("***********************************************");
    if(xb_selezionato < guardia_xb)
        puts("***SUPERATA GUARDIA XB***");
    int numero_gen_fatte = numero_generazioni_utente - numero_generazioni;
    printf("\nnumero generazioni fatte: %d",numero_gen_fatte);
    printf("\nmiglior XB:%lf\n", best_xb);
    puts("matrice V:");
    stampaMatrice(c, d, POP_NEW[indice_best]->V_p);
    puts("");
    stampaMatriceSuFile(c, d, POP_NEW[indice_best]->V_p, out_V);
    stampaMatriceSuFile(c, n, POP_NEW[indice_best]->U_p, out_U);
    puts("***********************************************");

    //GNUPLOT
    if ((d == 2 || d == 3) && attivaGnuPlot) {
        char *commandsForGnuplot[] = {"set key off","set term wxt 1","set title \"matrice X\"", "", "set term wxt 2","set key off", "set title \"matrice V\"", ""};
        if (d == 2) {
            commandsForGnuplot[3] = "plot 'x.dat'";
            commandsForGnuplot[7] = "plot 'v.dat'";
        } else if (d == 3) {
            commandsForGnuplot[3] = "splot 'x.dat'";
            commandsForGnuplot[7] = "splot 'v.dat'";
        }

        FILE * gnuplotPipe = popen("gnuplot -persistent", "w");

        for (i = 0; i < 8; i++) {
            fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);
            fflush(gnuplotPipe);
        }

    }



    return (0);
}

