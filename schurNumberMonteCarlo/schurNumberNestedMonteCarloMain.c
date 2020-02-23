//
//  schurNumberNestedMonteCarloMain.c
//  schurNumberMonteCarlo
//
//  Created by Weak Schur Number CS on 22/04/2019.
//  Ce fichier contient la fonction main qui permet de construire un exécutable.
//

#include "asprintf.h"
#include <stdio.h>
#include <unistd.h>
#include <libgen.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "schurNumberNestedMonteCarloHeader.h"

void usage(char *cmdname) {
    fprintf(stderr,
            "usage: %s [-h] [-l level] [-s simulnum] [-i iternum] [-m method] [-t time] [-pv] [-b filename] [-d filename] [-o filename] partnumber\n"\
            "\t-b filename : Specify a filename for an initial partition\n"\
            "\t-d filename : Create a file named filename where to save the generated numbers.\n"\
            "\t-o filename : Create a file named filename to save the best found partition.\n"\
            "\t-p : Print the found partition to stdout.\n"\
            "\t-v : Print mean and standard deviation to stdout.\n"\
            "\t-m method : Specify the method to use.\n"\
            "\t\t0 : Estimate the schur number using a simple action (integer by integer). This is the default method.\n"\
            "\t\t1 : Estimate the weak schur number using a simple action (integer by integer).\n"\
            "\t-l level: Specify the nesting level. By default level = 1.\n"\
            "\t-s simulnum: Specify the number of simulation at the top level. By default simulnum = 1.\n"\
            "\t-i iternum: Specify the number of iteration at level 0. By default iternum = 64.\n"\
            "\t-t time: Specify a time (in minutes) after which no more simulations are launched.\n"\
            "\t-h: Print usage message.\n",
            cmdname);
}

void printPartition(unsigned int p, unsigned long n, mp_limb_t **partition) {
    /*Affiche une partition.*/
    unsigned long limbn;
    mp_size_t limbsize;
    unsigned int i;
    mp_bitcnt_t j;
    mp_limb_t *set;
    mp_limb_t limb;
    limbsize = (n>>6) + 1;
    printf("Partition:\n");
    for (i=0; i<p; i++) {
        printf("\t");
        set = partition[i];
        for (limbn=0; limbn<limbsize; limbn++) {
            limb = set[limbn];
            for (j=0; j<mp_bits_per_limb; j++) {
                if (limb & ((mp_limb_t)1<<j)) {
                    printf(" %lu", limbn*mp_bits_per_limb + j + 1);
                }
            }
        }
        printf("\n");
    }
}

unsigned int makeInitialPartition(char *filename, partition_t *partitionbeginstruc) {
    /*Si filename != NULL, lit le fichier filename et initie partitionbeginstruc à la partition
        contenue dans le fichier. Si filename == NULL, partitionbeginstruc = {{1}}.
     Renvoie la taille de la partition en l'absence de problème, 0 sinon.*/
    unsigned int p;
    
    if (filename) {
        /*Débuter à partir d'une partition contenue dans le fichier filename.*/
        p = schurNumberScanPartitionFromFile(filename, partitionbeginstruc);
    } else {
        /*Prendre pour partition initiale {{1}}*/
        p = 1;
        partition_init(1, 1, partitionbeginstruc);
        partitionbeginstruc->p = 1;
        partitionbeginstruc->limbsize = 1;
        partitionbeginstruc->n = 1;
        **(partitionbeginstruc->partition) = (mp_limb_t)1;
        **(partitionbeginstruc->partitioninvert) = (mp_limb_t)1<<(mp_bits_per_limb - 1);
    }
    
    return p;
}

unsigned long countOccurence(unsigned long *distribution, size_t size, unsigned long n) {
    /*Renvoie le nombre d'occurences de n dans distribution, tableau de taille size.*/
    unsigned long count = 0;
    unsigned long i;
    unsigned long *ptr;
    
    ptr = distribution;
    
    for (i=0; i<size; i++) {
        if (*ptr == n) {
            count++;
        }
        ptr++;
    }
    
    return count;
}

char saveDistribution(char *filename, unsigned long *distribution, size_t size, unsigned long nmax) {
    /*Enregistre les n entiers de distribution dans un fichier au format tsv.
       Le format est entier_n\tnombre_d'occurences\n.
     Renvoie 1 si tout s'est bien passé, 0 sinon.*/
    FILE *fp;
    unsigned long n;
    unsigned long count, total;
    
    fp = fopen(filename, "w");
    if (fp) {
        n = nmax;
        total = 0;
        while (total < size) {
            count = countOccurence(distribution, size, n);  // Compte le nombre de fois qu'une partition de [1, n] a été obtenue
            total += count;
            if (count) {
                /*Au moins une partition sans-somme de [1, n] a été obtenue.*/
                fprintf(fp, "%lu\t%lu\r\n", n, count);
            }
            n--;
        }
        fclose(fp);
        return 1;
    }
    
    return 0;
}

void printStatistics(unsigned long *distribution, size_t size) {
    /*Calcule la valeur moyenne et l'écart-type grâce à la méthode de Welford. Calcule aussi un intervalle de confiance à 95%, i.e. le plus grand intervalle centré sur la moyenne I = [<x> - a; <x> + a] de sorte que P(I) ≤ 95%.
     
     Les résultats sont affichés dans stdout.*/
    unsigned long i;
    double mean, var, delta;
    unsigned long ninf, nsup;   // Borne de l'intervalle de confiance
    size_t icount;              // Nombre d'entiers appartenant à l'intervalle de confiance
    size_t count95;             // count * 95 %
    
    mean = 0;
    var = 0;
    for (i=0; i < size; i++) {
        delta = (double)distribution[i] - mean;
        mean += delta/(i+1);
        var += delta * ((double)distribution[i] - mean);
    }
    var = var / (size-1);
    
    ninf = floor(mean);
    nsup = ceil(mean);
    count95 = floor(size * 0.95);
    icount = countOccurence(distribution, size, ninf) + countOccurence(distribution, size, nsup);
    while (icount <= count95) {
        ninf--;
        nsup++;
        icount += countOccurence(distribution, size, ninf) + countOccurence(distribution, size, nsup);
    }
    ninf++;
    nsup--;
    
    printf("Moyenne : %f\nEcart-type : %f\nIntervalle de confiance : [%lu; %lu]\n", mean, sqrt(var), ninf, nsup);
}

unsigned long schurNumberNestedMonteCarlo(unsigned int p, unsigned long *narray, unsigned int level, unsigned int simulnum, unsigned int iternum, char method, mp_limb_t **sfpartitionbestglobal, partition_t *sfpartitionbeginstruc, clock_t time) {
    /*Fonction à lancer pour trouver une borne inférieure au nombre de Schur S(p).
     Elle alloue la mémoire nécessaire aux simulations et lance simulnum simulations de niveau level.
     Les simulations de niveau 0 seront effectués simulnum0 fois.
     
     La variable méthode permet de sélectionner la méthode à employer. La variable time permet de spécifier une durée au bout de laquelle plus aucune simulation n'est lancée.*/
    clock_t clock0;
    unsigned long nbest, nsimulated;
    unsigned int i, j, pbest;
    mp_size_t limballoc;
    partition_t sfpartitionstruc;
    mp_limb_t **sfpartitionbest;
    
    /*Initialisation*/
    limballoc = (1<<(2*p))>>mp_bits_per_limb;
    partition_init(p, mp_bits_per_limb * limballoc, &sfpartitionstruc);
    sfpartitionbest = calloc(p, sizeof(mp_limb_t *));       //Tableau servant au niveau en-dessous
    for (i=0; i<p; i++) {
        sfpartitionbest[i] = calloc(limballoc, sizeof(mp_limb_t));
    }
    
    nbest = 1;
    clock0 = clock();
    for (i=0; i<simulnum; i++) {
        /*Création de la première partition à partir de sfpartitionbegin*/
        partition_copy(&sfpartitionstruc, sfpartitionbeginstruc);
        
        /*Lancement de la simulation de niveau level*/
        switch (method) {
            case '1':
                nsimulated = schurNumberWeakSimpleMonteCarloLevelIteration(&sfpartitionstruc, level, sfpartitionbest, &pbest, iternum);
                break;
                
            default:
                nsimulated = schurNumberSimpleMonteCarloLevelIteration(&sfpartitionstruc, level, sfpartitionbest, &pbest, iternum);
                break;
        }
        narray[i] = nsimulated + 1;
        if (nsimulated > nbest) {
            nbest = nsimulated;
            for (j=0; j<pbest; j++) {
                mpn_copyd(sfpartitionbestglobal[j], sfpartitionbest[j], 1 + (nbest / GMP_NUMB_BITS));
            }
        }
        if (clock() - clock0 > time) {
            fprintf(stderr, "Process has lasted too much time. It has been interrupted after %u simulations.\n", i+1);
            i = simulnum;
        }
    }
    
    /*Nettoyage*/
    partition_unalloc(&sfpartitionstruc);
    for (i=0; i<p; i++) {
        free(sfpartitionbest[i]);
    }
    free(sfpartitionbest);
    
    return nbest + 1;
}

int main(int argc, const char * argv[]) {
    char optc;
    char *cmdname;
    unsigned int p, i;
    unsigned int level;
    unsigned int simulnum;
    unsigned int iternum;
    char printpartition;
    char statistics;
    char method;
    char problem;
    clock_t cputimemax;
    unsigned long *narray;
    unsigned long nmax;
    char *bfilename;
    char *dfilename;
    char *ofilename;
    mp_limb_t **sfpartitionbestglobal;
    partition_t partitionbeginstruc;
    
    /*Set variables to default*/
    cmdname = basename(argv[0]);
    level = 1;
    simulnum = 1;
    iternum = 64;
    printpartition = 0;
    statistics = 0;
    method = '0';
    problem = 0;
    cputimemax = ULLONG_MAX;
    bfilename = NULL;
    dfilename = NULL;
    ofilename = NULL;
    
    while ((optc = getopt(argc, argv, "hpvl:s:i:b:d:o:m:t:")) != -1) {
        /*Parse arguments*/
        switch (optc) {
                
            case 'b':
                asprintf(&bfilename, "%s", optarg);
                break;
                
            case 'd':
                asprintf(&dfilename, "%s", optarg);
                break;
                
            case 'o':
                asprintf(&ofilename, "%s", optarg);
                break;
                
            case 'h':
                usage(cmdname);
                break;
                
            case 'p':
                printpartition = 1;
                break;
            
            case 'v':
                statistics = 1;
                break;
                
            case 'm':
                method = *optarg;
                if (method != '0' && method != '1') {
                    problem = 1;
                    fprintf(stderr, "%s: -m %c: unknown method\n", cmdname, method);
                    usage(cmdname);
                }
                break;
                
            case 'l':
                level = atoi(optarg);
                break;
                
            case 's':
                simulnum = atoi(optarg);
                break;
                
            case 'i':
                iternum = atoi(optarg);
                break;
                
            case 't':
                cputimemax = atol(optarg) * CLOCKS_PER_SEC * 60;
                break;
                
            default:
                fprintf(stderr, "%s: -%c: invalid option\n", cmdname, optc);
                usage(cmdname);
                return 0;
                break;
        }
    }
    
    if (argc - optind <= 0) {
        fprintf(stderr, "%s: %s: invalid argument\n", cmdname, argv[optind]);
        usage(cmdname);
        return 0;
    }
    p = atoi(argv[optind]);
    if (p == 0) {
        fprintf(stderr, "%s: %s: invalid argument\n", cmdname, argv[optind]);
        usage(cmdname);
        return 0;
    }
    
    i = makeInitialPartition(bfilename, &partitionbeginstruc);
    if (!i) {
        /*Problème dans l'ouverture du fichier bfilename*/
        fprintf(stderr, "%s: unable to open file %s.\n", cmdname, bfilename);
        problem = 1;
    }
    if (i > p) {
        fprintf(stderr, "%s: the partition in file %s contains more sets (%u sets) than specified (%u sets). Unable to do any computations.\n", cmdname, bfilename, i, p);
        problem = 1;
    }
    if (bfilename) {
        free(bfilename);
    }
    if (dfilename) {
        i = saveDistribution(dfilename, NULL, 0, 0);
        if (!i) {
            fprintf(stderr, "%s: unable to open %s for writing.\n", cmdname, dfilename);
            free(dfilename);
            problem = 1;
        }
    }
    if (ofilename) {
        i = schurNumberPrintPartitionToFile(ofilename, 0, 1, NULL);
        if (!i) {
            fprintf(stderr, "%s: unable to open %s for writing.\n", cmdname, ofilename);
            free(ofilename);
            problem = 1;
        }
    }
    
    if (!problem) {
        /*Allocation des variables.*/
        narray = calloc(simulnum, sizeof(unsigned long));
        sfpartitionbestglobal = calloc(p, sizeof(mp_limb_t *)); //Contiendra la plus grande partition sans-somme trouvée.
        for (i=0; i<p; i++) {
            //sfpartitionbestglobal[i] = calloc(p, sizeof(mp_limb_t));
            sfpartitionbestglobal[i] = calloc((1<<(2*p))>>mp_bits_per_limb, sizeof(mp_limb_t));
        }
        
        nmax = schurNumberNestedMonteCarlo(p, narray, level, simulnum, iternum, method, sfpartitionbestglobal, &partitionbeginstruc, cputimemax);
        
        if (method == '0') {
            printf("Schur Number S(%u) ≥ %lu\n", p, nmax);
        } else {
            printf("Weak Schur Number WS(%u) ≥ %lu\n", p, nmax);
        }
        
        if (printpartition) {
            /*Afficher la meilleure partition sans-somme trouvée.*/
            printPartition(p, nmax, sfpartitionbestglobal);
        }
        
        if (statistics) {
            printStatistics(narray, simulnum);
        }
        
        if (dfilename) {
            i = saveDistribution(dfilename, narray, simulnum, nmax);
            if (!i) {
                fprintf(stderr, "%s: unable to write in file %s.\n", cmdname, dfilename);
            }
            free(dfilename);
        }
        
        if (ofilename) {
            i = schurNumberPrintPartitionToFile(ofilename, p, nmax, sfpartitionbestglobal);
            if (!i) {
                fprintf(stderr, "%s: unable to write in file %s.\n", cmdname, ofilename);
            }
            free(ofilename);
        }
        
        /*Déallocation des variables.*/
        free(narray);
        for (i=0; i<p; i++) {
            free(sfpartitionbestglobal[i]);
        }
        free(sfpartitionbestglobal);
        partition_unalloc(&partitionbeginstruc);
    }
    
    return 0;
}
