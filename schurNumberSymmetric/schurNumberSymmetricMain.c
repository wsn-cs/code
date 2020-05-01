//
//  schurNumberSymmetricMain.c
//  schurNumberSymmetric
//
//  Created by wsn-cs on 30/05/2019.
//

#include "schurNumberSymmetricHeader.h"
#include <stdio.h>
#include <unistd.h>
#include <libgen.h>
#include <time.h>

void usage(char *cmdname) {
    fprintf(stderr,
            "usage: %s [-h] [-l level] [-s simulnum] [-i iternum] [-m method] [-t time] [-apv] [-b filename] [-o filename] [-d mindepth] partnumber interval\n"\
            "\t-a: Perform a complete search to find the greatest depth.\n"\
            "\t-b filename : Specify a filename for an initial partition.\n"\
            "\t-o filename : Create a file named filename to save the best found partition.\n"\
            "\t-p : Print the found partition to stdout.\n"\
            "\t-d mindepth: Specify a minimal depth. If multiple depths are given, they are sorted and used to majorate each set.\n"\
            "\t-m method : Specify the method to use.\n"\
            "\t\t0 : Estimate the schur number using a simple action (integer by integer). This is the default method.\n"\
            "\t\t1 : Estimate the weak schur number using a simple action (integer by integer).\n"\
            "\t-t time: Specify a time (in minutes) after which no more simulations are launched.\n"\
            "\t-h: Print usage message.\n",
            cmdname);
}

void sort(unsigned long *array, size_t count) {
    unsigned long min, a;
    size_t i, imin, c;
    unsigned long *ptr;
    
    ptr = array;
    
    for (c = count; c > 0; c--) {
        
        min = ULONG_MAX;
        imin = 0;
        
        for (i=0; i < c; i++) {
            a = ptr[i];
            if (a < min) {
                imin = i;
                min = a;
            }
        }
        
        ptr[imin] = ptr[0];
        ptr[0] = min;
        
        ptr++;
    }
}

void printSymmetricPartition(unsigned int p, unsigned long n, mp_limb_t **partition) {
    /*Affiche une partition.*/
    unsigned long limbn;
    mp_size_t limbsize;
    unsigned int i;
    mp_bitcnt_t j;
    mp_limb_t *set;
    mp_limb_t limb;
    
    limbsize = ((n+1) >> 1) / mp_bits_per_limb + 1;
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
        for (limbn=limbsize; limbn>0; limbn--) {
            limb = set[limbn-1];
            for (j=mp_bits_per_limb; j>0; j--) {
                if (limb & ((mp_limb_t)1<<(j-1))) {
                    printf(" %lu", n - ((limbn-1)*mp_bits_per_limb + (j-1)));
                }
            }
        }
        printf("\n");
    }
}

void makeInitialPartition(unsigned int p, unsigned long n, partition_t *partitionstruc) {
    unsigned int i;
    mp_size_t limballoc;
    mp_limb_t **partition, **partitioninvert;
    
    limballoc = 1 + (((n+1) >> 1) / mp_bits_per_limb);
    
    partition = calloc(p, sizeof(mp_limb_t *));
    partitioninvert = calloc(p, sizeof(mp_limb_t *));
    for (i=0; i<p; i++) {
        partition[i] = calloc(limballoc, sizeof(mp_limb_t));
        partitioninvert[i] = calloc(limballoc, sizeof(mp_limb_t));
    }
    
    **partition = (mp_limb_t)1;
    (*partitioninvert)[limballoc - 1] = (mp_limb_t)1<<(mp_bits_per_limb-1);
    
    partitionstruc->pmax = p;
    partitionstruc->p = 1;
    partitionstruc->n = 1;
    partitionstruc->limballoc = limballoc;
    partitionstruc->limbsize = 1;
    partitionstruc->partition = partition;
    partitionstruc->partitioninvert = partitioninvert;
    partitionstruc->work0 = calloc(limballoc, sizeof(mp_limb_t));
    partitionstruc->work1 = calloc(limballoc, sizeof(mp_limb_t));
}

void partition_unalloc(partition_t *partitionstruc) {
    /*Libère tous les grands entiers associés aux ensembles de partitionstruc.*/
    unsigned i, pmax;
    mp_limb_t **partition;
    mp_limb_t **partitioninvert;
    
    pmax = partitionstruc->pmax;
    partition = partitionstruc->partition;
    partitioninvert = partitionstruc->partitioninvert;
    for (i=0; i<pmax; i++) {
        free(*partition);
        free(*partitioninvert);
        partition++;
        partitioninvert++;
    }
    free(partitionstruc->partition);
    free(partitionstruc->partitioninvert);
    free(partitionstruc->work0);
    free(partitionstruc->work1);
}

int main(int argc, const char * argv[]) {
    char optc;
    char *cmdname;
    char printpartition;
    char completesearch;
    char method;
    char problem;
    clock_t cputimemax;
    unsigned int p;
    unsigned long n, d;
    unsigned long *depths;
    size_t depthcount;
    char *bfilename;
    char *dfilename;
    char *ofilename;
    partition_t partitionstruc;
    
    /*Get sets'number*/
    p = atoi(argv[argc - 2]);
    if (p == 0) {
        fprintf(stderr, "%s: %s: invalid argument\n", cmdname, argv[argc - 2]);
        usage(cmdname);
        return 0;
    }
    
    /*Get interval*/
    n = atoi(argv[argc - 1]);
    if (n == 0) {
        fprintf(stderr, "%s: %s: invalid argument\n", cmdname, argv[argc - 2]);
        usage(cmdname);
        return 0;
    }
    
    /*Set variables to default*/
    cmdname = basename(argv[0]);
    printpartition = 0;
    depths = calloc(p, sizeof(unsigned long));
    depthcount = 0;
    method = '0';
    completesearch = 0;
    problem = 0;
    cputimemax = ULLONG_MAX;
    bfilename = NULL;
    dfilename = NULL;
    ofilename = NULL;
    
    while ((optc = getopt(argc, argv, "ahpvb:d:o:m:t:")) != -1 && !problem) {
        /*Parse arguments*/
        switch (optc) {
                
            case 'b':
                asprintf(&bfilename, "%s", optarg);
                break;
                
            case 'o':
                asprintf(&ofilename, "%s", optarg);
                break;
                
            case 'h':
                usage(cmdname);
                break;
                
            case 'd':
                if (depthcount == p) {
                    problem = 1;
                    fprintf(stderr, "%s: too many depths\n", cmdname);
                    usage(cmdname);
                } else {
                    depths[depthcount] = atol(optarg);
                    depthcount++;
                }
                break;
                
            case 'p':
                printpartition = 1;
                break;
                
            case 'a':
                completesearch = 1;
                break;
                
            case 'm':
                method = *optarg;
                if (method != '0' && method != '1') {
                    problem = 1;
                    fprintf(stderr, "%s: -m %c: unknown method\n", cmdname, method);
                    usage(cmdname);
                }
                break;
                
            case 't':
                cputimemax = atol(optarg) * CLOCKS_PER_SEC * 60;
                break;
                
            default:
                problem = 1;
                fprintf(stderr, "%s: -%c: invalid option\n", cmdname, optc);
                usage(cmdname);
                break;
        }
    }
    
    if (argc - optind <= 0) {
        fprintf(stderr, "%s: %s: invalid argument\n", cmdname, argv[optind]);
        usage(cmdname);
        return 0;
    }
    
    if (!problem) {
        
        sort(depths, p);
        makeInitialPartition(p, n, &partitionstruc);
        d = schurNumberSymmetricBranch(n, &partitionstruc, depths[p-1], completesearch);
        
        if (printpartition) {
            printSymmetricPartition(p, n, partitionstruc.partition);
        }
        
        printf("%lu\n", d);
        
        partition_unalloc(&partitionstruc);
    }
    
    /*Nettoyage*/
    free(depths);
    
    return 0;
}
/*
int main() {
    printf("%lu\n", schurNumberSymmetricImposedPartition(532, 6));
    return 0;
}*/

