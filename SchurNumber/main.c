//
//  main.c
//  SchurNumber
//
//  Created by rubis on 23/02/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include <unistd.h>
#include <time.h>
#include <libgen.h>

/*
#include "schurNumberIterative.c"
#include "schurNumberImprovedIterative.c"
#include "shur_number_bound.c"
#include "schur_number_unstack_bound.c"
#include "schur_number_unstack_memory.c"
//#include "schur_number_unorganized_memory.c"
#include "schurNumberWithPrint.c"
#include "schurNumberWeakImprovedIterative.c"*/

#include "schurNumberMethods.h"
#include "../schurNumberGlobal/schurNumberIO.h"
#include "../schurNumberGlobal/schurNumberThreads.h"

#ifdef schurNumberThreads_h

#define schurNumberLaunch(methodfunc, partitionstruc, action, nlimit, load_balancing_opt) schurNumberThreadsLaunch(partitionstruc, methodfunc, action, NULL, 0, load_balancing_opt)

#else

typedef unsigned long (*schur_number_method_t)(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit);
#define schurNumberLaunch(methodfunc, partitionstruc, action, nlimit, load_balancing_opt) methodfunc(partitionstruc, action, nlimit)

#endif

void usage(char *cmdname) {
    fprintf(stderr,
           "usage: %s [-abcehtu] [-p (1|a|b)] [-m method] set_number [begin_partition]\n"\
            "\t-a: Equivalent to -bucet\n"\
            "\t-b: Print number of sum-free partitions of [1, S(p)]\n"\
            "\t-u: Print number of sum-free partitions which are not prolongeable\n"\
            "\t-c: Print total number of tested partitions\n"\
            "\t-e: Print number of tested partitions for each thread\n"\
            "\t-m method: Select a method\n"\
            "\t-t: Print execution's time\n"\
            "\t-h: Print usage message\n"\
            "\t-p (1|a|b): Show 1 best partition, all not prolongeable partitions or only best found partitions\n",
           basename(cmdname));
}

#ifndef test

int main(int argc, const char * argv[]) {
    char c;
    char timeOption = 0;
    char bestSfPartitionNumbersOption = 0;
    char unprolongeableSfPartitionNumbersOption = 0;
    char testedPartitionNumberOption = 0;
    char threadPartitionNumberOption = 0;
    char method = 0;
    char print_range = 0;
    unsigned long *schurNumbers;
    mpz_t iternum;
    clock_t time0;
    clock_t time1;
    
    // Analyse des options
    while ((c = getopt(argc, argv, "abcehm:p:tu")) != -1) {
        switch (c) {
            case 'a':
                timeOption = 1;
                bestSfPartitionNumbersOption = 1;
                unprolongeableSfPartitionNumbersOption = 1;
                testedPartitionNumberOption = 1;
                threadPartitionNumberOption = 1;
                break;
                
            case 'b':
                bestSfPartitionNumbersOption = 1;
                break;
            
            case 'u':
                unprolongeableSfPartitionNumbersOption = 1;
                break;
            
            case 'c':
                testedPartitionNumberOption = 1;
                break;
                
            case 'e':
                threadPartitionNumberOption = 1;
                break;
                
            case 'h':
                usage(argv[0]);
                return 0;
                
            case 't':
                timeOption = 1;
                break;
            
            case 'm':
                method = optarg[0];
                break;
            
            case 'p':
                print_range = optarg[0];
                break;
                
            default:
                fprintf(stderr, "schurNumber: -%c: invalid option\n", c);
                usage(argv[0]);
                return 0;
        }
    }
    if (argc - optind <= 0) {
        fprintf(stderr, "schurNumber: %s: invalid argument\n", argv[optind]);
        usage(argv[0]);
        return 0;
    }
    unsigned long p = atol(argv[optind]);
    if (p == 0) {
        fprintf(stderr, "schurNumber: %s: invalid argument\n", argv[optind]);
        usage(argv[0]);
        return 0;
    }
    
    schurNumbers = calloc(p, sizeof(unsigned long));
    //partitionNumbers = calloc(p, sizeof(unsigned long));
    mpz_init_set_ui(iternum, 0);
    
    // Initialisation de l'action
    schur_number_action_t action_s;
    void (*actionfunc)(mp_limb_t **, unsigned long, struct schurNumberIOAction *);
    switch (print_range) {
        case '1':
            actionfunc = schurNumberSaveOnePartition;
            break;
            
        case 'a':
            actionfunc = schurNumberSaveAllPartition;
            break;
            
        case 'b':
            actionfunc = schurNumberSaveBestPartition;
            break;
            
        default:
            actionfunc = schurNumberDefaultAction;
            break;
    }
    schurNumberActionAlloc(&action_s, p, actionfunc);
    
    // Allocation de la partition
    schur_number_partition_t partition_s;
    
    mp_size_t limballoc = ((4 << p)>>6) + 1;
    schur_number_partition_alloc(&partition_s, limballoc, p);
    
    if (optind + 1 < argc) {
        // Récupérer la partition initiale depuis argv
        unsigned long p_init = argc - (optind + 1);  // Nombre d'ensembles passés en argument pour constituer la partition initiale
        if (p_init > p) {
            p_init = p;
        }
        char **set_str_ptr = &(argv[optind + 1]);
        
        unsigned long n_init = 0;
        
        for (unsigned long j = 0; j < p_init; j++) {
            
            unsigned long nmax = schurNumberGetSetMaximum(*set_str_ptr);
            if (n_init < nmax) {
                n_init = nmax;
            }
            
            schurNumberGetSet(*set_str_ptr, partition_s.partition[j], partition_s.partitioninvert[j], limballoc);
            set_str_ptr++;
        }
        
        partition_s.p = p_init;
        partition_s.n = n_init;
        partition_s.limbsize = (n_init >> 6) + 1;
    } else {
        // Initialisation à {{1}}
        partition_s.p = 1;
        partition_s.n = 1;
        ADD_POINT(*(partition_s.partition), 1);
        ADD_POINT(partition_s.partitioninvert[0], limballoc * mp_bits_per_limb - 1);
    }
    
    // Sélection de la méthode
    schur_number_method_t methodfunc;
    /*switch (method) {
        case '1':
            //schurNumberIterative1(p, schurNumbers, iternum);
            methodfunc = schurNumberIterative22;
            break;
        
        case '2':
            methodfunc = schurNumberIterative2;
            break;
        
        case '3':
            methodfunc = schurNumberIterWithUnstack;
            //schurNumberIterative21(p, schurNumbers, iternum);
            break;
        
        case '4':
            //schurNumberIterBound(p, schurNumbers, iternum);
            break;
            
        case '5':
            methodfunc = schurNumberIterBound2;
            break;
            
        case '6':
            methodfunc = schurNumberIterBoundUnstack;
            break;
            
        case '7':
            methodfunc = schurNumberIterBoundUnstack2;
            break;
            
        case '8':
            methodfunc = schurNumberIterUnstackMemory;
            break;
        
        case 'p':
            methodfunc = schurNumberIterWithPrint;
            break;
            
        case 'w':
            methodfunc = schurNumberWeakExhaustive;
            break;
            
        default:
            //methodfunc = schurNumberIterative;
            methodfunc = schurNumberExhaustive;
            break;
    }*/
    switch (method) {
        case '1':
            methodfunc = schurNumberExhaustive;
            break;
            
        case '4':
            methodfunc = schurNumberWeakBranchBound;
            break;
            
        case 'w':
            methodfunc = schurNumberWeakExhaustive;
            break;
            
        default:
            methodfunc = schurNumberExhaustive;
            break;
    }
    
    // Lancement du code
    time0 = clock();
    schurNumberLaunch(methodfunc, &partition_s, &action_s, mp_bits_per_limb * limballoc, threadPartitionNumberOption);
    time1 = clock();
    
    // Affichage des résultats
    
    /*for (i=0; i<p; i++) {
        printf("Schur Number S(%lu) = %lu\n", i+1, schurNumbers[i]);
    }*/
    free(schurNumbers);
    
    if (print_range) {
        schurNumberPrintPartitions(&action_s);
    }
    
    if (timeOption) {
        printf("Time: %lu CPU_time = %f seconds \n", time1 - time0, ((double)(time1 - time0)) / CLOCKS_PER_SEC);
    }
    
    if (testedPartitionNumberOption) {
        printf("Number of tested partitions: %lu\n", action_s.iter_num);
    }
    
    if (unprolongeableSfPartitionNumbersOption) {
        printf("Number of unprolongeable partitions: %lu\n", action_s.count_all);
    }
    
    if (bestSfPartitionNumbersOption) {
        printf("Number of maximal sized partitions: %lu\n", action_s.count_max);
    }
    
    // Nettoyage
    schur_number_partition_dealloc(&partition_s);
    schurNumberActionDealloc(&action_s);
    
    return 0;
}

#else

int main(int argc, const char * argv[]) {
    mp_limb_t limbs[2] = {0, 0};
    
    while (*limbs < 15) {
        (*limbs)++;
        unsigned long pos = mpn_rscan1(limbs, 128);
        printf("%lu : %lu\n", *limbs, pos);
    }
    
    for (unsigned long i=5; i<64; i++) {
        (*limbs) <<= 1;
        unsigned long pos = mpn_rscan1(limbs, 128);
        printf("%lu : %lu\n", *limbs, pos);
    }
    
    limbs[1] = 5;
    
    for (unsigned long i=3; i<64; i++) {
        limbs[1] <<= 1;
        unsigned long pos = mpn_rscan1(limbs, 128);
        printf("%lu : %lu\n", *limbs, pos);
    }
    
    return 0;
}

#endif
