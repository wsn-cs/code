//
//  main.c
//  schurNumberRecursive
//
//  Created by rubis on 01/12/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include <unistd.h>
#include <time.h>
#include <libgen.h>

#include "schurNumberConstrainedBuild.h"
#include "../schurNumberGlobal/schurNumberIO.h"
#include "../schurNumberGlobal/schurNumberThreads.h"

#ifdef schurNumberThreads_h
#define schurNumberLaunch(methodfunc, partitionstruc, action, constraint_partition, constraint_size, nlimit, load_balancing_opt) schurNumberThreadsLaunch(partitionstruc, methodfunc, action, constraint_partition, constraint_size, load_balancing_opt)
#else
#define schurNumberLaunch(methodfunc, partitionstruc, action, constraint_partition, constraint_size, nlimit, load_balancing_opt) methodfunc(partitionstruc, action, nlimit, constraint_partition, constraint_size)
#endif

typedef unsigned long (*schur_number_method_t)(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit, mp_limb_t **constraint_partition, mp_size_t constraint_size);

void usage(char *cmdname) {
    fprintf(stderr,
            "usage: %s [-abcehtu] [-p (1|a|b)] [-m method] set_number constraint_partition [begin_partition]\n"\
            "\t-a: Equivalent to -bucet\n"\
            "\t-b: Print number of partitions of maximal size\n"\
            "\t-u: Print number of partitions which are not prolongeable\n"\
            "\t-c: Print total number of tested partitions\n"\
            "\t-e: Print number of tested partitions for each thread\n"\
            "\t-m method: Select a method\n"\
            "\t-t: Print execution's time\n"\
            "\t-h: Print usage message\n"\
            "\t-p (1|a|b): Show 1 best partition, all not prolongeable partitions or only best found partitions\n",
            basename(cmdname));
}

int main(int argc, const char * argv[]) {
    char c;
    char timeOption = 0;
    char bestPartitionNumbersOption = 0;
    char unprolongeableSfPartitionNumbersOption = 0;
    char testedPartitionNumberOption = 0;
    char threadPartitionNumberOption = 0;
    char method = 0;
    char *print_range = NULL;
    clock_t time0;
    clock_t time1;
    
    // Analyse des options
    while ((c = getopt(argc, argv, "abcehm:p:tu")) != -1) {
        switch (c) {
            case 'a':
                timeOption = 1;
                bestPartitionNumbersOption = 1;
                unprolongeableSfPartitionNumbersOption = 1;
                testedPartitionNumberOption = 1;
                threadPartitionNumberOption = 1;
                break;
                
            case 'b':
                bestPartitionNumbersOption = 1;
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
                asprintf(&print_range, "%s", optarg);
                break;
                
            default:
                fprintf(stderr, "schurNumber: -%c: invalid option\n", c);
                usage(argv[0]);
                return 0;
        }
    }
    char **arg_ptr = argv + optind;
    int argc2 = argc - optind;
    
    if (argc2 <= 0) {
        fprintf(stderr, "schurNumber: %s: invalid argument\n", *arg_ptr);
        usage(argv[0]);
        return 0;
    }
    
    unsigned long p = atol(*arg_ptr);
    if (p == 0) {
        fprintf(stderr, "schurNumber: %s: invalid argument\n", *arg_ptr);
        usage(argv[0]);
        return 0;
    }
    argc2--;
    arg_ptr++;
    if (argc2 < p) {
        fprintf(stderr, "schurNumber: not enough set for the constraint partition\n");
        usage(argv[0]);
        return 0;
    }
    
    // Initialisation de l'action
    schur_number_action_t action_s;
    void (*actionfunc)(mp_limb_t **, unsigned long, struct schurNumberIOAction *);
    size_t part_count_limit = 0;
    
    if (print_range) {
        if (isdigit(*print_range)) {
            part_count_limit = atol(print_range);
        } else {
            switch (*print_range) {
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
        }
        free(print_range);
    } else {
        actionfunc = schurNumberDefaultAction;
    }
    
    schurNumberActionAlloc(&action_s, p, actionfunc);
    action_s.count_limit = part_count_limit;
    
    // Allocation de la partition de contrainte
    mp_limb_t **constraint_partition = calloc(sizeof(mp_limb_t *), p);
    
    mp_size_t constraint_size = (schurNumberGetPartition(p, arg_ptr, constraint_partition, NULL, 1) >> 6) + 1;
    
    arg_ptr += p;
    argc2 -= p;
    
    // Allocation de la partition
    schur_number_partition_t partition_s;
    
    mp_size_t limballoc = ((4 << p)>>6) + 1;
    schur_number_partition_alloc(&partition_s, limballoc, p);
    
    if (0 < argc2) {
        // Récupérer la partition initiale depuis argv
        unsigned long p_init = argc2;  // Nombre d'ensembles passés en argument pour constituer la partition initiale
        if (p_init > p) {
            p_init = p;
        }
        
        unsigned long n_init = 0;
        
        for (unsigned long j = 0; j < p_init; j++) {
            
            unsigned long nmax = schurNumberGetSetMaximum(*arg_ptr);
            if (n_init < nmax) {
                n_init = nmax;
            }
            
            schurNumberGetSet(*arg_ptr, partition_s.partition[j], partition_s.partitioninvert[j], limballoc);
            arg_ptr++;
        }
        
        partition_s.p = p_init;
        partition_s.n = n_init;
        partition_s.limbsize = (n_init >> 6) + 1;
    } else {
        // Initialisation à {{0}}
        //partition_s.p = 1;
        //partition_s.n = 0;
        ADD_POINT(*(partition_s.partition), 0);
        //ADD_POINT(partition_s.partitioninvert[0], limballoc * mp_bits_per_limb - 1);
    }
    
    // Sélection de la méthode
    schur_number_method_t methodfunc;
    switch (method) {
        case '1':
            methodfunc = schurNumberConstrainedBuild;
            break;
            
        default:
            methodfunc = schurNumberConstrainedBuild;
            break;
    }
    
    // Lancement du code
    time0 = clock();
    schurNumberLaunch(methodfunc, &partition_s, &action_s, constraint_partition, constraint_size, mp_bits_per_limb * limballoc, threadPartitionNumberOption);
    time1 = clock();
    
    // Affichage des résultats
    
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
    
    if (bestPartitionNumbersOption) {
        printf("Number of maximal sized partitions: %lu\n", action_s.count_max);
    }
    
    printf("Taille maximale : %lu\n", action_s.nmax);
    
    // Nettoyage
    schur_number_partition_dealloc(&partition_s);
    schurNumberActionDealloc(&action_s);
    
    return 0;
}
