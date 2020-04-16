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
#include <ctype.h>
#include <string.h>

#include "schurNumberConstrainedBuild.h"
//#include "../schurNumberGlobal/schurNumberThreads.h"

#ifdef schurNumberThreads_h
    #define schurNumberLaunch(methodfunc, partitionstruc, action, constraint_partition, constraint_size, nlimit, load_balancing_opt) schurNumberThreadsLaunch(partitionstruc, methodfunc, action, constraint_partition, constraint_size, load_balancing_opt)
#else
    #define schurNumberLaunch(methodfunc, partitionstruc, action, constraint_partition, constraint_size, nlimit, load_balancing_opt) do {\
            schurNumberSaveThreadRegister((action)->save);\
            methodfunc(partitionstruc, action, nlimit, constraint_partition, constraint_size);\
        } while(0)
#endif

typedef unsigned long (*schur_number_method_t)(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit, mp_limb_t **constraint_partition, mp_size_t constraint_size);

void usage(char *cmdname) {
    fprintf(stderr,
            "usage: %s [-abcehtu] [-p (a|b|num)] [-m method] set_number constraint_partition [begin_partition]\n"\
            "\t-a: Equivalent à -bucet\n"\
            "\t-b: Affichage du nombre de partitions de taille maximale\n"\
            "\t-u: Affichage du nombre de partitions non prolongeables\n"\
            "\t-c: Affichage du nombre total de partitions testées\n"\
            "\t-e: Affichage du nombre de partitions testées pour chaque thread\n"\
            "\t-f: Ensembles au format binaire (vecteur de 0 et 1)\n"\
            "\t-m method: Selection d'une méthode\n"\
            "\t\t1: Recherche totalement exhaustive (par défaut)\n"\
            "\t-t: Affichage du temps d'exécution\n"\
            "\t-h: Affichage de la syntaxe\n"\
            "\t-p (a|b|num): Affichage de certaines partitions trouvées\n"\
            "\t\ta: toutes les partitions non prolongeables\n"\
            "\t\tb: toutes les partitions de taille maximale\n"\
            "\t\tnum: au plus num partitions de taille maximale pour chaque thread\n",
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
    int format = SCHUR_NUMBER_ELEMENT_FORMAT;
    clock_t time0;
    clock_t time1;
    
    // Analyse des options
    while ((c = getopt(argc, argv, "abcefhm:p:tu")) != -1) {
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
                
            case 'f':
                format = SCHUR_NUMBER_NUMERIC_FORMAT;
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
            actionfunc = schurNumberSaveSomePartition;
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
    
    mp_size_t constraint_size = (schurNumberGetPartition(p, constraint_partition, NULL, 1, arg_ptr, format) >> 6) + 1;
    
    arg_ptr += p;
    argc2 -= p;
    
    // Allocation de la partition
    schur_number_partition_t partition_s;
    
    mp_size_t limballoc = PARTITION_2_LIMBSIZE(p);
    schur_number_partition_alloc(&partition_s, p);
    schur_number_partition_init(&partition_s, limballoc);
    
    if (0 < argc2) {
        // Récupérer la partition initiale depuis argv
        unsigned long p_init = argc2;  // Nombre d'ensembles passés en argument pour constituer la partition initiale
        if (p_init > p) {
            p_init = p;
        }
        
        unsigned long n_init = 0;
        
        for (unsigned long j = 0; j < p_init; j++) {
            
            unsigned long nmax = schurNumberGetSet(partition_s.partition[j], partition_s.partitioninvert[j], limballoc, *arg_ptr, format);
            if (n_init < nmax) {
                n_init = nmax;
            }
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
    
    // Création de la sauvegarde temporaire
    schur_number_intermediate_save_t save_str;
    schurNumberSaveAlloc(&save_str, p, partition_s.n);
    action_s.save = &save_str;
    
    // Lancement du code
    time0 = clock();
    //schurNumberLaunch(methodfunc, &partition_s, &action_s, constraint_partition, constraint_size, mp_bits_per_limb * limballoc, threadPartitionNumberOption);
    time1 = clock();

    
    // Destruction de la sauvegarde temporaire
    schurNumberSaveDealloc(&save_str);
    
    // Affichage des résultats
    
    if (print_range) {
        schurNumberPrintPartitions(&action_s);
    }
    
    if (timeOption) {
        printf("Durée: %lu CPU_time = %f seconds \n", time1 - time0, ((double)(time1 - time0)) / CLOCKS_PER_SEC);
    }
    
    if (testedPartitionNumberOption) {
        printf("Nombre de partitions testées: %lu\n", action_s.iter_num);
    }
    
    if (unprolongeableSfPartitionNumbersOption) {
        printf("Nombre de partitions non prolongeables: %lu\n", action_s.count_all);
    }
    
    if (bestPartitionNumbersOption) {
        printf("Nombre de partitions de taille maximale: %lu\n", action_s.count_max);
    }
    
    printf("Taille maximale : %lu\n", action_s.nmax);
    
    // Nettoyage
    schur_number_partition_dealloc(&partition_s);
    schurNumberActionDealloc(&action_s);
    
    return 0;
}
