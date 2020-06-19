//
//  main.c
//  schurNumberPuncturedInterval2
//
//  Created by rubis on 28/04/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include <unistd.h>
#include <time.h>
#include <libgen.h>
#include <ctype.h>

#include "schurNumberWeakPuncturedInterval.h"
#include "../schurNumberGlobal/schurNumberIOAction.h"
//#include "../schurNumberGlobal/schurNumberThreads.h"

#ifdef schurNumberThreads_h

#define schur_number_launch(methodfunc, partitionstruc, action, nlimit) schur_number_threads_launch(partitionstruc, methodfunc, action, NULL, 0)

#else

typedef unsigned long (*schur_number_method_t)(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit);
#define schur_number_launch(methodfunc, partitionstruc, action, nlimit) do {\
    schur_number_save_thread_register((action)->save);\
    methodfunc(partitionstruc, action, nlimit);\
    } while(0)

#endif

void usage(char *cmdname) {
    fprintf(stderr,
            "usage: %s [-abcehtu] [-p (a|b|num)] [-m method] set_number n1 n2 partition\n"\
            "\t-a: Equivalent à -bucet\n"\
            "\t-b: Affichage du nombre de partitions de taille maximale\n"\
            "\t-u: Affichage du nombre de partitions non prolongeables\n"\
            "\t-c: Affichage du nombre total de partitions testées\n"\
            "\t-e: Affichage du nombre de partitions testées pour chaque thread\n"\
            "\t-m method: Selection d'une méthode\n"\
            "\t-f: Ensembles au format binaire (vecteur de 0 et 1)\n"\
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
    struct timespec time0, time1;
    
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
    
    if (argc2 <= 1) {
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
    
    // Initialisation de l'action
    schur_number_action_t action_s;
    unsigned long (*actionfunc)(mp_limb_t **, unsigned long, struct schurNumberIOAction *);
    size_t part_count_limit = 0;
    
    if (print_range) {
        if (isdigit(*print_range)) {
            part_count_limit = atol(print_range);
            actionfunc = schur_number_save_some_partition;
        } else {
            switch (*print_range) {
                case 'a':
                    actionfunc = schur_number_save_all_partition;
                    break;
                    
                default:
                    actionfunc = schur_number_default_action;
                    break;
            }
        }
        free(print_range);
    } else {
        actionfunc = schur_number_default_action;
    }
    
    schur_number_action_alloc(&action_s, p, actionfunc);
    action_s.count_limit = part_count_limit;
    
    // Allocation de la partition
    schur_number_partition_t partition_s;
    
    mp_size_t limballoc = PARTITION_2_LIMBSIZE(p);
    schur_number_partition_alloc(&partition_s, p);
    schur_number_partition_init(&partition_s, limballoc);
    
    if (0 < argc2) {
        // Récupérer la partition initiale depuis argv
        unsigned long p_init = argc2;  // Nombre d'ensembles passés en argument pour constituer la partition initiale
        if (p_init > p) {
            // Si il y a plus d'ensembles que p
            p_init = p;
        }
        
        unsigned long n_init = 0;    // Taille de la partition initiale
        
        for (unsigned long j = 0; j < p_init; j++) {
            
            unsigned long nmax = schur_number_get_set(partition_s.partition[j], partition_s.partitioninvert[j], limballoc, *arg_ptr, format);
            if (n_init < nmax) {
                n_init = nmax;
            }
            arg_ptr++;
        }
        
        partition_s.p = p_init;
        partition_s.n = n_init;
        partition_s.limbsize = (n_init >> 6) + 1;
    } else {
        fprintf(stderr, "schurNumber: not enough set for the constraint partition\n");
        usage(argv[0]);
    }
    
    // Sélection de la méthode
    schur_number_method_t methodfunc;
    switch (method) {
        case '1':
            methodfunc = schur_number_weak_punctured_interval;
            break;
            
        case 'w':
            methodfunc = schur_number_weak_punctured_interval;
            break;
            
        default:
            methodfunc = schur_number_weak_punctured_interval;
            break;
    }
    
    // Création de la sauvegarde temporaire
    schur_number_intermediate_save_t save_str;
    schur_number_save_alloc(&save_str, p, partition_s.n);
    action_s.save = &save_str;
    
    // Lancement du code
    clock_gettime(CLOCK_MONOTONIC, &time0);
    
    schur_number_launch(methodfunc, &partition_s, &action_s, 2 * partition_s.n);
    
    clock_gettime(CLOCK_MONOTONIC, &time1);
    
    // Detruction de la sauvegarde temporaire
    schur_number_save_dealloc(&save_str);
    
    // Affichage des résultats
    
    if (print_range) {
        schur_number_action_print_partitions(&action_s);
    }
    
    if (timeOption) {
        printf("Durée : %f secondes\n", difftime(time1.tv_sec, time0.tv_sec) + (double)(time1.tv_nsec - time0.tv_nsec) / 1000000000.0);
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
    schur_number_action_dealloc(&action_s);
    
    return 0;
}

