//
//  schurNumberThreads.h
//  SchurNumber
//
//  Created by rubis on 12/02/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#ifndef schurNumberThreads_h
#define schurNumberThreads_h

#include <pthread.h>
#include <unistd.h>
#include <unistd.h>
#include "schurNumberIO.h"
#include "schurNumberPartitionStruc.h"

#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif

typedef unsigned long (*schur_number_method_t)(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit, mp_limb_t **constraint_partition, mp_size_t constraint_size);

struct schur_number_task_arg_struc {
    unsigned long p;            // Nombre de huches allouées
    unsigned long *nbest_ptr;   // Pointeur vers la taille de la meilleure partition trouvée
    
    size_t count;               // Nombre de partitions initiales
    size_t *current_index_ptr;  // Indice de la partition initiale à sélectionner
    schur_number_partition_t *partitionstruc_array; // Tableau des partitions initiales
    schur_number_action_t *action;                  // Actions à effectuer
    
    mp_size_t constraint_size;
    mp_limb_t **constraint_partition;
    
    schur_number_method_t func; // Fonction à exécuter
    
    pthread_mutex_t *mutex;     // Mutex protégeant le tableau
};

typedef struct schur_number_task_arg_struc schur_number_task_arg_t;

unsigned long schurNumberThreadsLaunch(schur_number_partition_t *partitionstruc, schur_number_method_t methodfunc, schur_number_action_t *action, mp_limb_t **constraint_partition, mp_size_t constraint_size, char load_balancing_stat);

#endif /* schurNumberThreads_h */
