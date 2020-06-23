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
#include "schurNumberIOAction.h"
#include "schurNumberPartitionStruc.h"
#include "schurNumberPartitionQueue.h"

#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif

typedef unsigned long (*schur_number_method_t)(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit, mp_limb_t **constraint_partition, mp_size_t constraint_size);

struct schur_number_task_arg_struc {
    unsigned long p;            // Nombre de huches allouées
    unsigned long *nbest_ptr;   // Pointeur vers la taille de la meilleure partition trouvée
    
    schur_number_partition_queue_t *partition_queue;// File des partitions initiales
    schur_number_partition_t *partition_struc;      // Structure de partition
    schur_number_action_t *action;                  // Actions à effectuer
    
    size_t *waiting_thread_count_ptr;               // Variable pointant vers le nombre de threads en attente
    
    unsigned long nlimit;       // Limite sur les partitions cherchées
    
    mp_size_t constraint_size;
    mp_limb_t **constraint_partition;
    
    schur_number_method_t func; // Fonction à exécuter
    
    pthread_mutex_t *mutex;     // Mutex protégeant le tableau
    pthread_cond_t *cond;       // Condition associée à la vacuité de la file
};

typedef struct schur_number_task_arg_struc schur_number_task_arg_t;


unsigned long schur_number_threads_launch(schur_number_partition_t *partitionstruc, schur_number_method_t methodfunc, schur_number_action_t *action, unsigned long nlimit, mp_limb_t **constraint_partition, mp_size_t constraint_size);

#endif /* schurNumberThreads_h */
