//
//  schurNumberThreads.c
//  SchurNumber
//
//  Created by rubis on 12/02/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberThreads.h"

size_t schurNumberPartitionPool(schur_number_partition_t *beginpartitionstruc, schur_number_partition_t **partitionstruc_array_ptr, schur_number_method_t methodfunc, mp_limb_t **constraint_partition, mp_size_t constraint_size) {
    /*Cette fonction crée une liste de partitions initiales à utiliser comme thread poll.
     Les partitions sont énumérées selon la méthode jusqu'à n = 4.
     Elle renvoie le nombre de partitions initiales, et alloue la mémoire associée à partitionstruc_array.*/
    unsigned long p = beginpartitionstruc->pmax;
    unsigned long n = beginpartitionstruc->n;
    mp_size_t limballoc = beginpartitionstruc->limballoc;
    unsigned long nalloc = limballoc * GMP_NUMB_BITS;
    
    schur_number_action_t action_s;
    #ifdef schurNumberConstrainedBuild_h
    schurNumberActionAlloc(&action_s, p, schurNumberConstrainedBuild);
    #else
    schurNumberActionAlloc(&action_s, p, schurNumberSaveBestPartition);
    #endif
    
    size_t count = 1;
    schur_number_partition_t *work_partitionstruc_array = beginpartitionstruc;
    char has_to_be_freed = 0;
    char has_improved = 1;
    
    while (count <= 2 * NUM_THREADS && has_improved) {
        // Recherche des partitions
        for (unsigned long i = 0; i < count; i++) {
            methodfunc(&(work_partitionstruc_array[i]), &action_s, n + 4, constraint_partition, constraint_size);
        }
        fflush(action_s.limbsize_stream);
        fflush(action_s.partition_stream);
        
        if (has_to_be_freed) {
            for (unsigned long i = 0; i < count; i++) {
                schur_number_partition_dealloc(&(work_partitionstruc_array[i]));
            }
            free(work_partitionstruc_array);
        }
        has_improved = (action_s.nmax == n+4);
        
        // Mise à jour des variables
        count = action_s.count;
        n = action_s.nmax;
        mp_size_t limbsize = (n >> 6) + 1;
        work_partitionstruc_array = calloc(sizeof(schur_number_partition_t), count);
        has_to_be_freed = 1;
        
        // Remplissage du tableau contenant les partitions
        mp_limb_t *set = action_s.partition_buffer;
        
        for (unsigned long i = 0; i < count; i++) {
            schur_number_partition_t *partitionstruc = &(work_partitionstruc_array[i]);
            schur_number_partition_alloc(partitionstruc, p);
            schur_number_partition_init(partitionstruc, limballoc);
            
            partitionstruc->n = n;
            partitionstruc->limbsize = limbsize;
            
            unsigned long j = 0;
            unsigned long last_non_zero_index = 0;
            while (j < p) {
                if (!mpn_zero_p(set, limbsize)) {
                    last_non_zero_index = j;
                }
                for (unsigned long k = 1; k <= n; k++) {
                    if (GET_POINT(set, k)) {
                        ADD_POINT(partitionstruc->partition[j], k);
                        ADD_POINT(partitionstruc->partitioninvert[j], nalloc - k);
                    }
                }
                set += limbsize;
                j++;
            }
            
            partitionstruc->p = last_non_zero_index + 1;
            set += (p - j) * limbsize;
        }
    }
    
    *partitionstruc_array_ptr = work_partitionstruc_array;
    
    schurNumberActionPrintPartitions(&action_s);
    schurNumberActionDealloc(&action_s);
    
    return count;
}

void schurNumberPartitionPoolDealloc(schur_number_partition_t **partitionstruc_array_ptr, size_t count) {
    /*Libère la mémoire associée à partitionstruc_array_ptr, qui contient count partitionstruc à libérer.*/
    for (unsigned long i = 0; i < count; i++) {
        schur_number_partition_dealloc(&(*partitionstruc_array_ptr)[i]);
    }
    free(*partitionstruc_array_ptr);
}

void schurNumberThreadTask(schur_number_task_arg_t *arg) {
    /*Lance la fonction sur une des partitions dans partitionstruc_array.*/
    
    // Initialisation des variables
    size_t count = arg->count;
    size_t *current_index_ptr = arg->current_index_ptr;
    schur_number_partition_t *partitionstruc_array = arg->partitionstruc_array;
    schur_number_action_t *action = arg->action;
    pthread_mutex_t *mutex = arg->mutex;
    
    mp_size_t constraint_size = arg->constraint_size;
    mp_limb_t **constraint_partition = NULL;
    if (arg->constraint_partition) {
        unsigned long p = arg->p;
        constraint_partition = calloc(sizeof(mp_limb_t *), p);
        
        for (unsigned long j = 0; j < p; j++) {
            constraint_partition[j] = calloc(sizeof(mp_limb_t), constraint_size);
            mpn_copyd(constraint_partition[j], arg->constraint_partition[j], constraint_size);
        }
    }
    
    unsigned long nbest = 0;
    unsigned long i;
    
    // Enregistrement du thread
    schur_number_intermediate_save_t *save = action->save;
    schurNumberSaveThreadRegister(save);

    // Lire l'indice courant
    pthread_mutex_lock(mutex);
    i = *current_index_ptr;
    *current_index_ptr = i + 1;
    pthread_mutex_unlock(mutex);
    
    while (i < count) {
        // Lance la procédure
        schur_number_partition_t *partitionstruc = &(partitionstruc_array[i]);
        
        schurNumberSaveNewExplorationRegister(save);
        
        unsigned long n;
        unsigned long nalloc = partitionstruc->limballoc * GMP_NUMB_BITS;
        
        n = arg->func(partitionstruc, action, nalloc, constraint_partition, constraint_size);
        
        if (n > nbest) {
            nbest = n;
        }
        
        // Lire l'indice courant
        pthread_mutex_lock(mutex);
        i = *current_index_ptr;
        *current_index_ptr = i + 1;
        unsigned long nbest_global = *(arg->nbest_ptr);
        if (nbest_global < nbest) {
            nbest_global = nbest;
            *(arg->nbest_ptr) = nbest_global;
        }
        pthread_mutex_unlock(mutex);
        
        if (nbest < nbest_global) {
            schur_number_action_t *action = arg->action;
            action->func(NULL, nbest_global, action);
        }
    }
    
    if (arg->constraint_partition) {
        unsigned long p = arg->p;
        
        for (unsigned long j = 0; j < p; j++) {
            free(constraint_partition[j]);
        }
        
        free(constraint_partition);
    }
}

unsigned long schurNumberThreadsLaunch(schur_number_partition_t *partitionstruc, schur_number_method_t methodfunc, schur_number_action_t *action, mp_limb_t **constraint_partition, mp_size_t constraint_size) {
    /*Cette fonction initialise plusieurs threads. La répartition se fait en spécifiant différentes partitions de départ pour n = 4.*/
    
    pthread_t threads[NUM_THREADS - 1];
    schur_number_action_t *actions[NUM_THREADS - 1];
    schur_number_task_arg_t arg_s_array[NUM_THREADS];
    schur_number_partition_t *partitionstruc_array = NULL;
    
    // Récupération des variables globales
    unsigned long p = partitionstruc->pmax;
    
    // Création des partitions initiales.
    size_t count = schurNumberPartitionPool(partitionstruc, &partitionstruc_array, methodfunc, constraint_partition, constraint_size); // Nombre de partitions dans partitionstruc_array
    size_t current_index = 0;                                                                   // Indice de la partition à traiter
    unsigned long nbest = 0;
    
    // Répercussion dans la sauvegarde
    schur_number_intermediate_save_t *save = action->save;
    schurNumberSavePartitionPoolRegister(save, count, partitionstruc_array->n);
    
    // Création de la mutex associée
    pthread_mutex_t mutex_s;
    pthread_mutex_init(&mutex_s, NULL);
    
    for (unsigned int i = 0; i < NUM_THREADS - 1; i++) {
        // Création de l'argument à appeler
        schur_number_task_arg_t *arg = &(arg_s_array[i]);
        arg->p = p;
        arg->nbest_ptr = &nbest;
        arg->count = count;
        arg->current_index_ptr = &current_index;
        arg->partitionstruc_array = partitionstruc_array;
        actions[i] = calloc(1, sizeof(schur_number_action_t));
        schurNumberActionAlloc(actions[i], p, action->func);
        actions[i]->count_limit = action->count_limit;
        arg->action = actions[i];
        arg->constraint_size = constraint_size;
        arg->constraint_partition = constraint_partition;
        arg->func = methodfunc;
        arg->mutex = &mutex_s;
        actions[i]->save = save;
        
        // Création du thread
        pthread_create(&threads[i], NULL, schurNumberThreadTask, arg);
    }
    
    // Création de l'argument à appeler
    schur_number_task_arg_t *arg = &(arg_s_array[NUM_THREADS - 1]);
    arg->p = p;
    arg->nbest_ptr = &nbest;
    arg->count = count;
    arg->current_index_ptr = &current_index;
    arg->partitionstruc_array = partitionstruc_array;
    arg->action = action;
    arg->constraint_size = constraint_size;
    arg->constraint_partition = constraint_partition;
    arg->func = methodfunc;
    arg->mutex = &mutex_s;
    
    // Exécution de la fonction
    schurNumberThreadTask(arg);
    
    // Regroupement
    for (unsigned int i = 0; i < NUM_THREADS - 1; i++) {
        pthread_join(threads[i], NULL);
    }
    
    //schurNumberActionGatherCopy(action, actions, NUM_THREADS - 1);
    schurNumberActionGatherNoCopy(action, actions, NUM_THREADS - 1);
    
//    for (unsigned int i = 0; i < NUM_THREADS - 1; i++) {
//        schurNumberActionDealloc(actions[i]);
//        free(actions[i]);
//    }
    
    // Nettoyage
    schurNumberPartitionPoolDealloc(&partitionstruc_array, count);
    pthread_mutex_destroy(&mutex_s);
    
    return 0;
}
