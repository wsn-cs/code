//
//  schurNumberThreads.c
//  SchurNumber
//
//  Created by rubis on 12/02/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberThreads.h"
//#include "../SchurNumber/schurNumberMethods.h"

size_t schur_number_initiate_partition_queue(schur_number_task_arg_t *arg) {
    /* Cette fonction remplit partition_queue avec au moins 2 * NUM_THREADS partitions, obtenus en cherchant successivement des partitions sans-sommes de taille n+4. */
    
    // Initialisation des variables
    unsigned long p = arg->p;
    schur_number_partition_queue_t *partition_queue = arg->partition_queue;
    schur_number_partition_t *partitionstruc = arg->partition_struc;
    schur_number_action_t *action = arg->action;
    
    unsigned long n = partitionstruc->n;
    size_t count = partition_queue->total_count;
    
    mp_size_t constraint_size = arg->constraint_size;
    mp_limb_t **constraint_partition = NULL;
    if (arg->constraint_partition) {
        constraint_partition = calloc(sizeof(mp_limb_t *), p);
        
        for (unsigned long j = 0; j < p; j++) {
            constraint_partition[j] = calloc(sizeof(mp_limb_t), constraint_size);
            mpn_copyd(constraint_partition[j], arg->constraint_partition[j], constraint_size);
        }
    }
    
    // Changer l'action
    schur_number_action_func_t action_func = action->func;
    schur_number_intermediate_save_t *save = action->save;
    action->func = schur_number_save_best_partition;
    action->save = NULL;
    
    char has_improved = 1;
    
    while (count < 2 * NUM_THREADS && has_improved) {
        unsigned long nbest = n;
        // Recherche des partitions
        while (schur_partition_queue_get_partition(partition_queue, partitionstruc)) {
            unsigned long nfound = arg->func(partitionstruc, action, n + 4, constraint_partition, constraint_size);
            if (nfound > nbest) {
                nbest = nfound;
            }
        }
        has_improved = (nbest > n);
        n = nbest;
        
        //schur_number_action_print_partitions(action);
        
        mp_size_t *limbsize_buffer;
        mp_limb_t *partition_buffer;
        count = schur_number_action_buffer_detach(action, &limbsize_buffer, &partition_buffer);
        count = schur_partition_queue_add_partitionarray_nocopy(partition_queue, partition_buffer, limbsize_buffer, count, 1);
    }
    // Rétablir l'action
    action->func = action_func;
    action->save = save;
    
    if (arg->constraint_partition) {
        unsigned long p = arg->p;
        
        for (unsigned long j = 0; j < p; j++) {
            free(constraint_partition[j]);
        }
        
        free(constraint_partition);
    }
    
    return count;
}

void schur_number_thread_task_weak(schur_number_task_arg_t *arg) {
    /* La procédure se déroule en trois phases.
     */
    
    // Initialisation des variables
    schur_number_partition_queue_t *queue = arg->partition_queue;
    size_t *waiting_thread_count_ptr = arg->waiting_thread_count_ptr;
    schur_number_partition_t *partitionstruc = arg->partition_struc;
    schur_number_action_t *action = arg->action;
    pthread_mutex_t *mutex = arg->mutex;
    pthread_cond_t *cond = arg->cond;
    
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
    
    
    // Enregistrement du thread
    schur_number_intermediate_save_t *save = action->save;
    schur_number_save_thread_register(save);
    
    char should_continue = 1;
    char old_flag = 1;
    
    while (should_continue) {
        // Obtenir une partition de la file
        pthread_mutex_lock(mutex);
        char flag = schur_partition_queue_get_partition(queue, partitionstruc);
        while (!flag) {
            // Plus de partitions dans le file
            if (*waiting_thread_count_ptr < NUM_THREADS - 1) {
                // Attendre le signal pour réessayer de prendre une partition ou s'arrêter
                (*waiting_thread_count_ptr)++;
                pthread_cond_wait(cond, mutex);
                (*waiting_thread_count_ptr)--;
                flag = schur_partition_queue_get_partition(queue, partitionstruc);
            } else {
                // Placer ses propres partitions dans la file au besoin
                if (old_flag == 1 && action->count) {
                    mp_size_t *limbsize_buffer;
                    mp_limb_t *partition_buffer;
                    size_t count = schur_number_action_buffer_detach(action, &limbsize_buffer, &partition_buffer);
                    schur_partition_queue_add_partitionarray_nocopy(queue, partition_buffer, limbsize_buffer, count, 2);
                } else {
                    schur_partition_queue_add_partitionarray_nocopy(queue, NULL, NULL, 0, -1);
                }
                // Envoyer un signal aux threads en attente
                pthread_cond_broadcast(cond);
            }
        }
        pthread_mutex_unlock(mutex);
        
        old_flag = flag;
        
        unsigned long nlimit = arg->nlimit;
        switch (flag) {
            case 1:
                // Effectuer une recherche limitée à nbest * 5/6
                nlimit = 5 * (nbest / 6);
                
            case 2:
            {   // Effectuer une recherche illimitée
                unsigned long n = arg->func(partitionstruc, action, nlimit, constraint_partition, constraint_size);
                
                if (n == partitionstruc->n) {
                    // La partition n'est pas prolongeable
                    action->func(partitionstruc->partition, n, action);
                }
                
                if (n > nbest) {
                    nbest = n;
                }
            }
                break;
                
            default:
                // Drapeau d'arrêt
                should_continue = 0;
                break;
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

void schur_number_thread_task(schur_number_task_arg_t *arg) {
    /* Lance la recherche sur une des partitions de la file. */
    
    // Initialisation des variables
    schur_number_partition_queue_t *queue = arg->partition_queue;
    schur_number_partition_t *partitionstruc = arg->partition_struc;
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
    
    // Enregistrement du thread
    schur_number_intermediate_save_t *save = action->save;
    schur_number_save_thread_register(save);
    
    char should_continue = 1;
    
    while (should_continue) {
        // Obtenir une partition de la file
        pthread_mutex_lock(mutex);
        char flag = schur_partition_queue_get_partition(queue, partitionstruc);
        /*while (!(flag = schur_partition_queue_get_partition(queue, partitionstruc))) {
            // Plus de partitions dans le file
            if (*waiting_thread_count_ptr < NUM_THREADS - 1) {
                // Attendre le signal pour réessayer de prendre une partition ou s'arrêter
                (*waiting_thread_count_ptr)++;
                pthread_cond_wait(cond, mutex);
                (*waiting_thread_count_ptr)--;
                flag = schur_partition_queue_get_partition(queue, partitionstruc);
            } else {
                // Placer ses propres partitions dans la file au besoin
                if (old_flag == 1 && action->count) {
                    mp_size_t *limbsize_buffer;
                    mp_limb_t *partition_buffer;
                    size_t count = schur_number_action_buffer_detach(action, &limbsize_buffer, &partition_buffer);
                    schur_partition_queue_add_partitionarray_nocopy(queue, partition_buffer, limbsize_buffer, count, 1);
                } else {
                    schur_partition_queue_add_partitionarray_nocopy(queue, NULL, NULL, 0, -1);
                }
                // Envoyer un signal aux threads en attente
                pthread_cond_broadcast(cond);
            }
        }*/
        pthread_mutex_unlock(mutex);
        
        unsigned long nlimit = arg->nlimit;
        switch (flag) {
            case 1:
            {   // Effectuer une recherche illimitée
                unsigned long n = arg->func(partitionstruc, action, nlimit, constraint_partition, constraint_size);
                
                if (n == partitionstruc->n) {
                    // La partition n'est pas prolongeable
                    action->func(partitionstruc->partition, n, action);
                }
                
                if (n > nbest) {
                    nbest = n;
                }
            }
                break;
                
            default:
                // Drapeau d'arrêt
                should_continue = 0;
                break;
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

unsigned long schur_number_threads_launch(schur_number_partition_t *partitionstruc, schur_number_method_t methodfunc, schur_number_action_t *action, unsigned long nlimit, mp_limb_t **constraint_partition, mp_size_t constraint_size) {
    /*Cette fonction initialise plusieurs threads. La répartition se fait en spécifiant différentes partitions de départ pour n = 4.*/
    
    pthread_t threads[NUM_THREADS - 1];
    schur_number_partition_t *partitions[NUM_THREADS - 1];
    schur_number_action_t *actions[NUM_THREADS - 1];
    schur_number_task_arg_t arg_s_array[NUM_THREADS];
    schur_number_partition_queue_t partition_queue_struc;
    size_t waiting_thread_count = 0;
    unsigned long nbest = 0;
    
    // Récupération des variables globales
    unsigned long p = partitionstruc->pmax;
    mp_size_t limballoc = partitionstruc->limballoc;
    
    // Création de la mutex associée
    pthread_mutex_t mutex_s;
    pthread_mutex_init(&mutex_s, NULL);
    
    // Création de la condition associée à la vacuité de la file
    pthread_cond_t cond_s;
    pthread_cond_init(&cond_s, NULL);
    
    // Création de l'argument à appeler par le thread principal
    schur_number_task_arg_t *main_arg = &(arg_s_array[NUM_THREADS - 1]);
    main_arg->p = p;
    main_arg->nbest_ptr = &nbest;
    main_arg->partition_queue = &partition_queue_struc;
    main_arg->partition_struc = partitionstruc;
    main_arg->action = action;
    main_arg->nlimit = nlimit;
    main_arg->constraint_size = constraint_size;
    main_arg->constraint_partition = constraint_partition;
    main_arg->func = methodfunc;
    main_arg->waiting_thread_count_ptr = &waiting_thread_count;
    main_arg->mutex = &mutex_s;
    main_arg->cond = &cond_s;
    
    // Création des partitions initiales
    schur_number_partition_queue_init(&partition_queue_struc);
    schur_partition_queue_add_partition_copy(&partition_queue_struc, partitionstruc, 1);
    size_t count = schur_number_initiate_partition_queue(main_arg); // Nombre de partitions dans partitionstruc_array
    
    // Répercussion dans la sauvegarde
    schur_number_intermediate_save_t *save = action->save;
    schur_number_save_partition_pool_register(save, count, partition_queue_struc.n);
    
    for (unsigned int i = 0; i < NUM_THREADS - 1; i++) {
        // Création de l'argument à appeler
        schur_number_task_arg_t *arg = &(arg_s_array[i]);
        arg->p = p;
        arg->nbest_ptr = &nbest;
        arg->partition_queue = &partition_queue_struc;
        
        partitions[i] = calloc(1, sizeof(schur_number_partition_t));
        schur_number_partition_alloc(partitions[i], p);
        schur_number_partition_init(partitions[i], limballoc);
        arg->partition_struc = partitions[i];
        
        actions[i] = calloc(1, sizeof(schur_number_action_t));
        schur_number_action_alloc(actions[i], p, action->func);
        actions[i]->count_limit = action->count_limit;
        actions[i]->nbest = action->nbest;
        arg->action = actions[i];
        
        arg->nlimit = nlimit;
        arg->constraint_size = constraint_size;
        arg->constraint_partition = constraint_partition;
        arg->func = methodfunc;
        arg->waiting_thread_count_ptr = &waiting_thread_count;
        arg->mutex = &mutex_s;
        arg->cond = &cond_s;
        actions[i]->save = save;
        
        // Création du thread
        pthread_create(&threads[i], NULL, schur_number_thread_task, arg);
    }
    
    // Exécution de la fonction
    schur_number_thread_task(main_arg);
    
    // Regroupement
    for (unsigned int i = 0; i < NUM_THREADS - 1; i++) {
        pthread_join(threads[i], NULL);
        schur_number_partition_dealloc(partitions[i]);
    }
    
    schur_number_action_gather_nocopy(action, actions, NUM_THREADS - 1);
    
    // Nettoyage
    schur_number_partition_queue_dealloc(&partition_queue_struc);
    pthread_mutex_destroy(&mutex_s);
    pthread_cond_destroy(&cond_s);
    
    return 0;
}
