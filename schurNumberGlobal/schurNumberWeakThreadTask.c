//
//  schurNumberWeakThreadTask.c
//  SchurNumber
//
//  Created by rubis on 12/07/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberThreads.h"
#include "../SchurNumber/schurNumberMethods.h"

void schur_number_thread_task_weak3(schur_number_task_arg_t *arg) {
    /* La procédure se déroule en trois phases.
     */
    
    // Initialisation des variables
    unsigned long p = arg->p;
    schur_number_partition_queue_t *queue = arg->partition_queue;
    size_t *waiting_thread_count_ptr = arg->waiting_thread_count_ptr;
    schur_number_partition_t *partitionstruc = arg->partition_struc;
    
    schur_number_action_t *action = arg->action;
    schur_number_sum_action_alloc(action);
    
    pthread_mutex_t *mutex = arg->mutex;
    pthread_cond_t *cond = arg->cond;
    
    mp_size_t constraint_size = arg->constraint_size;
    mp_limb_t **constraint_partition = NULL;
    if (arg->constraint_partition) {
        constraint_partition = calloc(sizeof(mp_limb_t *), p);
        
        for (unsigned long j = 0; j < p; j++) {
            constraint_partition[j] = calloc(sizeof(mp_limb_t), constraint_size);
            mpn_copyd(constraint_partition[j], arg->constraint_partition[j], constraint_size);
        }
    }
    
    // Détermination de bornes inférieures sur WS(i), 1 ≤ i ≤ p
    unsigned long nbest;
    unsigned long *wsinf = calloc(p, sizeof(unsigned long));
    *wsinf = 2;
    unsigned long pow3 = 1;
    for (unsigned long i = 1; i < p; i++) {
        wsinf[i] = 7 * pow3 + i;
        pow3 *= 3;
    }
    nbest = wsinf[p - 1];
    
    // Enregistrement du thread
    schur_number_intermediate_save_t *save = action->save;
    schur_number_save_thread_register(save);
    
    char should_continue = 1;
    partition_queue_flag_t old_flag = INITIAL_FLAG;
    unsigned long nlimit = arg->nlimit;
    
    // Statistiques
    unsigned long iternum1 = 0;
    unsigned long iternum2 = 0;
    unsigned long iternum3 = 0;
    unsigned long iternum4 = 0;
    
    while (should_continue) {
        // Obtenir une partition de la file
        pthread_mutex_lock(mutex);
        partition_queue_flag_t flag = schur_partition_queue_get_partition(queue, partitionstruc);
        while (flag == NULL_FLAG) {
            // Plus de partitions dans le file
            if (*waiting_thread_count_ptr < NUM_THREADS - 1) {
                // Attendre le signal pour réessayer de prendre une partition ou s'arrêter
                (*waiting_thread_count_ptr)++;
                pthread_cond_wait(cond, mutex);
                (*waiting_thread_count_ptr)--;
                // Ré-essayer d'acquérir une partition
                flag = schur_partition_queue_get_partition(queue, partitionstruc);
            } else {
                // Sonner la fin de la traque
                schur_partition_queue_add_partitionarray_nocopy(queue, NULL, NULL, 0, 0, STOP_FLAG);
                flag = STOP_FLAG;
                // Envoyer un signal aux threads en attente
                pthread_cond_broadcast(cond);
                break;
            }
        }
        pthread_mutex_unlock(mutex);
        
        old_flag = flag;
        
        schur_number_action_func_t action_func = action->func;
        size_t count = action->count;
        unsigned long iternum0 = action->iter_num;
        
        unsigned long i = p - 1;
        while (i > 0) {
            if (partitionstruc->n >= wsinf[--i] << 1) {
                break;
            }
            if (partitionstruc->n > (3 * wsinf[i]) >> 1) {
                // Sauter directement au placement de 2 * wsinf[i] + 1
                partitionstruc->n = wsinf[i] << 1;
                break;
            }
        }
        
        switch (flag) {
            case INITIAL_FLAG:
            {   // Effectuer une recherche limitée à nbest * 2/3
                //nlimit = ((2 * nbest) / 3) + 1;
                //nlimit = (nbest / 3) + (nbest / 6) + 4;
                nlimit = 12;
                action->nthreshold = nlimit;
                //action->func = schur_number_save_distinct_sum_partition;
                action->func = schur_number_save_threshold_partition;
                //action->func = schur_number_save_distinct_begin_partition;
                
                // Effectuer une recherche illimitée
                unsigned long n = arg->func(partitionstruc, action, nlimit, constraint_partition, constraint_size);
                iternum1 += action->iter_num - iternum0;
                //printf("Partitions : %3lu\tItérations : %3lu\n", action->count - count, action->iter_num - iternum0);
                
                if (n == partitionstruc->n) {
                    // La partition n'est pas prolongeable
                    action_func(partitionstruc->partition, n, action);
                }
                
                if (n > nbest) {
                    nbest = n;
                }
                
                action->func = action_func;
                action->count_all -= (action->count - count);
                // Ajouter ses propres partitions
                mp_size_t *limbsize_buffer;
                mp_limb_t *partition_buffer;
                count = schur_number_action_buffer_detach(action, &limbsize_buffer, &partition_buffer);
                pthread_mutex_lock(mutex);
                schur_partition_queue_add_partitionarray_nocopy(queue, partition_buffer, limbsize_buffer, count, 16, INTERMEDIATE_FLAG);
                // Envoyer un signal aux threads en attente
                pthread_cond_broadcast(cond);
                pthread_mutex_unlock(mutex);
                break;
            }
                
            case INTERMEDIATE_FLAG:
            {   // Effectuer une recherche limitée à nbest * 5/6
                //nlimit = ((2 * nbest) / 3) + ((2 * (nbest - ((2 * nbest) / 3))) / 3);
                nlimit = 18;
                action->nthreshold = nlimit;
                //action->func = schur_number_save_distinct_sum_partition;
                //action->func = schur_number_save_threshold_partition;
                action->func = schur_number_save_distinct_begin_partition;
                
                // Effectuer une recherche illimitée
                unsigned long n = arg->func(partitionstruc, action, nlimit, constraint_partition, constraint_size);
                iternum2 += action->iter_num - iternum0;
                
                if (n == partitionstruc->n) {
                    // La partition n'est pas prolongeable
                    action->func(partitionstruc->partition, n, action);
                }
                
                if (n > nbest) {
                    nbest = n;
                }
                
                action->func = action_func;
                action->count_all -= (action->count - count);
                // Ajouter ses propres partitions
                mp_size_t *limbsize_buffer;
                mp_limb_t *partition_buffer;
                count = schur_number_action_buffer_detach(action, &limbsize_buffer, &partition_buffer);
                pthread_mutex_lock(mutex);
                schur_partition_queue_add_partitionarray_nocopy(queue, partition_buffer, limbsize_buffer, count, nlimit, FINAL_FLAG);
                // Envoyer un signal aux threads en attente
                pthread_cond_broadcast(cond);
                pthread_mutex_unlock(mutex);
                break;
            }
                
            case FINAL_FLAG:
            {
                // Effectuer une recherche illimitée
                nlimit = arg->nlimit;
                unsigned long n = arg->func(partitionstruc, action, nlimit, constraint_partition, constraint_size);
                iternum4 += action->iter_num - iternum0;
                
                if (n == partitionstruc->n) {
                    // La partition n'est pas prolongeable
                    action_func(partitionstruc->partition, n, action);
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
    
    pthread_mutex_lock(mutex);
    printf("Thread %p:\n\tPhase 1: %lu test\n\tPhase 2: %lu test\n\tPhase 3: %lu test\n\tPhase 4: %lu test\n", pthread_self(), iternum1, iternum2, iternum3, iternum4);
    pthread_mutex_unlock(mutex);
    
}

void schur_number_thread_task_weak2(schur_number_task_arg_t *arg) {
    /* La procédure se déroule en trois phases.
     */
    
    // Initialisation des variables
    unsigned long p = arg->p;
    schur_number_partition_queue_t *queue = arg->partition_queue;
    size_t *waiting_thread_count_ptr = arg->waiting_thread_count_ptr;
    schur_number_partition_t *partitionstruc = arg->partition_struc;
    
    schur_number_action_t *action = arg->action;
    schur_number_sum_action_alloc(action);
    
    pthread_mutex_t *mutex = arg->mutex;
    pthread_cond_t *cond = arg->cond;
    
    mp_size_t constraint_size = arg->constraint_size;
    mp_limb_t **constraint_partition = NULL;
    if (arg->constraint_partition) {
        constraint_partition = calloc(sizeof(mp_limb_t *), p);
        
        for (unsigned long j = 0; j < p; j++) {
            constraint_partition[j] = calloc(sizeof(mp_limb_t), constraint_size);
            mpn_copyd(constraint_partition[j], arg->constraint_partition[j], constraint_size);
        }
    }
    
    unsigned long nbest;
    switch (p) {
        case 1:
            nbest = 2;
            break;
            
        default:
        {
            unsigned long pow3 = 1;
            unsigned long pdiff = p - 2;
            while (pdiff > 0) {
                pow3 *= 3;
                pdiff--;
            }
            nbest = 7 * pow3 + p - 1;
        }
            break;
    }
    
    
    // Enregistrement du thread
    schur_number_intermediate_save_t *save = action->save;
    schur_number_save_thread_register(save);
    
    char should_continue = 1;
    partition_queue_flag_t old_flag = INITIAL_FLAG;
    unsigned long nlimit = arg->nlimit;
    
    while (should_continue) {
        // Obtenir une partition de la file
        pthread_mutex_lock(mutex);
        partition_queue_flag_t flag = schur_partition_queue_get_partition(queue, partitionstruc);
        while (flag == NULL_FLAG) {
            // Plus de partitions dans le file
            
            if (old_flag != FINAL_FLAG && action->count) {
                // Placer ses propres partitions dans la file au besoin
                mp_size_t *limbsize_buffer;
                mp_limb_t *partition_buffer;
                size_t count = schur_number_action_buffer_detach(action, &limbsize_buffer, &partition_buffer);
                partition_queue_flag_t new_flag = STOP_FLAG;
                schur_partition_queue_add_partitionarray_nocopy(queue, partition_buffer, limbsize_buffer, count, nlimit, new_flag);
                // Envoyer un signal aux threads en attente
                pthread_cond_broadcast(cond);
            } else if (*waiting_thread_count_ptr < NUM_THREADS - 1) {
                // Attendre le signal pour réessayer de prendre une partition ou s'arrêter
                (*waiting_thread_count_ptr)++;
                pthread_cond_wait(cond, mutex);
                (*waiting_thread_count_ptr)--;
            } else {
                // Sonner la fin de la traque
                schur_partition_queue_add_partitionarray_nocopy(queue, NULL, NULL, 0, 0, STOP_FLAG);
                flag = STOP_FLAG;
                // Envoyer un signal aux threads en attente
                pthread_cond_broadcast(cond);
                break;
            }
            // Ré-essayer d'acquérir une partition
            flag = schur_partition_queue_get_partition(queue, partitionstruc);
            
            //            if (*waiting_thread_count_ptr < NUM_THREADS - 1) {
            //                // Attendre le signal pour réessayer de prendre une partition ou s'arrêter
            //                (*waiting_thread_count_ptr)++;
            //                pthread_cond_wait(cond, mutex);
            //                (*waiting_thread_count_ptr)--;
            //                // Ré-essayer d'acquérir une partition
            //                flag = schur_partition_queue_get_partition(queue, partitionstruc);
            //            } else {
            //                // Sonner la fin de la traque
            //                schur_partition_queue_add_partitionarray_nocopy(queue, NULL, NULL, 0, 0, STOP_FLAG);
            //                flag = STOP_FLAG;
            //                // Envoyer un signal aux threads en attente
            //                pthread_cond_broadcast(cond);
            //                break;
            //            }
        }
        pthread_mutex_unlock(mutex);
        
        if (old_flag != flag && action->count) {
            // Ajouter ses propres partitions à la file
            mp_size_t *limbsize_buffer;
            mp_limb_t *partition_buffer;
            size_t count = schur_number_action_buffer_detach(action, &limbsize_buffer, &partition_buffer);
            partition_queue_flag_t new_flag = STOP_FLAG;
            pthread_mutex_lock(mutex);
            schur_partition_queue_add_partitionarray_nocopy(queue, partition_buffer, limbsize_buffer, count, nlimit, new_flag);
            // Envoyer un signal aux threads en attente
            pthread_cond_broadcast(cond);
            pthread_mutex_unlock(mutex);
        }
        old_flag = flag;
        
        schur_number_action_func_t action_func = action->func;
        size_t count = action->count;
        switch (flag) {
            case INITIAL_FLAG:
            {   // Effectuer une recherche limitée à nbest * 2/3
                nlimit = (nbest / 2);
                action->nthreshold = nlimit;
                //action->func = schur_number_save_distinct_sum_partition;
                action->func = schur_number_save_threshold_partition;
                
                // Effectuer une recherche illimitée
                unsigned long n = arg->func(partitionstruc, action, nlimit, constraint_partition, constraint_size);
                
                if (n == partitionstruc->n) {
                    // La partition n'est pas prolongeable
                    action_func(partitionstruc->partition, n, action);
                }
                
                if (n > nbest) {
                    nbest = n;
                }
                
                action->func = action_func;
                action->count_all -= (action->count - count);
                break;
            }
                
            case FINAL_FLAG:
            {
                // Effectuer une recherche illimitée
                nlimit = arg->nlimit;
                unsigned long n = arg->func(partitionstruc, action, nlimit, constraint_partition, constraint_size);
                
                if (n == partitionstruc->n) {
                    // La partition n'est pas prolongeable
                    action_func(partitionstruc->partition, n, action);
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

void schur_number_thread_task_weak(schur_number_task_arg_t *arg) {
    /* La procédure se déroule en trois phases.
     */
    
    // Initialisation des variables
    unsigned long p = arg->p;
    schur_number_partition_queue_t *queue = arg->partition_queue;
    size_t *waiting_thread_count_ptr = arg->waiting_thread_count_ptr;
    schur_number_partition_t *partitionstruc = arg->partition_struc;
    
    schur_number_action_t *action = arg->action;
    schur_number_sum_action_alloc(action);
    
    pthread_mutex_t *mutex = arg->mutex;
    pthread_cond_t *cond = arg->cond;
    
    mp_size_t constraint_size = arg->constraint_size;
    mp_limb_t **constraint_partition = NULL;
    if (arg->constraint_partition) {
        constraint_partition = calloc(sizeof(mp_limb_t *), p);
        
        for (unsigned long j = 0; j < p; j++) {
            constraint_partition[j] = calloc(sizeof(mp_limb_t), constraint_size);
            mpn_copyd(constraint_partition[j], arg->constraint_partition[j], constraint_size);
        }
    }
    
    unsigned long nbest;
    switch (p) {
        case 1:
            nbest = 2;
            break;
            
        default:
        {
            unsigned long pow3 = 1;
            unsigned long pdiff = p - 2;
            while (pdiff > 0) {
                pow3 *= 3;
                pdiff--;
            }
            nbest = 7 * pow3 + p - 1;
        }
            break;
    }
    
    
    // Enregistrement du thread
    schur_number_intermediate_save_t *save = action->save;
    schur_number_save_thread_register(save);
    
    char should_continue = 1;
    partition_queue_flag_t old_flag = INITIAL_FLAG;
    unsigned long nlimit = arg->nlimit;
    
    // Statistiques
    unsigned long iternum1 = 0;
    unsigned long iternum2 = 0;
    unsigned long iternum3 = 0;
    unsigned long iternum4 = 0;
    
    while (should_continue) {
        // Obtenir une partition de la file
        pthread_mutex_lock(mutex);
        partition_queue_flag_t flag = schur_partition_queue_get_partition(queue, partitionstruc);
        while (flag == NULL_FLAG) {
            // Plus de partitions dans le file
            
            if (*waiting_thread_count_ptr < NUM_THREADS - 1) {
                // Attendre le signal pour réessayer de prendre une partition ou s'arrêter
                (*waiting_thread_count_ptr)++;
                pthread_cond_wait(cond, mutex);
                (*waiting_thread_count_ptr)--;
                // Ré-essayer d'acquérir une partition
                flag = schur_partition_queue_get_partition(queue, partitionstruc);
            } else {
                // Sonner la fin de la traque
                schur_partition_queue_add_partitionarray_nocopy(queue, NULL, NULL, 0, 0, STOP_FLAG);
                flag = STOP_FLAG;
                // Envoyer un signal aux threads en attente
                pthread_cond_broadcast(cond);
                break;
            }
        }
        pthread_mutex_unlock(mutex);
        old_flag = flag;
        
        schur_number_action_func_t action_func = action->func;
        size_t count = action->count;
        unsigned long iternum0 = action->iter_num;
        switch (flag) {
            case INITIAL_FLAG:
            {   // Effectuer une recherche limitée à nbest * 2/3
                //nlimit = ((2 * nbest) / 3) + 1;
                //nlimit = (nbest / 3) + (nbest / 6) + 4;
                nlimit = 12;
                action->nthreshold = nlimit;
                action->func = schur_number_save_distinct_sum_partition;
                //action->func = schur_number_save_threshold_partition;
                //action->func = schur_number_save_distinct_begin_partition;
                
                // Effectuer une recherche illimitée
                unsigned long n = arg->func(partitionstruc, action, nlimit, constraint_partition, constraint_size);
                iternum1 += action->iter_num - iternum0;
                //printf("Partitions : %3lu\tItérations : %3lu\n", action->count - count, action->iter_num - iternum0);
                
                if (n == partitionstruc->n) {
                    // La partition n'est pas prolongeable
                    action_func(partitionstruc->partition, n, action);
                }
                
                if (n > nbest) {
                    nbest = n;
                }
                
                action->func = action_func;
                action->count_all -= (action->count - count);
                // Ajouter ses propres partitions
                mp_size_t *limbsize_buffer;
                mp_limb_t *partition_buffer;
                count = schur_number_action_buffer_detach(action, &limbsize_buffer, &partition_buffer);
                pthread_mutex_lock(mutex);
                schur_partition_queue_add_partitionarray_nocopy(queue, partition_buffer, limbsize_buffer, count, 16, INTERMEDIATE_FLAG);
                // Envoyer un signal aux threads en attente
                pthread_cond_broadcast(cond);
                pthread_mutex_unlock(mutex);
                break;
            }
                
            case INTERMEDIATE_FLAG:
            {   // Effectuer une recherche limitée à nbest * 5/6
                //nlimit = ((2 * nbest) / 3) + ((2 * (nbest - ((2 * nbest) / 3))) / 3);
                nlimit = 17;
                action->nthreshold = nlimit;
                //action->func = schur_number_save_distinct_sum_partition;
                //action->func = schur_number_save_threshold_partition;
                action->func = schur_number_save_distinct_begin_partition;
                
                // Effectuer une recherche illimitée
                unsigned long n = arg->func(partitionstruc, action, nlimit, constraint_partition, constraint_size);
                iternum2 += action->iter_num - iternum0;
                
                if (n == partitionstruc->n) {
                    // La partition n'est pas prolongeable
                    action->func(partitionstruc->partition, n, action);
                }
                
                if (n > nbest) {
                    nbest = n;
                }
                
                action->func = action_func;
                action->count_all -= (action->count - count);
                // Ajouter ses propres partitions
                mp_size_t *limbsize_buffer;
                mp_limb_t *partition_buffer;
                count = schur_number_action_buffer_detach(action, &limbsize_buffer, &partition_buffer);
                pthread_mutex_lock(mutex);
                /*printf("Nombre : %lu\n", count);
                for (int i = 0; i < count; i++) {
                    mp_limb_t *part[3];
                    part[0] = partition_buffer + i * p * (*limbsize_buffer);
                    part[1] = part[0] + *limbsize_buffer;
                    part[2] = part[1] + *limbsize_buffer;
                    schur_number_dprint_partition(1, p, 23, part);
                }*/
                schur_partition_queue_add_partitionarray_nocopy(queue, partition_buffer, limbsize_buffer, count, nlimit+1, FINAL_FLAG);
                // Envoyer un signal aux threads en attente
                pthread_cond_broadcast(cond);
                pthread_mutex_unlock(mutex);
                break;
            }
                
                /*case FLAG_3:
                 {   // Effectuer une recherche limitée à nbest * 5/6
                 nlimit = ((2 * nbest) / 3) + ((2 * (nbest - ((2 * nbest) / 3))) / 3);
                 action->nthreshold = nlimit;
                 //action->func = schur_number_save_distinct_sum_partition;
                 //action->func = schur_number_save_threshold_partition;
                 action->func = schur_number_save_distinct_begin_partition;
                 
                 // Effectuer une recherche illimitée
                 unsigned long n = arg->func(partitionstruc, action, nlimit, constraint_partition, constraint_size);
                 iternum3 += action->iter_num - iternum0;
                 
                 if (n == partitionstruc->n) {
                 // La partition n'est pas prolongeable
                 action->func(partitionstruc->partition, n, action);
                 }
                 
                 if (n > nbest) {
                 nbest = n;
                 }
                 
                 action->func = action_func;
                 action->count_all -= (action->count - count);
                 // Ajouter ses propres partitions
                 mp_size_t *limbsize_buffer;
                 mp_limb_t *partition_buffer;
                 count = schur_number_action_buffer_detach(action, &limbsize_buffer, &partition_buffer);
                 pthread_mutex_lock(mutex);
                 schur_partition_queue_add_partitionarray_nocopy(queue, partition_buffer, limbsize_buffer, count, nlimit, FINAL_FLAG);
                 // Envoyer un signal aux threads en attente
                 pthread_cond_broadcast(cond);
                 pthread_mutex_unlock(mutex);
                 break;
                 }*/
                
            case FINAL_FLAG:
            {
                // Effectuer une recherche illimitée
                nlimit = arg->nlimit;
                unsigned long n = arg->func(partitionstruc, action, nlimit, constraint_partition, constraint_size);
                iternum4 += action->iter_num - iternum0;
                
                if (n == partitionstruc->n) {
                    // La partition n'est pas prolongeable
                    action_func(partitionstruc->partition, n, action);
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
    
    pthread_mutex_lock(mutex);
    printf("Thread %p:\n\tPhase 1: %lu test\n\tPhase 2: %lu test\n\tPhase 3: %lu test\n\tPhase 4: %lu test\n", pthread_self(), iternum1, iternum2, iternum3, iternum4);
    pthread_mutex_unlock(mutex);
    
}
