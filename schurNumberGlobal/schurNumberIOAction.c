//
//  schurNumberIOAction.c
//  SchurNumber
//
//  Created by rubis on 25/03/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberIOAction.h"

#if defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__APPLE__)

#include <sys/sysctl.h>
#define set_size_limit(size_limit_ptr) do {\
        int64_t memsize;\
        size_t len = sizeof(memsize);\
        sysctlbyname("hw.memsize", &memsize, &len, NULL, 0);\
        *(size_limit_ptr) = memsize / 16;\
    } while(0)

#elif defined(__linux__)

#include <sys/sysinfo.h>
#define set_size_limit(size_limit_ptr) do {\
        struct sysinfo info_s;\
        sysinfo(&info_s);\
        *(size_limit_ptr) = info_s.totalram / 16;\
    } while(0)

#else

#define set_size_limit(size_limit_ptr) do {\
        *(size_limit_ptr) = 16777216;\
    } while(0)

#endif

void schur_number_action_alloc(schur_number_action_t *action, unsigned long p, unsigned long (*func)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action)) {
    action->p = p;
    action->count = 0;
    action->nbest = 0;
    action->iter_num = 0;
    action->count_all = 0;
    action->count_max = 0;
    action->count_limit = 0;
    set_size_limit(&(action->size_limit));
    
    action->limbsize_buffer = NULL;
    action->limbsize_size = 0;
    action->limbsize_stream = open_memstream(&(action->limbsize_buffer), &(action->limbsize_size));
    
    action->partition_buffer = NULL;
    action->partition_size = 0;
    action->partition_stream = open_memstream(&(action->partition_buffer), &(action->partition_size));
    
    action->action_flag = SCHUR_NUMBER_DEFAULT;
    
    if ((func == schur_number_save_distinct_sum_partition) || (func == schur_number_save_distinct_restrictedsum_partition)) {
        action->sum_partition_stream = open_memstream(&(action->sum_partition_buffer), &(action->sum_partition_size));
        action->sorted_index_sum_partition_stream = open_memstream(&(action->sorted_index_sum_partition_buffer), &(action->sorted_index_sum_partition_size));
        action->work = calloc(PARTITION_2_LIMBSIZE(p), 4 * sizeof(mp_limb_t));
        action->limbsize = PARTITION_2_LIMBSIZE(p);
    } else {
        action->sum_partition_buffer = NULL;
        action->sum_partition_size = 0;
        action->sorted_index_sum_partition_buffer = NULL;
        action->sorted_index_sum_partition_size = 0;
        action->work = NULL;
    }
    
    action->func = func;
    
    action->save = NULL;
    
    action->count_gathered_actions = 0;
    action->gathered_actions = NULL;
}

void schur_number_action_dealloc(schur_number_action_t *action) {
    
    fclose(action->limbsize_stream);
    fclose(action->partition_stream);
    
    free(action->limbsize_buffer);
    free(action->partition_buffer);
    
    if (action->sum_partition_buffer) {
        fclose(action->sum_partition_stream);
        free(action->sum_partition_buffer);
        fclose(action->sorted_index_sum_partition_stream);
        free(action->sorted_index_sum_partition_buffer);
        free(action->work);
    }
    
    schur_number_action_t **actions = action->gathered_actions;
    if (actions) {
        unsigned long count = action->count_gathered_actions;
        for (unsigned long i = 0; i < count; i++) {
            schur_number_action_dealloc(actions[i]);
        }
        free(actions);
    }
}

void schur_number_action_gather_copy(schur_number_action_t *action_r, schur_number_action_t **actions, size_t n_actions) {
    /*Réunit les n_actions actions dans l'unique action_r, en copiant les multiples tampons dans celui de action_r.*/
    unsigned long nmax = action_r->nbest;
    unsigned long iter_num = action_r->iter_num;
    size_t count_all = action_r->count_all;
    size_t count_max = action_r->count_max;
    size_t count = action_r->count;
    FILE * limbsize_stream = action_r->limbsize_stream;
    FILE * partition_stream = action_r->partition_stream;
    schur_number_action_func_t func = action_r->func;
    
    for (unsigned long i = 0; i < n_actions; i++) {
        schur_number_action_t *action_s = actions[i];
        
        if (action_s->nbest > nmax) {
            nmax = action_s->nbest;
            count_max = 0;
        }
        
        if (action_s->nbest == nmax) {
            count_max += action_s->count_max;
        }
        
        count_all += action_s->count_all;
        iter_num += action_s->iter_num;
        
        if (func == schur_number_save_all_partition) {
            fflush(action_s->limbsize_stream);
            fflush(action_s->partition_stream);
            
            fwrite(action_s->limbsize_buffer, 1, action_s->limbsize_size, limbsize_stream);
            fwrite(action_s->partition_buffer, 1, action_s->partition_size, partition_stream);
            
            count += action_s->count;
        }
    }
    
    if (func == schur_number_save_best_partition || func == schur_number_save_some_partition) {
        
        if (action_r->nbest < nmax) {
            rewind(limbsize_stream);
            rewind(partition_stream);
            count = 0;
        }
        
        for (unsigned long i = 0; i < n_actions; i++) {
            schur_number_action_t *action_s = actions[i];
            
            if (action_s->nbest == nmax) {
                fflush(action_s->limbsize_stream);
                fflush(action_s->partition_stream);
                
                fwrite(action_s->limbsize_buffer, 1, action_s->limbsize_size, limbsize_stream);
                fwrite(action_s->partition_buffer, 1, action_s->partition_size, partition_stream);
                
                count += action_s->count;
            }
        }
    }
    
    action_r->nbest = nmax;
    action_r->count_all = count_all;
    action_r->count_max = count_max;
    action_r->iter_num = iter_num;
    action_r->count = count;
}

void schur_number_action_gather_nocopy(schur_number_action_t *action_r, schur_number_action_t **actions, size_t n_actions) {
    /*Réunit les n_actions actions_s dans l'unique action_r, sans copier les multiples tampons dans celui de action_r.
     Les actions_s sont automatiquement libérées.*/
    unsigned long nmax = action_r->nbest;
    schur_number_action_func_t func = action_r->func;
    
    size_t old_n_actions = action_r->count_gathered_actions; // Nombre de tampons déjà présents au sein de action_r
    
    // Déterminer le nbest global et le nombre d'actions l'atteignant
    for (unsigned long i = 0; i < n_actions; i++) {
        schur_number_action_t *action = actions[i];
        
        if (action->nbest > nmax) {
            nmax = action->nbest;
        }
    }
    
    schur_number_action_t **gathered_actions = realloc(action_r->gathered_actions, (old_n_actions + n_actions) * sizeof(schur_number_action_t*));
    
    if (gathered_actions) {
        
        if (func == schur_number_save_best_partition || func == schur_number_save_some_partition) {
            
            if (action_r->nbest < nmax) {
                // Remplacer le nbest de action_r
                schur_number_save_best_partition(NULL, nmax, action_r);
                
                for (size_t i = 0; i < old_n_actions; i++) {
                    schur_number_save_best_partition(NULL, nmax, actions[i]);
                }
            }
            
            for (size_t i = 0; i < n_actions; i++) {
                schur_number_action_t *action = actions[i];
                
                if (action->nbest < nmax) {
                    schur_number_save_best_partition(NULL, nmax, actions[i]);
                }
                
                gathered_actions[i + old_n_actions] = action;
            }
            
        } else {
            // Les cas schurNumberSaveDefault et schurNumberSaveAll
            for (size_t i = 0; i < n_actions; i++) {
                schur_number_action_t *action = actions[i];
                gathered_actions[i + old_n_actions] = action;
            }
        }
        
        action_r->gathered_actions = gathered_actions;
        action_r->count_gathered_actions += n_actions;
        
    } else {
        schur_number_action_gather_copy(action_r, actions, n_actions);
        
        for (unsigned long i = 0; i < n_actions; i++) {
            schur_number_action_dealloc(actions[i]);
        }
    }
}

size_t schurNumberPrintPartitionBuffer(unsigned long p, char *limbsize_buffer, char *partition_buffer, size_t count) {
    /*Affiche dans stdout les count partitions à p ensembles contenues dans partition_buffer, et de taille en limbes précisée dans limbsize_buffer.*/
    
    mp_size_t *limbsize_ptr = limbsize_buffer;
    mp_limb_t *set_ptr = partition_buffer;
    mp_limb_t **partition = calloc(sizeof(mp_limb_t *), p);
    
    size_t k;
    
    for (k = 0; k < count; k++) {
        
        mp_size_t limbsize = *limbsize_ptr;
        
        for (unsigned long j = 0; j < p; j++) {
            partition[j] = set_ptr;
            set_ptr += limbsize;
        }
        schurNumberPrintPartition(p, limbsize * mp_bits_per_limb - 1, partition);
        
        limbsize_ptr ++;
    }
    
    free(partition);
    
    return k;
}

size_t schur_number_action_print_partitions(schur_number_action_t *action) {
    /*Affiche toutes les partitions présentes dans action.*/
    unsigned long p = action->p;
    
    fflush(action->limbsize_stream);
    fflush(action->partition_stream);
    
    size_t total_count = schurNumberPrintPartitionBuffer(p, action->limbsize_buffer, action->partition_buffer, action->count);
    
    schur_number_action_t **actions = action->gathered_actions;
    for (size_t k = 0; k < action->count_gathered_actions; k++) {
        total_count += schur_number_action_print_partitions(actions[k]);
    }
    
    return total_count;
}

unsigned long schur_number_action_total_iterations(const schur_number_action_t *action) {
    /* Renvoie le nombre total d'itérations effectuées. */
    unsigned long iternum = action->iter_num;
    
    schur_number_action_t **actions = action->gathered_actions;
    for (unsigned long i = 0; i < action->count_gathered_actions; i++) {
        iternum += schur_number_action_total_iterations(actions[i]);
    }
    
    return iternum;
}

unsigned long schur_number_action_total_count_all(const schur_number_action_t *action) {
    /* Renvoie le nombre de partitions trouvées parmi toutes les actions. */
    unsigned long count = action->count_all;
    
    schur_number_action_t **actions = action->gathered_actions;
    for (unsigned long i = 0; i < action->count_gathered_actions; i++) {
        count += schur_number_action_total_count_all(actions[i]);
    }
    
    return count;
}

unsigned long schur_number_action_total_count_max(const schur_number_action_t *action) {
    /* Renvoie le nombre de partitions de taille maximale trouvées parmi toutes les actions. */
    unsigned long count = action->count_max;
    
    schur_number_action_t **actions = action->gathered_actions;
    for (unsigned long i = 0; i < action->count_gathered_actions; i++) {
        count += schur_number_action_total_count_max(actions[i]);
    }
    
    return count;
}

unsigned long schur_number_action_total_Nmax(const schur_number_action_t *action) {
    /* Renvoie la plus grande taille de partitions trouvée parmi toutes les actions. */
    unsigned long nmax = action->nbest;
    
    schur_number_action_t **actions = action->gathered_actions;
    for (unsigned long i = 0; i < action->count_gathered_actions; i++) {
        unsigned long n = schur_number_action_total_Nmax(actions[i]);
        if (n > nmax) {
            nmax = n;
        }
    }
    
    return nmax;
}
