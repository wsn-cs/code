//
//  schurNumberIOAction.c
//  SchurNumber
//
//  Created by rubis on 25/03/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberIOAction.h"


void schurNumberActionAlloc(schur_number_action_t *action, unsigned long p, void (*func)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action)) {
    action->p = p;
    action->count = 0;
    action->nmax = 0;
    action->iter_num = 0;
    action->count_all = 0;
    action->count_max = 0;
    action->count_limit = 0;
    
    //action->limbsize_buffer = NULL;
    //action->limbsize_size = 0;
    action->limbsize_stream = open_memstream(&(action->limbsize_buffer), &(action->limbsize_size));
    
    //action->partition_buffer = NULL;
    //action->partition_size = 0;
    action->partition_stream = open_memstream(&(action->partition_buffer), &(action->partition_size));
    
    action->func = func;
    
    action->n_buffers = 0;
    action->count_a = NULL;
    action->limbsize_buffer_a = NULL;
    action->limbsize_size_a = NULL;
    action->partition_buffer_a = NULL;
    action->partition_size_a = NULL;
}

void schurNumberActionDealloc(schur_number_action_t *action) {
    
    fclose(action->limbsize_stream);
    fclose(action->partition_stream);
    
    free(action->limbsize_buffer);
    free(action->partition_buffer);
    
    if (action->n_buffers > 0) {
        free(action->count_a);
        free(action->limbsize_buffer_a);
        free(action->limbsize_size_a);
        free(action->partition_buffer_a);
        free(action->partition_size_a);
    }
}

void schurNumberActionGatherCopy(schur_number_action_t *action_r, schur_number_action_t *actions_s, size_t n_actions) {
    /*Réunit les n_actions actions_s dans l'unique action_r, en copiant les multiples tampons dans celui de action_r.*/
    unsigned long nmax = action_r->nmax;
    unsigned long iter_num = action_r->iter_num;
    size_t count_all = action_r->count_all;
    size_t count_max = action_r->count_max;
    size_t count = action_r->count;
    FILE * limbsize_stream = action_r->limbsize_stream;
    FILE * partition_stream = action_r->partition_stream;
    void (*func)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) = action_r->func;
    
    for (unsigned long i = 0; i < n_actions; i++) {
        schur_number_action_t *action_s = &actions_s[i];
        
        if (action_s->nmax > nmax) {
            nmax = action_s->nmax;
            count_max = 0;
        }
        
        if (action_s->nmax == nmax) {
            count_max += action_s->count_max;
        }
        
        count_all += action_s->count_all;
        iter_num += action_s->iter_num;
        
        if (func == schurNumberSaveAllPartition) {
            fflush(action_s->limbsize_stream);
            fflush(action_s->partition_stream);
            
            fwrite(action_s->limbsize_buffer, 1, action_s->limbsize_size, limbsize_stream);
            fwrite(action_s->partition_buffer, 1, action_s->partition_size, partition_stream);
            
            count += action_s->count;
        }
    }
    
    if (func == schurNumberSaveBestPartition || func == schurNumberSaveSomePartition) {
        
        if (action_r->nmax < nmax) {
            rewind(limbsize_stream);
            rewind(partition_stream);
            count = 0;
        }
        
        for (unsigned long i = 0; i < n_actions; i++) {
            schur_number_action_t *action_s = &actions_s[i];
            
            if (action_s->nmax == nmax) {
                fflush(action_s->limbsize_stream);
                fflush(action_s->partition_stream);
                
                fwrite(action_s->limbsize_buffer, 1, action_s->limbsize_size, limbsize_stream);
                fwrite(action_s->partition_buffer, 1, action_s->partition_size, partition_stream);
                
                count += action_s->count;
            }
        }
    }
    
    action_r->nmax = nmax;
    action_r->count_all = count_all;
    action_r->count_max = count_max;
    action_r->iter_num = iter_num;
    action_r->count = count;
}

void schurNumberActionGatherNoCopy(schur_number_action_t *action_r, schur_number_action_t *actions_s, size_t n_actions) {
    /*Réunit les n_actions actions_s dans l'unique action_r, sans copier les multiples tampons dans celui de action_r.
     Les actions_s sont automatiquement libérées.*/
    unsigned long nmax = action_r->nmax;
    unsigned long iter_num = action_r->iter_num;
    size_t count_all = action_r->count_all;
    size_t count_max = action_r->count_max;
    size_t count = action_r->count;
    FILE * limbsize_stream = action_r->limbsize_stream;
    FILE * partition_stream = action_r->partition_stream;
    void (*func)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) = action_r->func;
    
    size_t n_bests_action = 0;              // Nombre d'actions ayant atteint nmax
    size_t n_buffers = action_r->n_buffers; // Nombre de tampons déjà présents au sein de action_r
    
    for (unsigned long i = 0; i < n_actions; i++) {
        schur_number_action_t *action_s = &actions_s[i];
        
        if (action_s->nmax > nmax) {
            nmax = action_s->nmax;
            count_max = 0;
            n_bests_action = 0;
        }
        
        if (action_s->nmax == nmax && action_s->count_max) {
            count_max += action_s->count_max;
            n_bests_action++;
        }
        
        count_all += action_s->count_all;
        iter_num += action_s->iter_num;
    }
    
    if (func == schurNumberSaveSomePartition && action_r->nmax < nmax) {
        schur_number_action_t *action_s = actions_s;
        
        while (action_s->nmax < nmax) {
            action_s++;
        }
        
        rewind(limbsize_stream);
        rewind(partition_stream);
        
        fflush(action_s->limbsize_stream);
        fflush(action_s->partition_stream);
        
        fwrite(action_s->limbsize_buffer, 1, action_s->limbsize_size, limbsize_stream);
        fwrite(action_s->partition_buffer, 1, action_s->partition_size, partition_stream);
        
        count = action_s->count;
    }
    
    if (func == schurNumberSaveAllPartition) {
        size_t i0 = n_buffers;
        n_buffers = i0 + n_actions;
        
        size_t *count_a = realloc(action_r->count_a, n_buffers);
        char **limbsize_buffer_a = realloc(action_r->limbsize_buffer_a, n_buffers);
        size_t *limbsize_size_a = realloc(action_r->limbsize_size_a, n_buffers);
        char **partition_buffer_a = realloc(action_r->partition_buffer_a, n_buffers);
        size_t *partition_size_a = realloc(action_r->partition_size_a, n_buffers);
        
        action_r->count_a = count_a;
        action_r->limbsize_buffer_a = limbsize_buffer_a;
        action_r->limbsize_size_a = limbsize_size_a;
        action_r->partition_buffer_a = partition_buffer_a;
        action_r->partition_size_a = partition_size_a;
        
        count_a += i0;
        limbsize_buffer_a += i0;
        limbsize_size_a += i0;
        partition_buffer_a += i0;
        partition_size_a += i0;
        
        for (size_t i = 0; i < n_actions; i++) {
            schur_number_action_t *action_s = &actions_s[i];
            
            fclose(action_s->limbsize_stream);
            fclose(action_s->partition_stream);
            
            *count_a = action_s->count;
            *limbsize_buffer_a = action_s->limbsize_buffer;
            *limbsize_size_a = action_s->limbsize_size;
            *partition_buffer_a = action_s->partition_buffer;
            *partition_size_a = action_s->partition_size;
            
            count_a++;
            limbsize_buffer_a++;
            limbsize_size_a++;
            partition_buffer_a++;
            partition_size_a++;
            
            action_s->limbsize_buffer = NULL;
            action_s->partition_buffer = NULL;
        }
    }
    
    if (func == schurNumberSaveBestPartition) {
        
        if (action_r->nmax < nmax) {
            rewind(limbsize_stream);
            rewind(partition_stream);
            count = 0;
            n_buffers = 0;
        }
        
        size_t i0 = n_buffers;
        n_buffers = i0 + n_bests_action;
        
        size_t *count_a = realloc(action_r->count_a, n_buffers);
        char **limbsize_buffer_a = realloc(action_r->limbsize_buffer_a, n_buffers);
        size_t *limbsize_size_a = realloc(action_r->limbsize_size_a, n_buffers);
        char **partition_buffer_a = realloc(action_r->partition_buffer_a, n_buffers);
        size_t *partition_size_a = realloc(action_r->partition_size_a, n_buffers);
        
        action_r->count_a = count_a;
        action_r->limbsize_buffer_a = limbsize_buffer_a;
        action_r->limbsize_size_a = limbsize_size_a;
        action_r->partition_buffer_a = partition_buffer_a;
        action_r->partition_size_a = partition_size_a;
        
        count_a += i0;
        limbsize_buffer_a += i0;
        limbsize_size_a += i0;
        partition_buffer_a += i0;
        partition_size_a += i0;
        
        for (unsigned long i = 0; i < n_actions; i++) {
            schur_number_action_t *action_s = &actions_s[i];
            
            fclose(action_s->limbsize_stream);
            fclose(action_s->partition_stream);
            
            if (action_s->nmax == nmax && actions_s->count) {
                
                *count_a = action_s->count;
                *limbsize_buffer_a = action_s->limbsize_buffer;
                *limbsize_size_a = action_s->limbsize_size;
                *partition_buffer_a = action_s->partition_buffer;
                *partition_size_a = action_s->partition_size;
                
                //action_s->limbsize_buffer = NULL;
                //action_s->partition_buffer = NULL;
                
                count_a++;
                limbsize_buffer_a++;
                limbsize_size_a++;
                partition_buffer_a++;
                partition_size_a++;
            } else {
                free(action_s->limbsize_buffer);
                free(action_s->partition_buffer);
            }
        }
    }
    
    action_r->nmax = nmax;
    action_r->count_all = count_all;
    action_r->count_max = count_max;
    action_r->iter_num = iter_num;
    action_r->count = count;
    
    action_r->n_buffers = n_buffers;
}

void schurNumberDefaultAction(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /*Met seulement à jour les indicateurs statistiques.*/
    
    if (partition) {
        action->count_all++;
    }
    
    if (n > action->nmax) {
        action->nmax = n;
        action->count_max = 0;
    }
    
    if (n == action->nmax && partition) {
        action->count_max++;
    }
}

void schurNumberSaveSomePartition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /*Compare n avec nmax. Si n ≥ nmax, le nmax est mis à jour et la partition remplace les précédentes.
     Cette action limite la quantité de partition grâce à count_limit.*/
    unsigned long  p = action->p;
    FILE *limbsize_stream = action->limbsize_stream;
    FILE *partition_stream = action->partition_stream;
    
    if (partition) {
        action->count_all++;
    }
    
    if (n > action->nmax) {
        /*Vider partitions*/
        rewind(limbsize_stream);
        rewind(partition_stream);
        action->nmax = n;
        action->count = 0;
        action->count_max = 0;
    }
    
    if (n == action->nmax && partition) {
        if (action->count < action->count_limit) {
            /*Ajouter la partition.*/
            mp_size_t limbsize = ((unsigned long)n>>6) + 1;
            
            fwrite(&limbsize, sizeof(mp_size_t), 1, limbsize_stream);
            
            for (unsigned long j = 0; j < p; j++) {
                fwrite(partition[j], sizeof(mp_limb_t), limbsize, partition_stream);
            }
            
            action->count ++;
        }
        action->count_max ++;
    }
}

void schurNumberSaveBestPartition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /*Compare n avec nmax. Si n ≥ nmax, le nmax est mis à jour et la partition est ajoutée aux partitions.*/
    unsigned long  p = action->p;
    FILE *limbsize_stream = action->limbsize_stream;
    FILE *partition_stream = action->partition_stream;
    
    if (partition) {
        action->count_all++;
    }
    
    if (n > action->nmax) {
        /*Vider partitions*/
        rewind(limbsize_stream);
        rewind(partition_stream);
        action->nmax = n;
        action->count = 0;
        action->count_max = 0;
    }
    
    if (n == action->nmax && partition) {
        /*Ajouter la partition.*/
        mp_size_t limbsize = ((unsigned long)n>>6) + 1;
        
        fwrite(&limbsize, sizeof(mp_size_t), 1, limbsize_stream);
        
        for (unsigned long j = 0; j < p; j++) {
            fwrite(partition[j], sizeof(mp_limb_t), limbsize, partition_stream);
        }
        
        action->count ++;
        action->count_max = action->count;
    }
}

void schurNumberSaveAllPartition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /*La partition est ajoutée aux autres.*/
    
    if (n > action->nmax) {
        action->nmax = n;
        action->count_max = 0;
    }
    
    if (partition) {
        unsigned long  p = action->p;
        FILE *limbsize_stream = action->limbsize_stream;
        FILE *partition_stream = action->partition_stream;
        mp_size_t limbsize = ((action->nmax)>>6) + 1;
        
        fwrite(&limbsize, sizeof(mp_size_t), 1, limbsize_stream);
        
        for (unsigned long j = 0; j < p; j++) {
            fwrite(partition[j], sizeof(mp_limb_t), limbsize, partition_stream);
        }
        
        action->count ++;
        if (n == action->nmax) {
            action->count_max ++;
        }
        action->count_all = action->count;
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

size_t schurNumberPrintPartitions(struct schurNumberIOAction *action) {
    /*Affiche toutes les partitions présentes dans action.*/
    unsigned long p = action->p;
    
    fflush(action->limbsize_stream);
    fflush(action->partition_stream);
    
    size_t total_count = schurNumberPrintPartitionBuffer(p, action->limbsize_buffer, action->partition_buffer, action->count);
    
    for (size_t k = 0; k < action->n_buffers; k++) {
        total_count += schurNumberPrintPartitionBuffer(p, action->limbsize_buffer_a[k], action->partition_buffer_a[k], action->count_a[k]);
    }
    
    return total_count;
}
