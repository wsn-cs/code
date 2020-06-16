//
//  schurNumberSaveDistinctSumPartition.c
//  SchurNumber
//
//  Created by rubis on 07/05/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberIOAction.h"
#include "schurNumberPartitionStruc.h"

static inline int find_sequential_index(mp_limb_t *sumset_array, mp_limb_t *sumset, mp_size_t limbsize, size_t *sorted_indexes, size_t size, unsigned long *i) {
    int cmp_r = 0;
    
    for (*i = 0; *i < size; (*i) ++) {
        cmp_r = mpn_cmp(&(sumset_array[limbsize * sorted_indexes[*i]]), sumset, limbsize);
        if (!cmp_r) {
            break;
        }
        if (cmp_r > 0) {
            (*i)--;
            break;
        }
    }
    
    return cmp_r;
}

static inline void find_dichotomic_index(mp_limb_t *sumset_array, mp_limb_t *sumset, mp_size_t limbsize, size_t *sorted_indexes, size_t size, size_t *i1_p, size_t *i2_p) {
    /* Cette fonction trouve les deux indices du tableau sumset_array encadrant sumset. */
    int cmp_r;
    *i1_p = 0;
    *i2_p = size;
    
    cmp_r = mpn_cmp(&(sumset_array[limbsize * sorted_indexes[0]]), sumset, limbsize);
    if (cmp_r >= 0) {
        *i2_p = 0;
        if (cmp_r) {
            *i1_p = -1;
        }
        return;
    }
    
    cmp_r = mpn_cmp(&(sumset_array[limbsize * sorted_indexes[size - 1]]), sumset, limbsize);
    if (cmp_r <= 0) {
        if (cmp_r) {
            *i1_p = size;
            *i2_p = size + 1;
        } else {
            *i1_p = size-1;
            *i2_p = *i1_p;
        }
        return;
    }
    
    while (*i2_p - *i1_p > 1) {
        unsigned long imedian = (*i1_p + *i2_p) >> 1;
        cmp_r = mpn_cmp(&(sumset_array[limbsize * sorted_indexes[imedian]]), sumset, limbsize);
        
        if (cmp_r <= 0) {
            *i1_p = imedian;
        }
        if (cmp_r >= 0) {
            *i2_p = imedian;
        }
    }
}

unsigned long schur_number_save_distinct_sum_partition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /* Ajoute la partition si sa somme n'a pas déjà été trouvée, et si n == nmax. */
    unsigned long  p = action->p;
    FILE *limbsize_stream = action->limbsize_stream;
    FILE *partition_stream = action->partition_stream;
    
    FILE *sum_partition_stream = action->sum_partition_stream;
    FILE *sorted_index_sum_partition_stream = action->sorted_index_sum_partition_stream;
    
    schur_number_intermediate_save_t *save = action->save;
    
    if (partition) {
        action->count_all++;
        
        if (save && !(action->count_all % SCHUR_NUMBER_SAVE_FREQUENCY)) {
            unsigned long nbest = schur_number_save_progression_update(save, n, partition);
            if (nbest > action->nmax) {
                schur_number_save_distinct_sum_partition(NULL, nbest, action);
            }
        }
    }
    
    if (n > action->nmax) {
        /*Vider partitions*/
        rewind(limbsize_stream);
        rewind(partition_stream);
        rewind(sum_partition_stream);
        rewind(sorted_index_sum_partition_stream);
        if (save) {
            action->nmax = schur_number_save_best_upgrade(save, n, partition);
        } else {
            action->nmax = n;
        }
        action->count = 0;
        action->count_max = 0;
    }
    
    if (n == action->nmax && partition && action->count < action->count_limit) {
        if (action->partition_size < action->size_limit) {
            mp_limb_t *work = action->work;
            
            /* Regarder si la somme de la partition a déjà été trouvée. */
            mp_size_t limbsize = ((unsigned long)n>>6) + 1;
            schur_number_sumset(work, partition[p - 1], partition[p - 1], 2 * limbsize, limbsize, 0, work + 2 * limbsize);
            
            fflush(sum_partition_stream);
            fflush(sorted_index_sum_partition_stream);
            
            size_t *sorted_index_sum_partition_buffer = action->sorted_index_sum_partition_buffer;
            
            size_t i1 = 0;
            size_t i2 = action->count;
            
            if (action->count > 0) {
                find_dichotomic_index(action->sum_partition_buffer, work, 2 * limbsize, sorted_index_sum_partition_buffer, action->count, &i1, &i2);
            } else {
                i2 = 1;
            }
            
            if (i1 != i2) {
                /* Somme pas encore présente, donc ajouter la partition. */
                fwrite(&limbsize, sizeof(mp_size_t), 1, limbsize_stream);
                
                for (unsigned long j = 0; j < p; j++) {
                    fwrite(partition[j], sizeof(mp_limb_t), limbsize, partition_stream);
                }
                
                /* Ajouter sa somme. */
                fwrite(work, sizeof(mp_limb_t), 2 * limbsize, sum_partition_stream);
                
                /* Placer son indice. */
                size_t index = action->count;
                for (size_t i = i1 + 1; i < action->count; i++) {
                    size_t old_index = sorted_index_sum_partition_buffer[i];
                    sorted_index_sum_partition_buffer[i] = index;
                    index = old_index;
                }
                fwrite(&index, sizeof(size_t), 1, sorted_index_sum_partition_stream);
                
                action->count ++;
            } else {
                /* Compter le nombre d'éléments de partition[p] et remplacer la partition déjà enregistrée si il lui est supérieur. */
                fflush(action->partition_stream);
                mp_limb_t *partition_buffer = action->partition_buffer;
                size_t index = sorted_index_sum_partition_buffer[i1];
                if (mpn_popcount(partition[p-1], limbsize) > mpn_popcount(&(partition_buffer[limbsize * (p * index + p - 1)]), limbsize)) {
                    for (unsigned long j = 0; j < p; j++) {
                        mpn_copyi(&(partition_buffer[limbsize * (p * index + j)]), partition[j], limbsize);
                    }
                }
            }
        }
        action->count_max ++;
    }
    
    return action->nmax;
}

unsigned long schur_number_save_distinct_restrictedsum_partition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /* Ajoute la partition si sa somme restreinte n'a pas déjà été trouvée, et si n == nmax. */
    unsigned long  p = action->p;
    FILE *limbsize_stream = action->limbsize_stream;
    FILE *partition_stream = action->partition_stream;
    
    FILE *sum_partition_stream = action->sum_partition_stream;
    FILE *sorted_index_sum_partition_stream = action->sorted_index_sum_partition_stream;
    
    schur_number_intermediate_save_t *save = action->save;
    
    if (partition) {
        action->count_all++;
        
        if (save && !(action->count_all % SCHUR_NUMBER_SAVE_FREQUENCY)) {
            unsigned long nbest = schur_number_save_progression_update(save, n, partition);
            if (nbest > action->nmax) {
                schur_number_save_distinct_sum_partition(NULL, nbest, action);
            }
        }
    }
    
    if (n > action->nmax) {
        /*Vider partitions*/
        rewind(limbsize_stream);
        rewind(partition_stream);
        rewind(sum_partition_stream);
        rewind(sorted_index_sum_partition_stream);
        if (save) {
            action->nmax = schur_number_save_best_upgrade(save, n, partition);
        } else {
            action->nmax = n;
        }
        action->count = 0;
        action->count_max = 0;
    }
    
    if (n == action->nmax && partition) {
        if (action->partition_size < action->size_limit) {
            mp_limb_t *work = action->work;
            
            /* Regarder si la somme de la partition a déjà été trouvée. */
            mp_size_t limbsize = action->limbsize;
            schur_number_restricted_sumset(work, partition[p - 1], partition[p - 1], 2 * limbsize, limbsize, 0, &(work[2 * limbsize]));
            //schurNumberWeakSumset2(work, partition[p - 1], 2 * limbsize, limbsize, work + 2 * limbsize);
            
            fflush(sum_partition_stream);
            fflush(sorted_index_sum_partition_stream);
            
            mp_limb_t *sorted_index_sum_partition_buffer = action->sorted_index_sum_partition_buffer;
            
            size_t i1 = 0;
            size_t i2 = action->count;
            
            if (action->count > 0) {
                find_dichotomic_index(action->sum_partition_buffer, work, 2 * limbsize, sorted_index_sum_partition_buffer, action->count, &i1, &i2);
            } else {
                i2 = 1;
            }
            
            if (i1 != i2) {
                /* Somme pas encore présente, donc ajouter la partition. */
                fwrite(&limbsize, sizeof(mp_size_t), 1, limbsize_stream);
                
                for (unsigned long j = 0; j < p; j++) {
                    fwrite(partition[j], sizeof(mp_limb_t), limbsize, partition_stream);
                }
                
                /* Ajouter sa somme. */
                fwrite(work, sizeof(mp_limb_t), 2 * limbsize, sum_partition_stream);
                
                /* Placer son indice. */
                size_t index = action->count;
                for (size_t i = i1 + 1; i < action->count; i++) {
                    size_t old_index = sorted_index_sum_partition_buffer[i];
                    sorted_index_sum_partition_buffer[i] = index;
                    index = old_index;
                }
                fwrite(&index, sizeof(size_t), 1, sorted_index_sum_partition_stream);
                
                action->count ++;
            } else {
                /* Compter le nombre d'éléments de partition[p] et remplacer la partition déjà enregistrée si il lui est supérieur. */
                fflush(action->partition_stream);
                mp_limb_t *partition_buffer = action->partition_buffer;
                size_t index = sorted_index_sum_partition_buffer[i1];
                if (mpn_popcount(partition[p-1], limbsize) > mpn_popcount(&(partition_buffer[limbsize * (p * index + p - 1)]), limbsize)) {
                    for (unsigned long j = 0; j < p; j++) {
                        mpn_copyi(&(partition_buffer[limbsize * (p * index + j)]), partition[j], limbsize);
                    }
                }
            }
        }
        action->count_max ++;
    }
    
    return action->nmax;
}

