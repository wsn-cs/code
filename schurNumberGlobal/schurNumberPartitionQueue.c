//
//  schurNumberPartitionQueue.c
//  SchurNumber
//
//  Created by rubis on 17/04/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberPartitionQueue.h"

void schur_number_partition_queue_init(schur_number_partition_queue_t *partition_queue) {
    /* Initialise la file. */
    partition_queue->pmax = 0;
    partition_queue->count = 0;
    partition_queue->current_index = 0;
    partition_queue->total_count = 0;
    
    partition_queue->limbsize_array = NULL;
    partition_queue->limbsize_ptr = NULL;
    partition_queue->partition_array = NULL;
    partition_queue->partition_ptr = NULL;
    
    partition_queue->flag = NULL_FLAG;
    
    partition_queue->child_queue = NULL;
    partition_queue->leaf_queue = partition_queue;
}

void schur_number_partition_queue_dealloc(schur_number_partition_queue_t *partition_queue) {
    /* Libère la mémoire associée aux éléments de la file, et nettoie recursivement ses enfants. */
    mp_size_t *limbsize_array = partition_queue->limbsize_array;
    mp_limb_t *partition_array = partition_queue->partition_array;
    
    if (limbsize_array) {
        free(limbsize_array);
        partition_queue->limbsize_array = NULL;
        partition_queue->limbsize_ptr = NULL;
    }
    
    if (partition_array) {
        free(partition_array);
        partition_queue->partition_array = NULL;
        partition_queue->partition_ptr = NULL;
    }
    
    if (partition_queue->child_queue) {
        schur_number_partition_queue_dealloc(partition_queue->child_queue);
        free(partition_queue->child_queue);
        partition_queue->child_queue = NULL;
    }
}

size_t schur_partition_queue_add_partition_copy(schur_number_partition_queue_t *partition_queue, schur_number_partition_t *partition_struc, partition_queue_flag_t flag) {
    /* Copie la partition dans la file. */
    schur_number_partition_queue_t *queue = partition_queue;
    size_t total_count = queue->total_count + 1;
    queue->total_count = total_count;
    /*while (queue->child_queue) {
        queue = queue->child_queue; // Récupérer le plus petit descendant de partition_queue
        queue->total_count ++;
    }*/
    schur_number_partition_queue_t *child = calloc(1, sizeof(schur_number_partition_queue_t));
    
    schur_number_partition_queue_t *leaf = queue->leaf_queue;
    queue->leaf_queue = child;
    queue = leaf;
    queue->leaf_queue = child;
    
    unsigned long pmax = partition_struc->pmax;
    mp_size_t limbsize = partition_struc->limbsize;
    
    child->pmax = pmax;
    child->n = partition_struc->n;
    
    child->count = 1;
    child->current_index = 0;
    child->total_count = 1;
    
    mp_size_t *limbsize_ptr = calloc(1, sizeof(mp_size_t));
    *limbsize_ptr = limbsize;
    child->limbsize_array = limbsize_ptr;
    child->limbsize_ptr = child->limbsize_array;
    
    mp_limb_t *partition = calloc(pmax, limbsize * sizeof(mp_limb_t));
    for (unsigned long i = 0; i < partition_struc->p; i++) {
        mpn_copyd(&partition[i * limbsize], partition_struc->partition[i], limbsize);
    }
    child->partition_array = partition;
    child->partition_ptr = child->partition_array;
    
    child->flag = flag;
    
    child->child_queue = NULL;
    child->leaf_queue = child;
    
    queue->child_queue = child;
    
    return total_count;
}

size_t schur_partition_queue_add_partitionarray_nocopy(schur_number_partition_queue_t *partition_queue, mp_limb_t *partition_array, mp_size_t *limbsize_array, size_t count, unsigned long n, partition_queue_flag_t flag) {
    /* Ajoute les count partitions contenues dans partition_array à la file. Celle-ci en devient propriétaire et est responsable de leur libération.
     Cette fonction renvoie le nouveau nombre total de partitions dans la file. */
    
    schur_number_partition_queue_t *queue = partition_queue;
    size_t total_count = queue->total_count + count;
    queue->total_count = total_count;
    /*while (queue->child_queue) {
        queue = queue->child_queue; // Récupérer le plus petit descendant de partition_queue
        queue->total_count += count;
    }*/
    
    schur_number_partition_queue_t *child = calloc(1, sizeof(schur_number_partition_queue_t));
    
    schur_number_partition_queue_t *leaf = queue->leaf_queue;
    queue->leaf_queue = child;
    queue = leaf;
    queue->leaf_queue = child;
    
    child->pmax = partition_queue->pmax;
    child->n = n;
    
    child->count = count;
    child->current_index = 0;
    child->total_count = count;
    
    child->limbsize_array = limbsize_array;
    child->limbsize_ptr = child->limbsize_array;
    
    child->partition_array = partition_array;
    child->partition_ptr = child->partition_array;
    
    child->flag = flag;
    
    child->child_queue = NULL;
    child->leaf_queue = child;
    
    queue->child_queue = child;
    
    return total_count;
}

partition_queue_flag_t schur_partition_queue_get_partition(schur_number_partition_queue_t *partition_queue, schur_number_partition_t *partition) {
    /* Fournit une partition à traiter et renvoie flag, ou 0 si il n'en reste plus. */
    size_t index = partition_queue->current_index;
    
    while (index >= partition_queue->count) {
        schur_number_partition_queue_t *child = partition_queue->child_queue;
        
        if (!child) {
            char flag = NULL_FLAG;
            if (partition_queue->flag == STOP_FLAG) {
                flag = STOP_FLAG;
            }
            return flag;
        }
        
        partition_queue->pmax = child->pmax;
        partition_queue->n = child->n;
        
        partition_queue->count = child->count;
        partition_queue->current_index = child->current_index;
        partition_queue->total_count = child->total_count;
        
        if (partition_queue->limbsize_array) {
            free(partition_queue->limbsize_array);
        }
        partition_queue->limbsize_array = child->limbsize_array;
        partition_queue->limbsize_ptr = child->limbsize_ptr;
        
        if (partition_queue->partition_array) {
            free(partition_queue->partition_array);
        }
        partition_queue->partition_array = child->partition_array;
        partition_queue->partition_ptr = child->partition_ptr;
        
        partition_queue->flag = child->flag;
        
        partition_queue->child_queue = child->child_queue;
        
        if (partition_queue->leaf_queue == child) {
            partition_queue->leaf_queue = partition_queue;
        }
        free(child);
        
        index = partition_queue->current_index;
    }
    
    unsigned long pmax = partition_queue->pmax;
    mp_size_t limballoc = partition->limballoc;
    
    schur_number_partition_set_empty(partition);
    
    unsigned long i = 0;
    mp_size_t limbsize = *partition_queue->limbsize_ptr;
    mp_limb_t *set = partition_queue->partition_ptr;
    mp_limb_t *set0 = set;      // set0 contiendra l'union des ensembles de la partition
    
    // Copier tous les ensembles dans la structure partition
    while (i < pmax && !mpn_zero_p(set, limbsize)) {
        mpn_copyd(partition->partition[i], set, limbsize);
        schur_number_set_revert(partition->partitioninvert[i], partition->partition[i], limballoc);
        mpn_ior_n(set0, set0, set, limbsize);  // Réaliser l'union des ensembles pour déterminer n
        i++;
        set += limbsize;
    }
    
    partition->p = i;
    unsigned long n0 = partition_queue->n;
    // Déterminer le plus grand entier n tel que [n0, n] soit recouvert
    mp_size_t q0 = n0 / GMP_NUMB_BITS;
    //schur_number_intersect_interval_n(set0, limbsize, limbsize0, n0 + 1);
    set0 += q0;
    mpn_com(set0, set0, limbsize - q0);
    *set0 &= (GMP_NUMB_MAX << (n0 % GMP_NUMB_BITS));
    if (mpn_zero_p(set0, limbsize - q0)) {
        partition->n = limbsize * GMP_NUMB_BITS - 1;
        partition->limbsize = limbsize;
    } else {
        partition->n = mpn_scan1(set0, n0 % GMP_NUMB_BITS) + q0 * GMP_NUMB_BITS - 1;
        partition->limbsize = INTEGER_2_LIMBSIZE(partition->n);
    }
    
    partition_queue->current_index = index + 1;
    partition_queue->limbsize_ptr ++;
    partition_queue->partition_ptr = set + limbsize * (pmax - i);
    partition_queue->total_count --;
    
    return partition_queue->flag;
}
