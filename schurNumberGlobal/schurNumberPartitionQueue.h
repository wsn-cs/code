//
//  schurNumberPartitionQueue.h
//  SchurNumber
//
//  Created by rubis on 17/04/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#ifndef schurNumberPartitionQueue_h
#define schurNumberPartitionQueue_h

#include "schurNumberPartitionStruc.h"

enum partition_queue_flag_enum {
    NULL_FLAG,              // Drapeau indiquant l'absence de partitions dans la file
    INITIAL_FLAG,
    INTERMEDIATE_FLAG,
    STOP_FLAG,              // Drapeau demandant l'arrêt des recherches
};

typedef enum partition_queue_flag_enum partition_queue_flag_t;

struct schur_number_partition_queue_struc {
    unsigned long pmax;         // Nombre de parties de chaque partition
    
    unsigned long n;            // Taille de l'intervalle
    
    size_t count;               // Nombre de partitions initiales dans partition_array
    size_t current_index;       // Indice de la partition initiale à sélectionner
    size_t total_count;         // Nombre totale de partitions dans la file
    mp_size_t *limbsize_array;  // Tableau contenant les nombres de limbes par ensembles de la partition
    mp_size_t *limbsize_ptr;    // Pointeur vers la taille actuelle
    mp_limb_t *partition_array; // Tableau des partitions initiales
    mp_limb_t *partition_ptr;   // Pointeur vers la partition actuelle
    
    partition_queue_flag_t flag;// Drapeau permettant de distinguer les phases d'une procédure
    
    struct schur_number_partition_queue_struc *child_queue; // Pointeur vers une file-fille
};

typedef struct schur_number_partition_queue_struc schur_number_partition_queue_t;

void schur_number_partition_queue_init(schur_number_partition_queue_t *partition_queue);
void schur_number_partition_queue_dealloc(schur_number_partition_queue_t *partition_queue);

size_t schur_partition_queue_add_partition_copy(schur_number_partition_queue_t *partition_queue, schur_number_partition_t *partition_struc, partition_queue_flag_t flag);
size_t schur_partition_queue_add_partitionarray_nocopy(schur_number_partition_queue_t *partition_queue, mp_limb_t *partition_array, mp_size_t *limbsize_array, size_t count, partition_queue_flag_t flag);

partition_queue_flag_t schur_partition_queue_get_partition(schur_number_partition_queue_t *partition_queue, schur_number_partition_t *partition);

#endif /* schurNumberPartitionQueue_h */
