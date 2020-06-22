//
//  schurNumberDefaultAction.c
//  SchurNumber
//
//  Created by rubis on 21/06/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberIOAction.h"

unsigned long schur_number_default_action(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /*Met seulement à jour les indicateurs statistiques.*/
    
    schur_number_intermediate_save_t *save = action->save;
    
    if (partition) {
        action->count_all++;
        
        if (save && !(action->count_all % SCHUR_NUMBER_SAVE_FREQUENCY)) {
            unsigned long nbest = schur_number_save_progression_update(save, n, partition);
            if (nbest > action->nbest) {
                schur_number_save_some_partition(NULL, nbest, action);
            }
        }
    }
    
    if (n > action->nbest) {
        if (save) {
            action->nbest = schur_number_save_best_upgrade(save, n, partition);
        } else {
            action->nbest = n;
        }
        action->count_max = 0;
    }
    
    if (n == action->nmax && partition) {
        action->count_max++;
    }
    
    return action->nmax;
}

unsigned long schur_number_save_some_partition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /*Compare n avec nbest. Si n ≥ nbest, le nbest est mis à jour et la partition remplace les précédentes.
     Cette action limite la quantité de partition grâce à count_limit.*/
    unsigned long  p = action->p;
    FILE *limbsize_stream = action->limbsize_stream;
    FILE *partition_stream = action->partition_stream;
    
    schur_number_intermediate_save_t *save = action->save;
    
    if (partition) {
        action->count_all++;
        
        if (save && !(action->count_all % SCHUR_NUMBER_SAVE_FREQUENCY)) {
            unsigned long nbest = schur_number_save_progression_update(save, n, partition);
            if (nbest > action->nbest) {
                schur_number_save_some_partition(NULL, nbest, action);
            }
        }
    }
    
    if (n > action->nbest) {
        /*Vider partitions*/
        rewind(limbsize_stream);
        rewind(partition_stream);
        if (save) {
            action->nbest = schur_number_save_best_upgrade(save, n, partition);
        } else {
            action->nbest = n;
        }
        action->count = 0;
        action->count_max = 0;
    }
    
    if (n == action->nbest && partition) {
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
    
    return action->nbest;
}

unsigned long schur_number_save_best_partition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /*Compare n avec nbest. Si n ≥ nbest, le nbest est mis à jour et la partition est ajoutée aux partitions.*/
    unsigned long  p = action->p;
    FILE *limbsize_stream = action->limbsize_stream;
    FILE *partition_stream = action->partition_stream;
    
    schur_number_intermediate_save_t *save = action->save;
    
    if (partition) {
        action->count_all++;
        
        if (save && !(action->count_all % SCHUR_NUMBER_SAVE_FREQUENCY)) {
            unsigned long nbest = schur_number_save_progression_update(save, n, partition);
            if (nbest > action->nbest) {
                schur_number_save_some_partition(NULL, nbest, action);
            }
        }
    }
    
    if (n > action->nbest) {
        /*Vider partitions*/
        rewind(limbsize_stream);
        rewind(partition_stream);
        if (save) {
            action->nbest = schur_number_save_best_upgrade(save, n, partition);
        } else {
            action->nbest = n;
        }
        action->count = 0;
        action->count_max = 0;
    }
    
    if (n == action->nbest && partition) {
        if (action->partition_size < action->size_limit) {
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
    
    return action->nbest;
}

unsigned long schur_number_save_all_partition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /*La partition est ajoutée aux autres.*/
    
    schur_number_intermediate_save_t *save = action->save;
    
    if (n > action->nbest) {
        if (save) {
            action->nbest = schur_number_save_best_upgrade(save, n, partition);
        } else {
            action->nbest = n;
        }
        action->count_max = 0;
    }
    
    if (partition) {
        if (action->partition_size < action->size_limit) {
            /* Ajouter la partition. */
            unsigned long  p = action->p;
            FILE *limbsize_stream = action->limbsize_stream;
            FILE *partition_stream = action->partition_stream;
            mp_size_t limbsize = ((action->nbest)>>6) + 1;
            
            fwrite(&limbsize, sizeof(mp_size_t), 1, limbsize_stream);
            
            for (unsigned long j = 0; j < p; j++) {
                fwrite(partition[j], sizeof(mp_limb_t), limbsize, partition_stream);
            }
            
            action->count ++;
        }
        if (n == action->nbest) {
            action->count_max ++;
        }
        action->count_all ++;
        
        if (save && !(action->count_all % SCHUR_NUMBER_SAVE_FREQUENCY)) {
            unsigned long nbest = schur_number_save_progression_update(save, n, partition);
            if (nbest > action->nbest) {
                schur_number_save_some_partition(NULL, nbest, action);
            }
        }
    }
    
    return action->nbest;
}
