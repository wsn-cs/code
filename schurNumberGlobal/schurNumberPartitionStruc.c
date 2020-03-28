//
//  schurNumberPartitionStruc.c
//  schurNumberRecursive
//
//  Created by rubis on 07/02/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberPartitionStruc.h"

void schur_number_partition_alloc(schur_number_partition_t *partitionstruc, mp_size_t limballoc, unsigned long p) {
    /*Initialise partitionstruc.*/
    
    mp_limb_t **partition = calloc(p, sizeof(mp_limb_t *));
    mp_limb_t **partitioninvert = calloc(p, sizeof(mp_limb_t *));
    for (unsigned long j = 0; j < p; j++) {
        //partition[j] = calloc(limballoc, sizeof(mp_limb_t));
        partitioninvert[j] = calloc(2 * limballoc, sizeof(mp_limb_t));
        partition[j] = partitioninvert[j] + limballoc;
    }
    
    partitionstruc->p = 0;
    partitionstruc->pmax = p;
    partitionstruc->n = 0;
    partitionstruc->limballoc = limballoc;
    partitionstruc->limbsize = 1;
    partitionstruc->partition = partition;
    partitionstruc->partitioninvert = partitioninvert;
}

void schur_number_partition_dealloc(schur_number_partition_t *partitionstruc) {
    /*Libère partitionstruc.*/
    unsigned long p = partitionstruc->pmax;
    
    mp_limb_t **partition = partitionstruc->partition;
    mp_limb_t **partitioninvert = partitionstruc->partitioninvert;
    
    for (unsigned long j = 0; j < p; j++) {
        free(partitioninvert[j]);
    }
    free(partition);
    free(partitioninvert);
}
