//
//  schurNumberPartitionStruc.h
//  schurNumberRecursive
//
//  Created by rubis on 07/02/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#ifndef schurNumberPartitionStruc_h
#define schurNumberPartitionStruc_h

#include <stdlib.h>
#include <gmp.h>

#define ADD_POINT(set, x) set[(x) / mp_bits_per_limb] |= ((unsigned long)1 << ((x) % mp_bits_per_limb))

#define DELETE_POINT(set, x) set[(x) / mp_bits_per_limb] ^= ((unsigned long)1 << ((x) % mp_bits_per_limb))

#define GET_POINT(set, x) (set[(x) / mp_bits_per_limb] & ((unsigned long)1 << ((x) % mp_bits_per_limb)))

struct schur_number_partition_struc {
    unsigned long p;     // Nombre de huches non vides
    unsigned long pmax;  // Nombre de huches allouées
    
    unsigned long n;    // Entier courant
    mp_size_t limballoc;// Nombre de limbes alloués à chaque huche
    mp_size_t limbsize; // Nombre de limbes utilisés par chaque huche
    
    mp_limb_t **partition;
    mp_limb_t **partitioninvert;
};

typedef struct schur_number_partition_struc schur_number_partition_t;

void schur_number_partition_alloc(schur_number_partition_t *partitionstruc, mp_size_t limballoc, unsigned long p);

void schur_number_partition_dealloc(schur_number_partition_t *partitionstruc);

#endif /* schurNumberPartitionStruc_h */
