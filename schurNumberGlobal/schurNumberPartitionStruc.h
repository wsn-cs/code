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

#define ADD_POINT(set, x) (set)[(x) / GMP_NUMB_BITS] |= ((unsigned long)1 << ((x) % GMP_NUMB_BITS))

#define DELETE_POINT(set, x) (set)[(x) / GMP_NUMB_BITS] ^= ((unsigned long)1 << ((x) % GMP_NUMB_BITS))

#define GET_POINT(set, x) ((set)[(x) / GMP_NUMB_BITS] & ((unsigned long)1 << ((x) % GMP_NUMB_BITS)))

#define PARTITION_2_LIMBSIZE(p) (((1 << (2 * (p) - 1)) >> 6) + 1)

#define INTEGER_2_LIMBSIZE(n) ((n) / GMP_NUMB_BITS) + 1

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

void schur_number_partition_alloc(schur_number_partition_t *partitionstruc, unsigned long p);
void schur_number_partition_init(schur_number_partition_t *partitionstruc, mp_size_t limballoc);

void schur_number_partition_dealloc(schur_number_partition_t *partitionstruc);

unsigned long schurNumberPartitionSetString(schur_number_partition_t *partition, char **str, size_t str_size, int format);

void schurNumberSumset(mp_limb_t *r_set, mp_limb_t *set1, mp_limb_t *set2, mp_size_t limbsize, unsigned long x);
void schurNumberWeakSumset(mp_limb_t *r_set, mp_limb_t *set1, mp_limb_t *set2, mp_size_t limbsize, unsigned long x);

void schur_number_set_revert(mp_limb_t *r_set, mp_limb_t *set, mp_size_t limbsize);

static inline void schur_number_translation(mp_limb_t *r_set, const mp_limb_t *set, mp_size_t limbsize, unsigned long n) {
    /* Calcule set + n et place le résultat dans r_set. */
    const mp_limb_t *work = set;
    
    while (n > 0) {
        unsigned int shift = n % GMP_NUMB_BITS;
        if (!shift) {
            shift = GMP_NUMB_BITS - 1;
        }
        
        mpn_lshift(r_set, work, limbsize, shift);     // r_set = work + shift
        work = r_set;
        
        n -= shift;
    }
}

static inline void schur_number_ntranslation(mp_limb_t *r_set, const mp_limb_t *set, mp_size_t limbsize, unsigned long n) {
    /* Calcule set - n et place le résultat dans r_set. */
    const mp_limb_t *work = set;
    
    while (n > 0) {
        unsigned int shift = n % GMP_NUMB_BITS;
        if (!shift) {
            shift = GMP_NUMB_BITS - 1;
        }
        
        mpn_rshift(r_set, work, limbsize, shift);     // r_set = work - shift
        work = r_set;
        
        n -= shift;
    }
}

static inline void schur_number_setinterval_1(mp_limb_t *set, mp_size_t limbsize, unsigned long n) {
    /* Place l'intervalle [1, n-1] = 2^n - 1 dans set. */
    mp_size_t blockingsize = (n >> 6) + 1;
    
    mpn_zero(set, limbsize);
    *set = (mp_limb_t)1;
    if (blockingsize > 1) {
        mpn_neg(set, set, blockingsize - 1);
    }
    set[blockingsize - 1] = ((mp_limb_t)1 << (n % GMP_NUMB_BITS)) - 1;
}

static inline void schur_number_setinterval_s(mp_limb_t *set, mp_size_t limbsize, unsigned long n) {
    /* Place l'intervalle [n+1, limbsize * GMP_NUMB_BITS - 1] = -2^n dans set. */
    mp_size_t blockingsize = (n >> 6) + 1;
    
    mpn_zero(set, limbsize);
    mp_limb_t *set0 = set + limbsize - blockingsize;
    *set0 = (mp_limb_t)1 << (GMP_NUMB_BITS - (n % GMP_NUMB_BITS));
    mpn_neg(set0, set0, blockingsize);
}

#endif /* schurNumberPartitionStruc_h */
