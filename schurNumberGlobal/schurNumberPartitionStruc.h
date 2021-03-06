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

#define DELETE_POINT(set, x) (set)[(x) / GMP_NUMB_BITS] &= ~((unsigned long)1 << ((x) % GMP_NUMB_BITS))

#define GET_POINT(set, x) ((set)[(x) / GMP_NUMB_BITS] & ((unsigned long)1 << ((x) % GMP_NUMB_BITS)))

#define PARTITION_2_LIMBSIZE(p) (((1 << (2 * (p) - 1)) / GMP_NUMB_BITS) + 1)

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

void schur_number_partition_set_empty(schur_number_partition_t *partitionstruc);

unsigned long schurNumberPartitionSetString(schur_number_partition_t *partition, char **str, size_t str_size, int format);

void schur_number_sumset(mp_limb_t *r_set, mp_limb_t *set1, mp_limb_t *set2, mp_size_t r_limbsize, mp_size_t limbsize, unsigned long x, mp_limb_t *work);
void schur_number_restricted_sumset(mp_limb_t *r_set, mp_limb_t *set1, mp_limb_t *set2, mp_size_t r_limbsize, mp_size_t limbsize, unsigned long x, mp_limb_t *work);

void schur_number_set_revert(mp_limb_t *r_set, mp_limb_t *set, mp_size_t limbsize);

static inline void schur_number_translation(mp_limb_t *r_set, const mp_limb_t *s_set, mp_size_t limbsize, unsigned long n) {
    /* Calcule s_set + n et place le résultat dans r_set. Cela suppose que r_set possède lui aussi limbsize limbes.*/
    
    // Division euclidenne n = q * GMP_NUMB_BITS + r
    unsigned long q = n / GMP_NUMB_BITS;
    unsigned int r = n % GMP_NUMB_BITS;
    
    // s_set + n = (s_set + q * GMP_NUMB_BITS) + r
    if (limbsize > q) {
        mpn_lshift(r_set + q, s_set, limbsize - q, r);
    } else {
        mpn_zero(r_set, limbsize);
    }
}

static inline void schur_number_ntranslation(mp_limb_t *r_set, const mp_limb_t *s_set, mp_size_t limbsize, unsigned long n) {
    /* Calcule set - n et place le résultat dans r_set. */
    
    // Division euclidenne n = q * GMP_NUMB_BITS + r
    unsigned long q = n / GMP_NUMB_BITS;
    unsigned int r = n % GMP_NUMB_BITS;
    
    // s_set - n = (s_set - q * GMP_NUMB_BITS) - r
    mpn_zero(r_set + limbsize - q, q);
    mpn_rshift(r_set, s_set + q, limbsize - q, r);
}

static inline void schur_number_setinterval_1(mp_limb_t *set, mp_size_t limbsize, mp_size_t blockingsize, unsigned long n) {
    /* Place l'intervalle [1, n-1] = 2^n - 1 dans set. */
    //mp_size_t blockingsize = (n >> 6) + 1;
    
    mpn_zero(set, limbsize);
    *set = (mp_limb_t)1;
    if (blockingsize > 1) {
        mpn_neg(set, set, blockingsize - 1);
    }
    set[blockingsize - 1] = ((mp_limb_t)1 << (n % GMP_NUMB_BITS)) - 1;
}

static inline void schur_number_intersect_interval_0(mp_limb_t *set, mp_size_t limbsize, mp_size_t blockingsize, unsigned long n) {
    /* Intersecte set avec l'intervalle [0, n], sous l'hypothèse que blockingsize <= limbsize. */
    //mp_size_t blockingsize = (n >> 6) + 1;
    if (blockingsize < limbsize) {
        mpn_zero(set + blockingsize, limbsize - blockingsize);
    }
    set[blockingsize - 1] &= (GMP_NUMB_MAX >> ((GMP_NUMB_BITS - n - 1) % GMP_NUMB_BITS));
}

static inline void schur_number_intersect_interval_n(mp_limb_t *set, mp_size_t limbsize, mp_size_t blockingsize, unsigned long n) {
    /* Intersecte set avec l'intervalle [n, +∞[. */
    //mp_size_t blockingsize = (n >> 6) + 1;
    if (blockingsize > 1) {
        mpn_zero(set, blockingsize - 1);
    }
    set[blockingsize - 1] &= (GMP_NUMB_MAX << (n % GMP_NUMB_BITS));
}

static inline void schur_number_setinterval_s(mp_limb_t *set, mp_size_t limbsize, mp_size_t blockingsize, unsigned long n) {
    /* Place l'intervalle [n+1, limbsize * GMP_NUMB_BITS - 1] = -2^n dans set. Cette fonction ne marche pas si n = 1. */
    //mp_size_t blockingsize = (n >> 6) + 1;
    
    mpn_zero(set, limbsize);
    if (n > 1) {
        ADD_POINT(set, GMP_NUMB_BITS * limbsize - n + 1);
        mp_limb_t *set0 = set + limbsize - blockingsize;
        mpn_neg(set0, set0, blockingsize);
    }
}

#endif /* schurNumberPartitionStruc_h */
