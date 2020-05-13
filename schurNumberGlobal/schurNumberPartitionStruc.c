//
//  schurNumberPartitionStruc.c
//  schurNumberRecursive
//
//  Created by rubis on 07/02/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberPartitionStruc.h"

void schur_number_partition_alloc(schur_number_partition_t *partitionstruc, unsigned long p) {
    /* Alloue partitionstruc, c'est-à-dire alloue le tableau pour placer des pointeurs vers chaque ensemble de la partition. */
    
    mp_limb_t **partition = calloc(p, sizeof(mp_limb_t *));
    mp_limb_t **partitioninvert = calloc(p, sizeof(mp_limb_t *));
    
    for (unsigned long j = 0; j < p; j++) {
        partitioninvert[j] = NULL;
        partition[j] = NULL;
    }
    
    partitionstruc->p = 0;
    partitionstruc->pmax = p;
    partitionstruc->n = 0;
    partitionstruc->limballoc = 0;
    partitionstruc->limbsize = 0;
    partitionstruc->partition = partition;
    partitionstruc->partitioninvert = partitioninvert;
}

void schur_number_partition_init(schur_number_partition_t *partitionstruc, mp_size_t limballoc) {
    /* Initialise partitionstruc, c'est-à-dire alloue l'espace nécessaire à chaque ensemble de la partition et met à jour les variables adéquates. */
    
    // partition et partitioninvert ont normalement déjà été allouées
    unsigned long p = partitionstruc->pmax;
    mp_limb_t **partition = partitionstruc->partition;
    mp_limb_t **partitioninvert = partitionstruc->partitioninvert;
    
    for (unsigned long j = 0; j < p; j++) {
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
        if (partitioninvert[j]) {
            free(partitioninvert[j]);
        }
    }
    free(partition);
    free(partitioninvert);
}

static inline void schur_number_translation2(mp_limb_t *r_set, const mp_limb_t *s_set, mp_size_t r_limbsize, mp_size_t s_limbsize, unsigned long n) {
    /* Calcule s_set + n de taille s_limbsize et place le résultat dans r_set de taille r_limbsize. */
    
    unsigned int shift = n % GMP_NUMB_BITS;
    if (!shift) {
        shift = GMP_NUMB_BITS - 1;
    }
    mp_limb_t limb_overflow = mpn_lshift(r_set, s_set, s_limbsize, shift);     // r_set = work + shift
    if (r_limbsize > s_limbsize) {
        r_set[s_limbsize] = limb_overflow;
    }
    n -= shift;
    
    while (n > 0) {
        shift = n % GMP_NUMB_BITS;
        if (!shift) {
            shift = GMP_NUMB_BITS - 1;
        }
        
        mpn_lshift(r_set, r_set, r_limbsize, shift);     // r_set = work + shift
        
        n -= shift;
    }
}

void schurNumberSumset(mp_limb_t *r_set, mp_limb_t *set1, mp_limb_t *set2, mp_size_t r_limbsize, mp_size_t limbsize, unsigned long x) {
    /* Cette fonction calcule set1 + set2 - x et le place dans r_set, qui doit pouvoir contenir tous les éléments. */
    unsigned long nsize = limbsize * GMP_NUMB_BITS;
    
    mp_limb_t *work = calloc(r_limbsize, sizeof(mp_limb_t));
    
    for (unsigned long n = 0; n < nsize; n++) {
        // Parcourir les éléments de set1
        if (GET_POINT(set1, n)) {
            if (n < x) {
                // Calculer set2 - (x - n)
                schur_number_ntranslation(work, set2, limbsize, x - n);
                
            } else {
                // Calculer set2 + (n - x)
                schur_number_translation2(work, set2, r_limbsize, limbsize, n - x);
            }
            // Ajouter work à r_set
            mpn_ior_n(r_set, r_set, work, r_limbsize);
        }
    }
    free(work);
}

void schurNumberWeakSumset(mp_limb_t *r_set, mp_limb_t *set1, mp_limb_t *set2, mp_size_t r_limbsize, mp_size_t limbsize, unsigned long x) {
    /* Cette fonction calcule set1 + set2 - x et le place dans r_set, qui doit pouvoir contenir tous les éléments. */
    unsigned long nsize = limbsize * GMP_NUMB_BITS;
    
    mp_limb_t *work = calloc(r_limbsize, sizeof(mp_limb_t));
    
    for (unsigned long n = 0; n < nsize; n++) {
        // Parcourir les éléments de set1
        if (GET_POINT(set1, n)) {
            if (n < x) {
                // Calculer set2 - (x - n)
                schur_number_ntranslation(work, set2, limbsize, x - n);
                
                if (2 * (x - n) < nsize) {
                    DELETE_POINT(work, 2 * (x - n));
                }
                
            } else {
                // Calculer set2 + (n - x)
                schur_number_translation2(work, set2, r_limbsize, limbsize, n - x);
                
                if (2 * (n - x) < nsize) {
                    DELETE_POINT(work, 2 * (n - x));
                }
                
            }
            // Ajouter work à r_set
            mpn_ior_n(r_set, r_set, work, r_limbsize);
        }
    }
    free(work);
}

void schur_number_set_revert(mp_limb_t *r_set, mp_limb_t *set, mp_size_t limbsize) {
    /* Cette fonction renverse set et place le résultat dans r_set. */
    for (mp_size_t i = 0; i < limbsize; i++) {
        mp_limb_t v = set[i];
        // swap odd and even bits
        v = ((v >> 1) & (mp_limb_t)0x5555555555555555) | ((v & (mp_limb_t)0x5555555555555555) << 1);
        // swap consecutive pairs
        v = ((v >> 2) & (mp_limb_t)0x3333333333333333) | ((v & (mp_limb_t)0x3333333333333333) << 2);
        // swap nibbles ...
        v = ((v >> 4) & (mp_limb_t)0x0F0F0F0F0F0F0F0F) | ((v & (mp_limb_t)0x0F0F0F0F0F0F0F0F) << 4);
        // swap bytes
        v = ((v >> 8) & (mp_limb_t)0x00FF00FF00FF00FF) | ((v & (mp_limb_t)0x00FF00FF00FF00FF) << 8);
        // swap 2-byte long pairs
    #if GMP_NUMB_BITS == 32
        v = ( v >> 16             ) | ( v               << 16);
    #else
        v = ((v >> 16) & (mp_limb_t)0x0000FFFF0000FFFF) | ((v & (mp_limb_t)0x0000FFFF0000FFFF) << 16);
        v = ( v >> 32                                 ) | ( v                                  << 32);
    #endif
        r_set[limbsize - i - 1] = v;
    }
    mpn_lshift(r_set, r_set, limbsize, 1);
}
