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

void schurNumberSumset(mp_limb_t *r_set, mp_limb_t *set1, mp_limb_t *set2, mp_size_t limbsize, unsigned long x) {
    /* Cette fonction calcule set1 + set2 - x et le place dans r_set, qui doit pouvoir contenir tous les éléments. */
    unsigned long nsize = limbsize * mp_bits_per_limb;
    
    mp_limb_t *work = calloc(limbsize, sizeof(mp_limb_t));
    
    for (unsigned long n = 0; n < nsize; n++) {
        // Parcourir les éléments de set1
        if (GET_POINT(set1, n)) {
            if (n < x) {
                // Calculer set2 - (x - n)
                schur_number_ntranslation(work, set2, limbsize, x - n);
                
            } else {
                // Calculer set2 + (n - x)
                schur_number_translation(work, set2, limbsize, n - x);
            }
            // Ajouter work1 à r_set
            mpn_ior_n(r_set, r_set, work, limbsize);
        }
    }
    free(work);
}

void schurNumberWeakSumset(mp_limb_t *r_set, mp_limb_t *set1, mp_limb_t *set2, mp_size_t limbsize, unsigned long x) {
    /* Cette fonction calcule set1 + set2 - x et le place dans r_set, qui doit pouvoir contenir tous les éléments. */
    unsigned long nsize = limbsize * mp_bits_per_limb;
    
    mp_limb_t *work = calloc(limbsize, sizeof(mp_limb_t));
    
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
                schur_number_translation(work, set2, limbsize, n - x);
                
                if (2 * (n - x) < nsize) {
                    DELETE_POINT(work, 2 * (n - x));
                }
                
            }
            // Ajouter work1 à r_set
            mpn_ior_n(r_set, r_set, work, limbsize);
        }
    }
    free(work);
}
