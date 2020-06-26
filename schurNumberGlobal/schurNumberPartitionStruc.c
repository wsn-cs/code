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
    
    mp_limb_t *memory = calloc(p, 2 * limballoc * sizeof(mp_limb_t));
    for (unsigned long j = 0; j < p; j++) {
        partitioninvert[j] = memory + 2 * j * limballoc;
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
    mp_limb_t **partition = partitionstruc->partition;
    mp_limb_t **partitioninvert = partitionstruc->partitioninvert;
    
    if (*partitioninvert) {
        free(*partitioninvert);
    }
    
    free(partition);
    free(partitioninvert);
}

void schur_number_partition_set_empty(schur_number_partition_t *partitionstruc) {
    /* Rend vide tous les ensembles de la partition, puis adjoint 0 au premier ensemble. */
    unsigned long p = partitionstruc->pmax;
    mp_size_t limballoc = partitionstruc->limballoc;
    mp_limb_t **partitioninvert = partitionstruc->partitioninvert;
    mp_limb_t **partition = partitionstruc->partition;
    
    mpn_zero(*partitioninvert, 2 * p * limballoc);
    ADD_POINT(partition[0], 0);
    
    partitionstruc->p = 0;
    partitionstruc->n = 0;
    partitionstruc->limbsize = 1;
}


void schur_number_sumset(mp_limb_t *r_set, mp_limb_t *set1, mp_limb_t *set2, mp_size_t r_limbsize, mp_size_t limbsize, unsigned long x, mp_limb_t *work) {
    /* Cette fonction calcule set1 + set2 - x et le place dans r_set, qui doit pouvoir contenir tous les éléments.
     Le pointeur work doit pointer vers un tableau pouvant contenir r_limbsize. */
    unsigned long nsize = limbsize * GMP_NUMB_BITS;
    
    if (x) {
        mpn_copyd(work, set2, limbsize);
        mpn_zero(work + limbsize, r_limbsize - limbsize);
        
        unsigned long nlast = x + 1;    // Dernier point de set1 trouvé, de sorte que work = set2 - (x - nlast + 1)
        
        // Parcourir d'abord les éléments de set1 strictement inférieurs à x
        for (unsigned long n = 0; n < x; n++) {
            if (GET_POINT(set1, n)) {
                // Calculer set2 - (x - n + 1) = set2 - (x - nlast + 1) - (nlast - n) = work - (nlast - n)
                schur_number_ntranslation(work, work, r_limbsize, nlast - n);
                nlast = n;
                
                // Ajouter work à r_set
                mpn_ior_n(r_set, r_set, work, r_limbsize);
            }
        }
    }
    
    mpn_copyd(work, set2, limbsize);
    mpn_zero(work + limbsize, r_limbsize - limbsize);
    
    unsigned long nlast = x;    // Dernier point de set1 trouvé, de sorte que work = set2 + (nlast - x)
    
    // Parcourir ensuite les éléments de set1 supérieurs à x
    for (unsigned long n = x; n < nsize; n++) {
        if (GET_POINT(set1, n)) {
            // Calculer set2 + (n - x) = set2 + (n - nlast) + (nlast - x) = work + (nlast - x)
            schur_number_translation(work, work, r_limbsize, n - nlast);
            nlast = n;
            
            // Ajouter work à r_set
            mpn_ior_n(r_set, r_set, work, r_limbsize);
        }
    }
}


void schur_number_restricted_sumset(mp_limb_t *r_set, mp_limb_t *set1, mp_limb_t *set2, mp_size_t r_limbsize, mp_size_t limbsize, unsigned long x, mp_limb_t *work) {
    /* Cette fonction calcule set1 + set2 - x et le place dans r_set, qui doit pouvoir contenir tous les éléments.
     Le pointeur work doit pointer vers un tableau pouvant contenir r_limbsize. */
    unsigned long nsize = limbsize * GMP_NUMB_BITS;
    unsigned long r_nsize = r_limbsize * GMP_NUMB_BITS;
    
    if (x) {
        mpn_copyd(work, set2, limbsize);
        mpn_zero(work + limbsize, r_limbsize - limbsize);
        
        unsigned long nlast = x + 1;    // Dernier point de set1 trouvé, de sorte que work = set2 - (x - nlast + 1)
        
        // Parcourir les éléments de set1 strictement inférieurs à x
        for (unsigned long n = x; n > 0; n--) {
            if (GET_POINT(set1, n)) {
                // Calculer set2 - (x - n + 1) = set2 - (x - nlast + 1) - (nlast - n) = work - (nlast - n)
                schur_number_ntranslation(work, work, r_limbsize, nlast - n);
                nlast = n;
                
                if (2 * (x - n) < nsize) {
                    DELETE_POINT(work, 2 * (x - n));
                }
                // Ajouter work à r_set
                mpn_ior_n(r_set, r_set, work, r_limbsize);
            }
        }
    }
    
    mpn_copyd(work, set2, limbsize);
    mpn_zero(work + limbsize, r_limbsize - limbsize);
    
    unsigned long nlast = x;    // Dernier point de set1 trouvé, de sorte que work = set2 + (nlast - x)
    
    // Parcourir les éléments de set1 supérieurs à x
    for (unsigned long n = x; n < nsize; n++) {
        if (GET_POINT(set1, n)) {
            // Calculer set2 + (n - x) = set2 + (n - nlast) + (nlast - x) = work + (nlast - x)
            schur_number_translation(work, work, r_limbsize, n - nlast);
            nlast = n;
            
            if (2 * (n - x) < r_nsize) {
                DELETE_POINT(work, 2 * (n - x));
            }
            // Ajouter work à r_set
            mpn_ior_n(r_set, r_set, work, r_limbsize);
        }
    }
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
