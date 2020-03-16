//
//  schurNumberWeakInterval.c
//  schurNumberPuncturedInterval
//
//  Created by rubis on 16/03/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberWeakInterval.h"

void schur_number_translation(mp_limb_t *r_set, mp_limb_t *set, mp_size_t limbsize, unsigned long n) {
    /* Calcule set + n et place le résultat dans r_set. */
    mp_limb_t *work = set;
    
    while (n > 0) {
        unsigned int shift = n % mp_bits_per_limb;
        if (!shift) {
            shift = mp_bits_per_limb - 1;
        }
        
        mpn_lshift(r_set, work, limbsize, shift);     // r_set = work + shift
        work = r_set;
        
        n -= shift;
    }
}

void schur_number_ntranslation(mp_limb_t *r_set, mp_limb_t *set, mp_size_t limbsize, unsigned long n) {
    /* Calcule set - n et place le résultat dans r_set. */
    mp_limb_t *work = set;

    while (n > 0) {
        unsigned int shift = n % mp_bits_per_limb;
        if (!shift) {
            shift = mp_bits_per_limb - 1;
        }
        
        mpn_rshift(r_set, work, limbsize, shift);     // r_set = work - shift
        work = r_set;
        
        n -= shift;
    }
}

void schurNumberRestrictedSumset(mp_limb_t *r_set, mp_limb_t *set, mp_size_t limbsize, unsigned long x) {
    /* Cette fonction calcule set +¨ set - x et le place dans r_set, qui doit pouvoir contenir tous les éléments. */
    unsigned long nsize = limbsize * mp_bits_per_limb;
    
    mp_limb_t *work1 = calloc(sizeof(mp_limb_t), 2 * limbsize);
    mp_limb_t *work2 = calloc(sizeof(mp_limb_t), limbsize);
    
    //schurNumberPrintSet(nsize, set);
    
    for (unsigned long n = 0; n < nsize; n++) {
        // Parcourir les éléments de set
        if (GET_POINT(set, n)) {
            if (n < x) {
                // Calculer set - (x - n)
                mpn_copyd(work2, set, limbsize);
                DELETE_POINT(work2, n);
                
                schur_number_ntranslation(work1, work2, limbsize, x - n);
                
            } else {
                // Calculer set + (n - x)
                mpn_copyd(work2, set, limbsize);
                DELETE_POINT(work2, n);
                
                schur_number_translation(work1, work2, limbsize, n - x);
            }
            // Ajouter work1 à r_set
            //schurNumberPrintSet(nsize, work1);
            mpn_ior_n(r_set, r_set, work1, 2 * limbsize);
        }
    }
    free(work1);
    free(work2);
}

unsigned long schurNumberPuncturedInterval(schur_number_partition_t *partitionstruc, struct schurNumberIOAction *action, unsigned long nlimit, mp_limb_t **constraint_partition, mp_size_t constraint_size, unsigned long depth) {
    /* Cette fonction génère des intervalles percés à partir de la partition de contrainte fournie. Elle lui adjoint un nouvel ensemble et tente de
     l'agrandir jusqu'à ce que les sommes restreintes soient telles qu'elles forment un intervalle de longueur très grande.*/
    
    // Initialisation de la partition à calculer
    mp_limb_t **partition = partitionstruc->partition;
    mp_limb_t **partitioninvert = partitionstruc->partitioninvert;
    mp_size_t limballoc = partitionstruc->limballoc;        // Nombre de limbes alloué à chaque ensemble de partition
    mp_size_t limbsize = partitionstruc->limbsize;          // Nombre de limbes utilisés par les ensembles de partition
    unsigned long nsize = mp_bits_per_limb * limbsize;      // Plus grand entier pouvant être contenu dans limbsize limbes
    unsigned long nalloc = mp_bits_per_limb * limballoc;    // Plus grand entier pouvant être contenu dans limballoc limbes
    
    mp_limb_t *work1 = calloc(sizeof(mp_limb_t), 2 * limballoc);
    mp_limb_t *work2 = calloc(sizeof(mp_limb_t), constraint_size);
    
    // Initialisation des variables
    unsigned long n0 = partitionstruc->n;   // Taille de la partition initiale
    unsigned long n = n0 + 1;               // Taille de l'intervalle+1
    unsigned long nbest = n0;               // Taille de la plus grande partition trouvée
    unsigned long p = partitionstruc->pmax; // Nombre de huches total
    unsigned long i = p;                    // Huche où placer l'entier suivant
    unsigned char notSumFree = 1;
    unsigned char is_new_branch = 1;
    
    // Initialisation du tableau contenant la pile des sommes restreintes de l'ensemble p
    mp_size_t sumlimballoc = 2 * limballoc;     // Nombre de limbes alloués à chaque somme de l'ensemble p
    mp_limb_t *sumsets = calloc(sizeof(mp_limb_t) * sumlimballoc, 3 * depth);   // Tableau contenant la pile des sommes des ensembles p
    mp_limb_t *sumset = sumsets;                // Somme restreinte de l'ensemble p, sommet de la pile sumsets
    //unsigned long *nlimits = calloc(sizeof(unsigned long), 3 * d);
    //*nlimits = nlimit;
    
    // Construction des sommes de constraint_partition
    mp_limb_t **constraint_partition_sum = calloc(sizeof(mp_limb_t *), p - 1);
    
    for (unsigned long j = 0; j < p - 1; j++) {
        mp_limb_t *work0 = calloc(sizeof(mp_limb_t), 2 * constraint_size);
        schurNumberRestrictedSumset(work0, constraint_partition[j], constraint_size, depth);
        constraint_partition_sum[j] = work0;
    }
    
    unsigned long iter_num = action->iter_num;  // Nombre d'itérations
    
    while (n > n0) {
        while (notSumFree && i > 0) {
            // Regarder si n peut être adjoint à l'ensemble i
            iter_num ++;
            i -= notSumFree;
            
            if (i == p-1) {
                // Regarder si n appartient à sumset

                if (n > depth) {
                    notSumFree = !!GET_POINT(sumset, n - depth);
                } else {
                    notSumFree = 0;
                }
                
            } else {
                // Regarder si n + d appartient à constraint_partition[i] + constraint_partition[i]
                notSumFree = !!GET_POINT(constraint_partition_sum[i], n);
                
                // Calcul de work1 = (n+1) - partition[i] = partitioninvert[i] - (nsize - n)
                schur_number_ntranslation(work1, partitioninvert[i] + (limballoc - limbsize), limbsize, nsize - n);
                
                // Intersection
                mpn_and_n(work2, work1, constraint_partition[i], constraint_size);
                notSumFree |= !mpn_zero_p(work2, constraint_size);
            }
        }
        
        if (notSumFree || n > nlimit) {
            /*n ne peut être placé nul part*/
            if (is_new_branch) {
                /*Ecrire la partition*/
                action->func(partition, n-1, action);
                if (n > nbest) {
                    nbest = n;
                }
            }
            is_new_branch = 0;
            
            /*Retirer n-1 de la partition*/
            n--;
            i = 0;
            while (i < p && !GET_POINT(partition[i], n)) {
                i++;
            }
            
            if (i == p-1) {
                sumset -= sumlimballoc;
                //nlimit = nlimits[n];
            }
            
            DELETE_POINT(partition[i], n);
            DELETE_POINT(partitioninvert[i], nalloc - n);
            if (nsize >= n + mp_bits_per_limb) {
                limbsize--;
                nsize -= mp_bits_per_limb;
            }
            
        } else {
            /*n est placé dans l'ensemble i*/
            
            if (i == p - 1) {
                // Calculer partition[i] + n + 2d et la placer dans la pile sumsets
                schur_number_translation(work1, partition[i], limbsize, n);
                
                // Réaliser l'union de sumset avec work1
                mpn_ior_n(sumset + sumlimballoc, sumset, work1, sumlimballoc);
                sumset += sumlimballoc;
                
                // Déterminer la longueur du plus grand intervalle de sumset
                /*mpn_lshift(work1, sumset, sumlimballoc, 1);
                mpn_and_n(work1, sumset, work1, sumlimballoc);  // Les bits mis correspondent aux bornes des intervalles
                
                mp_bitcnt_t pos2 = 0;
                mp_bitcnt_t pos_max = 0;
                mp_bitcnt_t len_max = 0;
                
                while (!mpn_zero_p(work1, sumlimballoc)) {
                    mp_bitcnt_t pos1 = mpn_scan1(work1, pos2);
                    DELETE_POINT(work1, pos1);
                    pos2 = mpn_scan1(work1, pos1);
                    DELETE_POINT(work1, pos2);
                    mp_bitcnt_t len = pos2 - pos1;
                    if (len > len_max) {
                        len_max = len;
                        pos_max = pos1;
                    }
                }
                
                if (len_max >= d) {
                    nlimit = pos_max + d;
                    nlimits[n] = nlimit;
                } else {
                    nlimits[n] = nlimit;
                }*/
            }
            
            ADD_POINT(partition[i], n);
            ADD_POINT(partitioninvert[i], nalloc - n);
            n++;
            if (nsize < n) {
                limbsize++;
                nsize += mp_bits_per_limb;
            }
            i = p;
            is_new_branch = 1;
        }
        notSumFree = 1;
    }
    
    // Nettoyage
    free(work1);
    free(work2);
    
    free(sumsets);
    
    for (unsigned long j = 0; j < p-1; j++) {
        free(constraint_partition_sum[j]);
    }
    free(constraint_partition_sum);
    
    action->iter_num = iter_num;
    
    return nbest - 1;
}
