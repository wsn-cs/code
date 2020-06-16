//
//  schurNumberWeakPuncturedInterval.c
//  schurNumberPuncturedInterval2
//
//  Created by rubis on 28/04/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberWeakPuncturedInterval.h"

unsigned long schur_number_weak_punctured_interval(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit) {
    /* Cette fonction essaie de combler l'intervalle [n0+1, nlimit] après avoir fixé une partition A de [1, n0] en respectant le caractère faiblement sans-somme. */
    
    // Initialisation de la partition à calculer
    mp_limb_t **partition = partitionstruc->partition;
    mp_limb_t **partitioninvert = partitionstruc->partitioninvert;  // Tableau contenant les ensembles (-A ⋃ A)
    mp_size_t limballoc = partitionstruc->limballoc;            // Nombre de limbes alloué à chaque ensemble de partition
    mp_size_t limballoc2 = 2 * partitionstruc->limballoc;       // Nombre de limbes alloué à chaque ensemble de partitioninvert
    unsigned long nalloc = mp_bits_per_limb * limballoc;        // Plus grand entier pouvant être contenu dans limballoc limbes
    
    // Initialisation des ensembles intermédiaires
    mp_limb_t *work1 = calloc(limballoc2, sizeof(mp_limb_t));
    
    // Initialisation des variables
    const unsigned long n0 = partitionstruc->n;       // Taille de la partition initiale
    unsigned long n = n0 + 1;                   // Taille de l'intervalle
    unsigned long nbest = n0;                   // Taille de la plus grande partition trouvée
    unsigned long p = partitionstruc->pmax;     // Nombre de huches
    unsigned long i = p;                        // Huche où placer l'entier suivant
    char notSumFree = 1;
    char is_new_branch = 1;
    
    //schurNumberPrintPartition(p, nalloc, partition);
    
    unsigned long iter_num = action->iter_num;  // Nombre d'itérations
    
    // Itération jusqu'à énumérer toutes les partitions sans-sommes de [n0 + 1, nlimit] à p huches
    while (n > n0) {
        // Placer n dans une des huches en conservant le caractère faiblement sans-somme
        while (notSumFree && i > 0) {
            // Tester si l'ensemble obtenu en ajoutant n à la huche i est faiblement sans-somme
            i--;
            iter_num++;
            
            // Calculer n + (-A ⋃ A)
            schur_number_translation(work1, partitioninvert[i], limballoc2, n);
            
            // Faire l'intersection avec A
            mp_limb_t *work0 = work1 + limballoc;   // Pointeur vers A ⋂ [0, +∞[
            
            if (!(n & (mp_limb_t)1)) {
                // n est pair
                DELETE_POINT(work0, n >> 1);
            }
            mpn_and_n(work0, work0, partition[i], limballoc);
            
            notSumFree = !mpn_zero_p(work0, limballoc);
            
            //i += notSumFree;
        }
        
        if (notSumFree || n > nlimit) {
            // L'entier n n'a pu être placé dans aucune huche
            
            if (is_new_branch) {
                // Agir sur la partition
                if (n > nlimit) {
                    //schurNumberPrintPartition(p, nalloc, partition);
                    action->func(partition, n-1, action);
                }
                if (n > nbest) {
                    nbest = n;
                }
            }
            is_new_branch = 0;
            
            // Déterminer la huche contenant n-1
            n--;
            i = p-1;
            while (!GET_POINT(partition[i], n)) {
                i--;
            }
            // Retirer n de la huche
            DELETE_POINT(partition[i], n);
            DELETE_POINT(partitioninvert[i], nalloc - n);
            
            //i++;
            
        } else {
            // Adjoindre n à la huche i puis incrémenter n
            ADD_POINT(partitioninvert[i], nalloc - n);
            ADD_POINT(partition[i], n);
            n++;
            i = p;
            is_new_branch = 1;
        }
        notSumFree = 1;
    }
    
    // Nettoyage
    free(work1);
    
    action->iter_num = iter_num;
    
    return nbest;
}
