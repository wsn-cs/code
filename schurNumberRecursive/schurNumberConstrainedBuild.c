//
//  schurNumberConstrainedBuild.c
//  schurNumberRecursive
//
//  Created by rubis on 08/01/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberConstrainedBuild.h"

unsigned long schurNumberConstrainedBuild(schur_number_partition_t *partitionstruc, struct schurNumberIOAction *action, unsigned long nlimit, mp_limb_t **constraint_partition, mp_size_t constraint_size) {
    /*
     Cette fonction trouve toutes les partitions B tels que B - B ne rencontre pas A = constraint_partition, en partant de la partition initiale contenue dans partitionstruc. Chaque ensemble de constraint_partition possède constraint_size limbes.
     
     La partition B est représentée comme un tableau de grands entiers partition.
     Chaque grand entier représente un ensemble grâce à ses bits: si le bit k vaut 1,
     alors k appartient à l'ensemble; sinon il n'y appartient pas.
     
     Pour tester si il est possible d'ajouter n+1 à l'ensemble B sans-somme relativement à A,
     on effectue un & entre les premiers bits de A et les derniers bits de B décalés de n+1.
     Si le résultat est 0, il est possible d'ajouter n+1.
     
     Lorsqu'un entier n n'a pu être inséré dans aucune huches, on retire le dernier entier ajouté.
     */
    
    // Initialisation de la partition à calculer
    mp_limb_t **partition = partitionstruc->partition;
    mp_limb_t **partitioninvert = partitionstruc->partitioninvert;
    mp_size_t limballoc = partitionstruc->limballoc;    // Nombre de limbes alloué à chaque ensemble de sfpartition
    mp_size_t limbsize = partitionstruc->limbsize;      // Nombre de limbes utilisés par les ensembles de sfpartition
    unsigned long nsize = mp_bits_per_limb * limbsize;      // Plus grand entier pouvant être contenu dans limbsize limbes
    unsigned long nalloc = mp_bits_per_limb * limballoc;    // Plus grand entier pouvant être contenu dans limballoc limbes
    
    mp_limb_t *work1 = calloc(sizeof(mp_limb_t), limballoc);
    mp_limb_t *work2 = calloc(sizeof(mp_limb_t), constraint_size);
    
    // Initialisation des variables
    unsigned long n0 = partitionstruc->n;   // Taille de la partition initiale
    unsigned long n = n0 + 1;               // Taille de l'intervalle+1
    unsigned long nbest = n0;               // Taille de la plus grande partition trouvée
    unsigned long i = 0;                    // Huche où placer l'entier suivant
    unsigned long p = partitionstruc->pmax; // Nombre de huches total
    unsigned char notSumFree = 1;
    unsigned char is_new_branch = 1;
    
    unsigned long iter_num = action->iter_num;  // Nombre d'itérations
    
    while (n > n0) {
        while (notSumFree && i < p) {
            /*Regarder si n peut être adjoint à l'ensemble i*/
            iter_num ++;
            /*Calcul de work1 = (n+1) - partition[i] = partitioninvert[i] - (nsize - n).*/
            unsigned long nrem = nsize - n;
            mp_limb_t *work0 = partitioninvert[i] + (limballoc - limbsize);
            
            while (nrem > 0) {
                unsigned int shift = nrem % mp_bits_per_limb;
                if (!shift) {
                    shift = mp_bits_per_limb - 1;
                }
                
                mpn_rshift(work1, work0, limbsize, shift);     // work1 = work0 - shift
                work0 = work1;
                
                nrem -= shift;
            }
            
            /*Intersection*/
            mpn_and_n(work2, work1, constraint_partition[i], constraint_size);
            notSumFree = !mpn_zero_p(work2, constraint_size);
            i += notSumFree;
        }
        
        if (notSumFree || n > nlimit) {
            /*n ne peut être placé nul part*/
            if (is_new_branch) {
                /*Ecrire la partition*/
                nbest = action->func(partition, n-1, action);
            }
            is_new_branch = 0;
            
            /*Retirer n-1 de la partition*/
            n--;
            i = 0;
            while (i < p && !GET_POINT(partition[i], n)) {
                i++;
            }
            DELETE_POINT(partition[i], n);
            DELETE_POINT(partitioninvert[i], nalloc - n);
            if (nsize >= n + mp_bits_per_limb) {
                limbsize--;
                nsize -= mp_bits_per_limb;
            }
            i++;
            
        } else {
            /*n est placé dans l'ensemble i*/
            ADD_POINT(partition[i], n);
            ADD_POINT(partitioninvert[i], nalloc - n);
            n++;
            if (nsize < n) {
                limbsize++;
                nsize += mp_bits_per_limb;
            }
            i = 0;
            is_new_branch = 1;
        }
        notSumFree = 1;
    }
    
    // Nettoyage
    free(work1);
    free(work2);
    
    action->iter_num = iter_num;
    
    return nbest;
}
