//
//  schurNumberWeakBranchBound.c
//  SchurNumber
//
//  Created by rubis on 16/02/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberMethods.h"

unsigned long schur_number_weak_branch_bound(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit) {
    /*
     Cette fonction calcule successivement les nombres de Schur faibles WS(p) pour p<= pmax, en partant de la partition initiale contenue dans partitionstruc.
     Elle se limite à explorer les partitions de taille <= nlimit.
     
     La partition est représentée comme un tableau de grands entiers sfpartition.
     Chaque grand entier représente un ensemble grâce à ses bits: si le bit k vaut 1,
     alors k appartient à l'ensemble; sinon il n'y appartient pas.
     
     Pour tester si il est possible d'ajouter n+1 à l'ensemble P sans-somme,
     on effectue un & entre les [n/2] premiers bits et les [n/2] derniers.
     Si le résultat est 0, il est possible d'ajouter n+1.
     
     Lorsqu'un entier n n'a pu être inséré dans aucune huches, on revient au sup sur les huches
     du plus petit n'≥[n/2]+1 tel que n' et n-n' appartiennent simultanément à cette huche,
     c'est-à-dire le plus petit couple (n', n - n') empêchant le placement de n dans un des ensembles.
     */
    
    // Initialisation de la partition à calculer
    mp_limb_t **partition = partitionstruc->partition;
    mp_limb_t **partitioninvert = partitionstruc->partitioninvert;
    mp_size_t limballoc = partitionstruc->limballoc;    // Nombre de limbes alloué à chaque ensemble de sfpartition
    mp_size_t limbsize = partitionstruc->limbsize;      // Nombre de limbes utilisés par les ensembles de sfpartition
    unsigned long nsize = GMP_NUMB_BITS * limbsize - 1;      // Plus grand entier pouvant être contenu dans limbsize limbes
    unsigned long nalloc = GMP_NUMB_BITS * limballoc - 1;    // Plus grand entier pouvant être contenu dans limballoc limbes
    
    // Initialisation des ensembles intermédiaires
    mp_limb_t *work1 = calloc(sizeof(mp_limb_t), 2 * limballoc);
    mp_limb_t *work2 = work1 + limballoc;
    
    // Initialisation des variables
    unsigned long n0 = partitionstruc->n;   // Taille de la partition initiale
    unsigned long n = n0;                   // Taille de l'intervalle
    unsigned long nbest = n0;               // Taille de la plus grande partition trouvée
    unsigned long nblocking = 1;            // Plus petit entier ≥ [n/2]+1 tel que nblocking et (n+1) - nblocking soient dans la même huche
    unsigned long i = 0;                    // Huche où placer l'entier suivant
    unsigned long p = partitionstruc->p;    // Nombre de huches non vides
    unsigned long pmax = partitionstruc->pmax; // Nombre de huches total
    char notSumFree = 1;
    char is_new_branch = 1;
    
    unsigned long iter_num = action->iter_num;  // Nombre d'itérations
    action->limballoc = limballoc;
    
    // Itération jusqu'à énumérer toutes les partions sans-somme à au plus pmax huches
    while (n >= n0) {
        //schurNumberPrintPartition(p, n, partition);
        // Placer n+1 dans une des huches en conservant le caractère faiblement sans-somme
        while (notSumFree && i < p) {
            // Tester si l'ensemble obtenu en ajoutant n+1 à la huche i est faiblement sans-somme
            iter_num++;
            mpn_copyi(work1, partitioninvert[i] + (limballoc - limbsize), limbsize);
            
            // Retirer éventuellement (n + 1)/2
            if (n % 2 && GET_POINT(work1, nsize - (n>>1))) {
                DELETE_POINT(work1, nsize - (n >> 1));
            }
            
            // Calculer (n+1) - huche i = (nsize + 1 - partition[i]) - (nsize - n) = partitioninvert[i] - (nsize - n) en effectuant une succession de décalage vers la droite
            schur_number_ntranslation(work2, work1, limbsize, nsize - n);
            
            // Intersecter la somme avec la huche initiale
            mpn_and_n(work1, work2, partition[i], limbsize);
            notSumFree = !mpn_zero_p(work1, limbsize);
            
            if (notSumFree) {
                // Déterminer nblocking en cas de blocage
                unsigned long j = mpn_rscan1(work1, nsize) - 1;
                if (j > nblocking) {
                    nblocking = j;
                }
            }
            
            i += notSumFree;
        }
        
        if (notSumFree || n >= nlimit) {
            // L'entier n+1 n'a pu être placé dans aucune huche
            if (p < pmax && i <= p && n < nlimit) {
                // Remplir une nouvelle huche
                ADD_POINT(partitioninvert[p], nalloc - n);
                n++;
                ADD_POINT(partition[p], n);
                i = 0;
                if (n > nsize) {
                    limbsize++;
                    nsize += GMP_NUMB_BITS;
                }
                p++;
                is_new_branch = 1;
                nblocking = 1;
            } else {
                
                if (is_new_branch) {
                    // Agir sur la partition
                    action->func(partition, n, action);
                    if (n > nbest) {
                        nbest = n;
                    }
                    if (n >= nlimit) {
                        nblocking = n;
                    }
                    //schurNumberPrintPartition(p, n, partition);
                }
                is_new_branch = 0;
                
                /*if (!nblocking) {
                    nblocking = n;
                }*/
                
                // Déterminer la huche contenant nblocking
                unsigned long iblocking = 0;
                while (iblocking < p && !GET_POINT(partition[iblocking], nblocking)) {
                    iblocking++;
                }
                
                // Revenir à une partition de [1, nblocking-1] par intersection
                mp_size_t blockinglimbsize = ((nblocking - 1) >> 6) + 1;        // Nombre de limbes nécessaires pour contenir nblocking
                unsigned long nblockinginvert = nalloc - nblocking + 2;
                mp_size_t invert_blockinglimbsize = (nblockinginvert / GMP_NUMB_BITS) + 1;
                for (unsigned long i = 0; i < p; i++) {
                    schur_number_intersect_interval_n(partitioninvert[i], limballoc, invert_blockinglimbsize, nblockinginvert);
                    schur_number_intersect_interval_0(partition[i], limballoc, blockinglimbsize, nblocking - 1);
                }
                
                while (0 < p && mpn_zero_p(partition[p - 1], limbsize)) {
                    // Supprimer la dernière huche
                    p--;
                }
                
                n = nblocking - 1;
                limbsize = blockinglimbsize;
                nsize = GMP_NUMB_BITS * blockinglimbsize - 1;
                
                i = iblocking + 1;
                nblocking = n;
            }
        } else {
            // Adjoindre n+1 à la huche i et incrémenter n
            ADD_POINT(partitioninvert[i], nalloc - n);
            n++;
            ADD_POINT(partition[i], n);
            i = 0;
            if (n >= nsize) {
                limbsize++;
                nsize += GMP_NUMB_BITS;
            }
            is_new_branch = 1;
            nblocking = 1;
        }
        notSumFree = 1;
    }
    
    // Nettoyage
    free(work1);
    
    action->iter_num = iter_num;
    
    return nbest;
}
