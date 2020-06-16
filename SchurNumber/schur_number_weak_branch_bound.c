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
     du plus petit n'≥[n/2]+1 tel que n' et n-n' appartiennent simultanément à cette huche.
     */
    
    // Initialisation de la partition à calculer
    mp_limb_t **partition = partitionstruc->partition;
    mp_limb_t **partitioninvert = partitionstruc->partitioninvert;
    mp_size_t limballoc = partitionstruc->limballoc;    // Nombre de limbes alloué à chaque ensemble de sfpartition
    mp_size_t limbsize = partitionstruc->limbsize;      // Nombre de limbes utilisés par les ensembles de sfpartition
    unsigned long nsize = GMP_NUMB_BITS * limbsize - 1;      // Plus grand entier pouvant être contenu dans limbsize limbes
    unsigned long nalloc = GMP_NUMB_BITS * limballoc - 1;    // Plus grand entier pouvant être contenu dans limballoc limbes
    
    // Initialisation des ensembles intermédiaires
    mp_limb_t *work1 = calloc(sizeof(mp_limb_t), limballoc);
    mp_limb_t *work2 = calloc(sizeof(mp_limb_t), limballoc);
    
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
            
            // Calculer (n+1) - huche i = (nsize + 1 - wsfpartitioninvert[i]) - (nsize - n) en effectuant une succession de décalage vers la droite
            schur_number_ntranslation(work1, work1, limbsize, nsize - n);
            
            // Intersecter la somme avec la huche initiale
            mpn_and_n(work2, work1, partition[i], limbsize);
            notSumFree = !mpn_zero_p(work2, limbsize);
            
            if (notSumFree) {
                // Déterminer nblocking en cas de blocage
                /*printf("Ensembles : ");
                schurNumberPrintSet(n, partition[i]);
                printf("\n");
                schurNumberPrintSet(nsize, work2);
                printf("\n");*/
                unsigned long j = mpn_rscan1(work2, nsize) - 1;
                //printf("Blocage %lu pour %lu\n", j, n+1);
                /*if (j > nsize) {
                    unsigned long jj = mpn_rscan1(work2, nsize);
                    jj = mpn_rscan1(work2, nsize);
                    jj--;
                    printf("Blocage %lu pour %lu\n", jj, n+1);
                }*/
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
                
                // Revenir à une partition de [1, nblocking-1] grâce à deux masques
                
                mp_size_t blockinglimbsize = (n >> 6) + 1;          // Nombre de limbes nécessaires pour contenir nblocking
                
                schur_number_setinterval_1(work1, limbsize, blockinglimbsize, nblocking); // work1 = [1, nblocking-1]
                schur_number_setinterval_s(work2, limbsize, blockinglimbsize, nblocking); // work2 = [nsize - nblocking, nsize]
                
                /*printf("Masques : ");
                schurNumberPrintSet(STDOUT_FILENO, nsize, work1);
                printf("\n");
                schurNumberPrintSet(STDOUT_FILENO, nalloc, work2);
                printf("\n");*/
                
                //schurNumberPrintPartition(p, n, partition);
                //schurNumberPrintPartition(p, nalloc, partitioninvert);
                
                // Appliquer le masque
                for (unsigned long i = 0; i < p; i++) {
                    mpn_and_n(partition[i], work1, partition[i], limbsize);
                    mpn_and_n(partitioninvert[i] + (limballoc - limbsize), work2, partitioninvert[i] + (limballoc - limbsize), limbsize);
                }
                /*schurNumberPrintPartition(p, n, partition);
                schurNumberPrintPartition(p, nalloc, partitioninvert);*/
                
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
            if (n > nsize) {
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
    free(work2);
    
    action->iter_num = iter_num;
    
    return nbest;
}
