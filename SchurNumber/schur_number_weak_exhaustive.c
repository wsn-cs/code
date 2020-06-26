//
//  schurNumberWeakExhaustive.c
//  SchurNumber
//
//  Created by rubis on 10/02/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberMethods.h"

unsigned long schur_number_weak_exhaustive(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit) {
    /*
     Cette fonction calcule successivement les nombres de Schur faibles WS(p) pour p<= pmax, en partant de la partition initiale contenue dans partitionstruc. Elle remplit le tableau nbests, et compte le nombre d'itérations dans iternum.
     
     La partition est représentée comme un tableau de grands entiers sfpartition.
     Chaque grand entier représente un ensemble grâce à ses bits: si le bit k vaut 1,
     alors k appartient à l'ensemble; sinon il n'y appartient pas.
     
     Pour tester si il est possible d'ajouter n+1 à l'ensemble P sans-somme,
     on effectue un & entre les [n/2] premiers bits et les [n/2] derniers.
     Si le résultat est 0, il est possible d'ajouter n+1.
     
     Lorsqu'un entier n n'a pu être inséré dans aucune huches, on retire le dernier entier ajouté.
     */
    
    // Initialisation de la partition à calculer
    mp_limb_t **wsfpartition = partitionstruc->partition;
    mp_limb_t **wsfpartitioninvert = partitionstruc->partitioninvert;
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
    unsigned long i = 0;                    // Huche où placer l'entier suivant
    unsigned long p = partitionstruc->p;    // Nombre de huches non vides
    unsigned long pmax = partitionstruc->pmax; // Nombre de huches total
    char notSumFree = 1;
    char is_new_branch = 1;
    
    unsigned long iter_num = action->iter_num;  // Nombre d'itérations
    
    // Itération jusqu'à énumérer toutes les partions sans-somme à au plus pmax huches
    while (n >= n0) {
        // Placer n+1 dans une des huches en conservant le caractère faiblement sans-somme
        while (notSumFree && i < p) {
            // Tester si l'ensemble obtenu en ajoutant n+1 à la huche i est faiblement sans-somme
            iter_num++;
            mpn_copyi(work1, wsfpartitioninvert[i] + (limballoc - limbsize), limbsize);
            
            // Retirer éventuellement (n + 1)/2
            if (n % 2 && GET_POINT(work1, nsize - (n>>1))) {
                DELETE_POINT(work1, nsize - (n >> 1));
            }
            
            // Calculer (n+1) - huche i = (nsize + 1 - wsfpartitioninvert[i]) - (nsize - n) en effectuant une succession de décalage vers la droite
            schur_number_ntranslation(work2, work1, limbsize, nsize - n);
            
            // Intersecter la somme avec la huche initiale
            mpn_and_n(work1, work2, wsfpartition[i], limbsize);
            notSumFree = !mpn_zero_p(work1, limbsize);
            
            i += notSumFree;
        }
        
        if (notSumFree || n >= nlimit) {
            // L'entier n+1 n'a pu être placé dans aucune huche
            if (p < pmax && i <= p && n < nlimit) {
                // Remplir une nouvelle huche
                ADD_POINT(wsfpartitioninvert[p], nalloc - n);
                n++;
                ADD_POINT(wsfpartition[p], n);
                i = 0;
                if (n > nsize) {
                    limbsize++;
                    nsize += GMP_NUMB_BITS;
                }
                p++;
                is_new_branch = 1;
            } else {
                
                if (is_new_branch) {
                    // Agir sur la partition
                    action->func(wsfpartition, n, action);
                    if (n > nbest) {
                        nbest = n;
                    }
                }
                is_new_branch = 0;
                
                // Déterminer la huche contenant n
                i = 0;
                while (!GET_POINT(wsfpartition[i], n)) {
                    i++;
                }
                // Retirer n de la huche puis décrémenter n
                DELETE_POINT(wsfpartition[i], n);
                
                if (i == (p-1) && mpn_zero_p(wsfpartition[i], limbsize)) {
                    // Supprimer la dernière huche
                    p--;
                }
                
                n--;
                DELETE_POINT(wsfpartitioninvert[i], nalloc - n);
                
                if (n + GMP_NUMB_BITS <= nsize) {
                    limbsize--;
                    nsize -= GMP_NUMB_BITS;
                }
                i++;
            }
            notSumFree = 1;
        } else {
            // Adjoindre n+1 à la huche i et incrémenter n
            ADD_POINT(wsfpartitioninvert[i], nalloc - n);
            n++;
            ADD_POINT(wsfpartition[i], n);
            i = 0;
            if (n > nsize) {
                limbsize++;
                nsize += GMP_NUMB_BITS;
            }
            is_new_branch = 1;
            notSumFree = 1;
        }
        
    }
    
    // Nettoyage
    free(work1);
    
    action->iter_num = iter_num;
    
    return nbest;
}
