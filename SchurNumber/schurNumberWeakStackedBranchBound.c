//
//  schurNumberWeakStackedBranchBound.c
//  SchurNumber
//
//  Created by rubis on 07/04/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberMethods.h"
#include "gmp-impl.h"
#include "longlong.h"

unsigned long schurNumberWeakStackedBranchBound(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit) {
    /*
     Cette fonction calcule successivement les nombres de Schur S(p) pour p<= pmax, en partant de la partition initiale contenue dans partitionstruc.
     Elle se limite à explorer les partitions de taille <= nlimit.
     
     La partition est représentée comme un tableau de grands entiers sfpartition.
     Chaque grand entier représente un ensemble grâce à ses bits: si le bit k vaut 1,
     alors k appartient à l'ensemble; sinon il n'y appartient pas.
     
     Pour tester si il est possible d'ajouter n+1 à l'ensemble P sans-somme,
     on effectue un & entre les [n/2] premiers bits et les [n/2] derniers.
     Si le résultat est 0, il est possible d'ajouter n+1.
     
     Lorsqu'un entier n n'a pu être inséré dans aucune huches, on revient au sup sur les huches
     du plus petit n'≥[n/2]+1 tel que n' et n-n' appartiennent simultanément à cette huche.
     
     De plus, les sommes faibles de la partition sont conservées en mémoire dans le tableau sumpartition, de sorte qu'à chaque itération, il est possible de réaliser leur intersection pour borner la taille de la partition faiblement sans-somme en construction.
     */
    
    // Initialisation de la partition à calculer
    mp_limb_t **partition = partitionstruc->partition;
    mp_size_t limballoc = partitionstruc->limballoc;            // Nombre de limbes alloué à chaque ensemble de partition
    mp_size_t limbsize = partitionstruc->limbsize;              // Nombre de limbes utilisés par les ensembles de partition
    unsigned long nsize = mp_bits_per_limb * limbsize - 1;      // Plus grand entier pouvant être contenu dans limbsize limbes
    unsigned long nalloc = mp_bits_per_limb * limballoc - 1;    // Plus grand entier pouvant être contenu dans limballoc limbes
    
    // Initialisation des ensembles intermédiaires
    mp_limb_t *work1 = calloc(sizeof(mp_limb_t), limballoc);
    
    // Initialisation des variables
    const unsigned long n0 = partitionstruc->n;       // Taille de la partition initiale
    unsigned long n = n0;                       // Taille de l'intervalle
    unsigned long nbest = n0;                   // Taille de la plus grande partition trouvée
    unsigned long nblocking = 1;                // Plus petit entier ≥ [n/2]+1 tel que nblocking et (n+1) - nblocking soient dans la même huche
    unsigned long i = 0;                        // Huche où placer l'entier suivant
    unsigned long p = partitionstruc->p;        // Nombre de huches non vides
    unsigned long pmax = partitionstruc->pmax;  // Nombre de huches total
    char notSumFree = 1;
    char is_new_branch = 1;
    
    // Initialisation du tableau contenant la pile des sommes restreintes de la partition
    mp_limb_t **sumpartition = calloc(sizeof(mp_limb_t *), pmax);
    mp_limb_t **sums_ptr = calloc(sizeof(mp_limb_t *), pmax);
    
    unsigned long nmax = nalloc;    // Plus grand entier pouvant être atteint par cette partition
    
    unsigned long *cardinals = calloc(sizeof(unsigned long), pmax);     // Tableau contenant les cardinaux des ensembles de partition
    char *popcounts = calloc(limballoc, pmax);                          // Tableau dont l'élément j * limballoc + k compte le nombre de bit dans le limbe k de partition[j]
    
    for (unsigned long j = 0; j < pmax; j++) {
        sumpartition[j] = calloc(sizeof(mp_limb_t) * limballoc, nlimit);    // Tableau contenant la pile des sommes de l'ensembles j
        sums_ptr[j] = sumpartition[j];                                      // Pointeur vers le sommet de la pile des sommes de l'ensembles j
        schurNumberSumset(sums_ptr[j], partition[j], partition[j], limballoc, 0);
        
        unsigned long cardinal = 0;
        for (unsigned long k = 0; k < limballoc; k++) {
            mp_bitcnt_t poplimb;
            popc_limb(poplimb, *(partition[j] + k));
            popcounts[j * limballoc + k] = poplimb;
            cardinal += poplimb;
        }
        cardinals[j] = cardinal;
    }
    
    unsigned long iter_num = action->iter_num;  // Nombre d'itérations
    
    // Itération jusqu'à énumérer toutes les partions sans-somme à au plus pmax huches
    while (n >= n0) {
        // Placer n+1 dans une des huches en conservant le caractère faiblement sans-somme
        
        if (p == pmax) {
            // Effectuer l'intersection des sommes pour obtenir un majorant sur la plus grande partition sans-somme
            mpn_copyd(work1, *sums_ptr, limballoc);
            for (unsigned long j = 1; j < p; j++) {
                mpn_and_n(work1, work1, sums_ptr[j], limballoc);
            }
            nmax = nalloc;
            if (!mpn_zero_p(work1, limballoc)) {
                nmax = mpn_scan1(work1, limballoc);
            }
        }
        
        while (notSumFree && i < p) {
            // Regarder si n+1 appartient à la somme de l'ensemble i
            
            iter_num++;
            notSumFree = !!GET_POINT(sums_ptr[i], n+1);
            
            if (notSumFree) {
                // Déterminer nblocking en cas de blocage
                // Il s'agit du plus petit m tel que la somme de l'ensemble i à l'étape m rencontre n
                
                mp_limb_t *blockingsumset = sums_ptr[i] - limballoc;
                
                unsigned long c = cardinals[i] - 1;     // Numéro de l'entier m dans l'ensemble i, commençant à 0
                
                while (c > 0 && GET_POINT(blockingsumset, n+1)) {
                    c--;
                    blockingsumset -= limballoc;
                }
                
                // Déterminer de façon décroissante le limbe contenant l'entier bloquant
                mp_limb_t *limb_ptr = partition[i] + limbsize - 1;
                char *count_ptr = &(popcounts[i * limballoc + limbsize - 1]);
                unsigned long m = cardinals[i] - *count_ptr;
                while (m > c) {
                    limb_ptr --;
                    count_ptr --;
                    m -= *count_ptr;
                }
                mp_limb_t blockinglimb = *limb_ptr;
                
                // Déterminer la position au sein du limbe de l'entier bloquant
                while (m < c) {
                    blockinglimb &= blockinglimb - 1;
                    m++;
                }
                count_trailing_zeros(m, blockinglimb);
                m += (limb_ptr - partition[i]) * GMP_NUMB_BITS;
                
                if (m > nblocking) {
                    nblocking = m;
                }
            }
            
            i += notSumFree;
        }
        
        if (notSumFree || n >= nlimit || nmax <= nbest) {
            // L'entier n+1 n'a pu être placé dans aucune huche
            
            if (p < pmax && i <= p && n < nlimit) {
                // Remplir une nouvelle huche
                
                n++;
                ADD_POINT(partition[p], n);
                cardinals[p] = 1;
                popcounts[p * limballoc + limbsize - 1] = 1;

                sums_ptr[p] = sumpartition[p];
                
                i = 0;
                if (n > nsize) {
                    limbsize++;
                    nsize += mp_bits_per_limb;
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
                }
                is_new_branch = 0;
                
                // Déterminer la huche contenant nblocking
                unsigned long iblocking = 0;
                while (iblocking < p && !GET_POINT(partition[iblocking], nblocking)) {
                    iblocking++;
                }
                
                // Revenir à une partition de [1, nblocking-1]
                
                // Masquer tous les bits au-delà de nblocking en construisant un masque
                mpn_zero(work1, limbsize);
                
                mp_size_t blockinglimbsize = (nblocking / mp_bits_per_limb) + 1;    // Nombre de limbes nécessaires pour contenir nblocking
                
                *work1 = (mp_limb_t)1;
                
                mpn_neg(work1, work1, blockinglimbsize);                            // Attribue 1 à tous les bits < nblockinglimbsize * 64
                
                work1[blockinglimbsize - 1] >>= mp_bits_per_limb - (nblocking % mp_bits_per_limb);
                
                // Appliquer le masque
                for (unsigned long i = 0; i < p; i++) {
                    mpn_and_n(partition[i], work1, partition[i], limbsize);
                }
                
                // Dépiler les sommes
                for (unsigned long i = 0; i < p; i++) {
                    mp_bitcnt_t cardinal = mpn_popcount(partition[i], blockinglimbsize);
                    
                    if (cardinal < cardinals[i]) {
                        mp_bitcnt_t limb_diff = (cardinals[i] - cardinal - 1) * limballoc;
                        sums_ptr[i] -= limb_diff;
                        mpn_zero(sums_ptr[i], limb_diff + limballoc);
                        sums_ptr[i] -= (!!cardinal) * limballoc;
                        
                        cardinals[i] = cardinal;
                        
                        char popcnt;
                        popc_limb(popcnt, partition[i][blockinglimbsize - 1]);
                        popcounts[i * limballoc + blockinglimbsize - 1] = popcnt;
                        
                        for (unsigned long k = blockinglimbsize; k < limbsize; k++) {
                            popcounts[i * limballoc + k] = 0;
                        }
                        
                    }
                    
                }
                
                while (0 < p && mpn_zero_p(partition[p - 1], limbsize)) {
                    // Supprimer la dernière huche
                    p--;
                }
                
                n = nblocking - 1;
                limbsize = blockinglimbsize;
                nsize = mp_bits_per_limb * blockinglimbsize - 1;
                
                i = iblocking + 1;
                nblocking = n;
            }
        } else {
            // Adjoindre n+1 à la huche i et incrémenter n
            n++;
            
            mp_limb_t *sumset = sums_ptr[i] + limballoc;
            schur_number_translation(sumset, partition[i], limbsize, n);
            mpn_ior_n(sumset, sumset, sums_ptr[i], limballoc);
            sums_ptr[i] = sumset;
            
            ADD_POINT(partition[i], n);
            
            cardinals[i]++;
            popcounts[i * limballoc + limbsize - 1] ++;
            
            i = 0;
            if (n > nsize) {
                limbsize++;
                nsize += mp_bits_per_limb;
            }
            is_new_branch = 1;
            nblocking = n;
        }
        notSumFree = 1;
    }
    
    // Nettoyage
    free(work1);
    
    free(cardinals);
    free(popcounts);
    
    for (unsigned long j = 0; j < pmax; j++) {
        free(sumpartition[j]);
    }
    free(sumpartition);
    free(sums_ptr);
    
    action->iter_num = iter_num;
    
    return nbest;
}
