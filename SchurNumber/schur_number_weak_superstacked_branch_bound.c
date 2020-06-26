//
//  schur_number_weak_superstacked_branch_bound.c
//  SchurNumber
//
//  Created by rubis on 21/06/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberMethods.h"
#include "gmp-impl.h"
#include "longlong.h"

static inline unsigned long find_nblocking(unsigned long x, mp_limb_t *sumset, mp_size_t limbsize, unsigned long count) {
    /* Renvoie le plus petit indice idx tel que x appartienne à sumset[idx], sachant que sumset est un tableau de count grands entiers de taille limbsize. */
    
    if (GET_POINT(sumset, x)) {
        return 0;
    }
    
    unsigned long i0 = 0;
    unsigned long i1 = count - 1;
    
    while (i1 - i0 > 1) {
        unsigned long imedian = (i1 + i0) >> 1;
        
        if (GET_POINT((sumset + imedian * limbsize), x)) {
            i1 = imedian;
        } else {
            i0 = imedian;
        }
    }
    
    return i1;
}

unsigned long schur_number_weak_superstacked_branch_bound(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit) {
    /*
     Cette fonction calcule successivement les nombres de Schur S(p) pour p<= pmax, en partant de la partition initiale contenue dans partitionstruc.
     Elle se limite à explorer les partitions de taille <= nlimit.
     
     La partition est représentée comme un tableau de grands entiers sfpartition.
     Chaque grand entier représente un ensemble grâce à ses bits: si le bit k vaut 1,
     alors k appartient à l'ensemble; sinon il n'y appartient pas.
     
     Pour tester si il est possible d'ajouter n+1 à l'ensemble P sans-somme,
     on regarde si n+1 appartient à P+P, qui est mis à jour à chaque ajout d'un nouvel entier.
     Les valeurs successives de P+P sont suavegardées dans une pile.
     
     Lorsqu'un entier n n'a pu être inséré dans aucune huches, on revient au sup sur les huches
     du plus petit n'≥[n/2]+1 tel que n' et n-n' appartiennent simultanément à cette huche.
     
     De plus, les sommes de la partition sont conservées en mémoire dans le tableau sumpartition, de sorte qu'à chaque itération, il est possible de réaliser leur intersection pour borner la taille de la partition sans-somme en construction.
     */
    
    // Initialisation de la partition à calculer
    mp_limb_t **partition = partitionstruc->partition;
    mp_limb_t **partitioninvert = partitionstruc->partitioninvert;
    mp_size_t limballoc = partitionstruc->limballoc;            // Nombre de limbes alloué à chaque ensemble de partition
    mp_size_t limbsize = partitionstruc->limbsize;              // Nombre de limbes utilisés par les ensembles de partition
    unsigned long nsize = GMP_NUMB_BITS * limbsize - 1;      // Plus grand entier pouvant être contenu dans limbsize limbes
    unsigned long nalloc = GMP_NUMB_BITS * limballoc - 1;    // Plus grand entier pouvant être contenu dans limballoc limbes
    
    // Initialisation des ensembles intermédiaires
    mp_limb_t *work1 = calloc(2, limballoc * sizeof(mp_limb_t));
    mp_limb_t *work2 = work1 + limballoc;
    
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
    
    unsigned long *cardinals = calloc(sizeof(unsigned long), pmax);     // Tableau contenant les cardinaux des ensembles de partition
    
    for (unsigned long j = 0; j < pmax; j++) {
        sumpartition[j] = calloc(sizeof(mp_limb_t) * limballoc, nlimit);    // Tableau contenant la pile des sommes de l'ensembles j
        sums_ptr[j] = sumpartition[j];                                      // Pointeur vers le sommet de la pile des sommes de l'ensembles j
        schur_number_restricted_sumset(sums_ptr[j], partition[j], partition[j], limballoc, limballoc, 0, work1);
        
        cardinals[j] = mpn_popcount(partition[j], limballoc);
    }
    
    unsigned long *setmin = calloc(pmax, sizeof(unsigned long));    // Tableau contenant les plus petits éléments de chaque partie
    unsigned long sum_min = 0;                                      // Somme des plus petits éléments, utile pour avoir une borne inférieure de nbest
    for (unsigned long j = 0; j < p; j++) {
        setmin[j] = mpn_scan1(partition[j], 1);
        sum_min += setmin[j];
    }
    //nbest = (n0 << (pmax - p)) + sum_min;
    if (pmax > 2) {
        unsigned long pow3 = 1;
        unsigned long pdiff = pmax - 2;
        while (pdiff > 0) {
            pow3 *= 3;
            pdiff--;
        }
        nbest = action->func(NULL, 7 * pow3 + pmax - 1, action);
    }
    //nbest = action->func(NULL, (n << (pmax - p)) + sum_min, action);
    
    // Initialisation du tableau contenant la pile des co-sommes Cj = intersection des Ak + Ak, k ≠ j
    // Pour plus d'efficacité, la co-somme est scindée en deux parties = les k < j  et les k > j
    mp_limb_t **cosumpartition_inf = calloc(sizeof(mp_limb_t *), 2 * pmax);
    mp_limb_t **cosumpartition_sup = cosumpartition_inf + pmax;
    mp_limb_t **cosums_inf_ptr = calloc(sizeof(mp_limb_t *), 2 * pmax);
    mp_limb_t **cosums_sup_ptr = cosums_inf_ptr + pmax;
    
    // Construction initale des co-sommes inférieures
    cosumpartition_inf[0] = calloc(sizeof(mp_limb_t) * limballoc, 1);
    cosums_inf_ptr[0] = cosumpartition_inf[0];
    mpn_com(cosums_inf_ptr[0], cosums_inf_ptr[0], limballoc);
    
    for (unsigned long j = 1; j < pmax; j++) {
        cosumpartition_inf[j] = calloc(sizeof(mp_limb_t) * limballoc, nlimit);
        cosums_inf_ptr[j] = cosumpartition_inf[j];
        
        mpn_and_n(cosums_inf_ptr[j], cosums_inf_ptr[j - 1], sums_ptr[j - 1], limballoc);
    }
    
    // Construction initale des co-sommes supérieures
    cosumpartition_sup[pmax - 1] = calloc(sizeof(mp_limb_t) * limballoc, 1);
    cosums_sup_ptr[pmax - 1] = cosumpartition_sup[pmax - 1];
    mpn_com(cosums_sup_ptr[pmax - 1], cosums_sup_ptr[pmax - 1], limballoc);
    
    for (unsigned long j = pmax - 1; j > 0; j--) {
        cosumpartition_sup[j - 1] = calloc(sizeof(mp_limb_t) * limballoc, nlimit);
        cosums_sup_ptr[j - 1] = cosumpartition_sup[j - 1];
        
        mpn_and_n(cosums_sup_ptr[j - 1], cosums_sup_ptr[j], sums_ptr[j], limballoc);
    }
    
    
    unsigned long iter_num = action->iter_num;  // Nombre d'itérations
    
    // Itération jusqu'à énumérer toutes les partions sans-sommes à au plus pmax huches
    while (n >= n0) {
        // Placer n+1 dans une des huches en conservant le caractère sans-somme
        
        if (p == pmax) {
            // Gérer les contraintes dûes aux co-sommes
            unsigned long nmax = nalloc;        // Plus grand entier pouvant être atteint par cette partition
            
            char has_improved = 1;
            
            while (has_improved) {
                has_improved = 0;
                
                for (unsigned long j = 0; j < p; j++) {
                    mpn_and_n(work1, cosums_inf_ptr[j], cosums_sup_ptr[j], limballoc);
                    
                    // Supprimer les entiers de [1, n]
                    /*mpn_zero(work2, limballoc);
                    schur_number_setinterval_1(work2, limballoc, ((n+1) / GMP_NUMB_BITS) + 1, n+1);
                    mpn_andn_n(work1, work1, work2, limballoc);*/
                    schur_number_intersect_interval_n(work1, limballoc, limbsize, n + 1);
                    
                    // Supprimer les entiers de [nbest, +∞[
                    //schur_number_setinterval_s(work2, limballoc, (nbest / GMP_NUMB_BITS) + 1, nbest);
                    /*mpn_zero(work2, limballoc);
                    ADD_POINT(work2, nbest);
                    mpn_neg(work2, work2, limballoc);
                    mpn_andn_n(work1, work1, work2, limballoc);*/
                    schur_number_intersect_interval_0(work1, limballoc, INTEGER_2_LIMBSIZE(nbest - 1), nbest + 1);
                    
                    if (!mpn_zero_p(work1, limballoc)) {
                        while (!mpn_zero_p(work1, limballoc)) {
                            unsigned long x = mpn_scan1(work1, n + 1);
                            
                            // Regarder si x peut être ajouté
                            if (GET_POINT(sums_ptr[j], x)) {
                                if (nmax > x) {
                                    nmax = x;
                                }
                                break;
                            } else {
                                if (!GET_POINT(partition[j], x)) {
                                    has_improved = 1;
                                }
                                // Modifier la somme
                                mpn_zero(work2, limballoc);
                                schur_number_translation(work2, partition[j], limballoc, x);
                                mpn_ior_n(sums_ptr[j], sums_ptr[j], work2, limballoc);
                                
                                mpn_zero(work2, limballoc);
                                schur_number_ntranslation(work2, partitioninvert[j], limballoc, nalloc - x + 1);
                                mpn_ior_n(sums_ptr[j], sums_ptr[j], work2, limballoc);
                                
                                // Ajouter le point
                                ADD_POINT(partition[j], x);
                                ADD_POINT(partitioninvert[j], nalloc - x + 1);
                            }
                            DELETE_POINT(work1, x);
                        }
                        
                        
                        // Modifier les co-sommes inférieures
                        for (unsigned long k = j; k < p - 1; k++) {
                            mpn_and_n(cosums_inf_ptr[k + 1], cosums_inf_ptr[k], sums_ptr[k], limballoc);
                        }
                        
                        // Modifier les co-sommes supérieures
                        for (unsigned long k = j; k > 0; k--) {
                            mpn_and_n(cosums_sup_ptr[k - 1], cosums_sup_ptr[k], sums_ptr[k], limballoc);
                        }
                    }
                }
            }
            
            if (nmax <= nbest) {
                notSumFree = 1;
                i = p;
                nblocking = n;
            }
        }
        
        while (notSumFree && i < p) {
            // Regarder si n+1 appartient à la somme de l'ensemble i
            
            iter_num++;
            notSumFree = !!GET_POINT(sums_ptr[i], n + 1);
            
            if (notSumFree) {
                // Déterminer nblocking en cas de blocage
                // Il s'agit du plus petit m tel que la somme de l'ensemble i à l'étape m rencontre n
                
                unsigned long m = setmin[i];
                
                if (m < n0) {
                    m = n0;
                }
                
                m += find_nblocking(n + 1, sumpartition[i], limballoc, n - m + 1);
                
                if (m > nblocking) {
                    nblocking = m;
                }
            }
            
            i += notSumFree;
        }
        
        if (notSumFree || n >= nlimit) {
            // L'entier n+1 n'a pu être placé dans aucune huche
            
            if (p < pmax && i <= p && n < nlimit) {
                // Remplir une nouvelle huche
                ADD_POINT(partitioninvert[i], nalloc - n);
                n++;
                ADD_POINT(partition[p], n);
                cardinals[p] = 1;
                setmin[p] = n;
                
                // Mettre à jour la somme
                sums_ptr[p] = sumpartition[p];
                
                // Mettre à jour les autres sommes
                for (unsigned long j = 0; j < p; j++) {
                    mp_limb_t *sumset = sums_ptr[j] + limballoc;
                    mpn_copyi(sumset, sums_ptr[j], limballoc);
                    sums_ptr[j] = sumset;
                }
                
                // Recopier les co-sommes inférieures pour 0 < j < p
                for (unsigned long j = 1; j < p; j++) {
                    mp_limb_t *cosumset = cosums_inf_ptr[j] + limballoc;
                    mpn_copyi(cosumset, cosums_inf_ptr[j], limballoc);
                    cosums_inf_ptr[j] = cosumset;
                }
                
                if (p == pmax - 1) {
                    // Inutile de mettre à jour les co-sommes supérieures p > j >= 0 puisque toutes ces co-sommes doivent être nulles
                    for (unsigned long j = p; j > 0; j--) {
                        mpn_and_n(cosums_sup_ptr[j - 1], cosums_sup_ptr[j], sums_ptr[j], limballoc);
                    }
                }
                
                //nbest = action->func(NULL, (n << (pmax - p)) + sum_min, action);
                /*if (nbest < (n << (pmax - p)) + sum_min) {
                    nbest = (n << (pmax - p)) + sum_min;
                }*/
                sum_min += n;
                
                i = 0;
                if (n > nsize) {
                    limbsize++;
                    nsize += GMP_NUMB_BITS;
                }
                p++;
                is_new_branch = 1;
                nblocking = n0;
                
            } else {
                
                if (is_new_branch) {
                    // Agir sur la partition
                    action->work = work1;
                    action->sumset = sums_ptr[pmax - 1];
                    action->setmin = setmin[pmax - 1];
                    action->limballoc = limballoc;
                    nbest = action->func(partition, n, action);
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
                
                // Intersecter partition et partitioninvert
                mp_size_t blockinglimbsize = ((nblocking - 1) / GMP_NUMB_BITS) + 1;    // Nombre de limbes nécessaires pour contenir nblocking-1
                unsigned long nblockinginvert = nalloc - nblocking + 2;
                mp_size_t invert_blockinglimbsize = (nblockinginvert / GMP_NUMB_BITS) + 1;
                for (unsigned long i = 0; i < p; i++) {
                    schur_number_intersect_interval_n(partitioninvert[i], limballoc, invert_blockinglimbsize, nblockinginvert);
                    schur_number_intersect_interval_0(partition[i], limballoc, blockinglimbsize, nblocking - 1);
                }
                
                if (p == pmax) {
                    // Dépiler les co-sommes supérieures
                    mp_bitcnt_t limb_diff;
                    mp_size_t is_not_empty;
                    
                    if (nblocking > setmin[p - 1]) {
                        limb_diff = (n - nblocking) * limballoc;
                        is_not_empty = 1;
                    } else {
                        limb_diff = (n - setmin[p - 1]) * limballoc;
                        is_not_empty = 0;
                    }
                    
                    for (unsigned long j = 0; j < p - 1; j++) {
                        cosums_sup_ptr[j] -= limb_diff;
                        mpn_zero(cosums_sup_ptr[j], limb_diff + limballoc);
                        cosums_sup_ptr[j] -= is_not_empty * limballoc;
                    }
                }
                
                // Dépiler les sommes et co-sommes
                for (unsigned long i = 0; i < p; i++) {
                    // Il faut déterminer le nombre d'éléments de [nblocking, n] contenu dans l'ensemble i
                    mp_bitcnt_t cardinal = mpn_popcount(partition[i], blockinglimbsize);    // Nombre d'éléments dans A_i ∩ [1, nblocking-1]
                    
                    mp_bitcnt_t limb_diff;
                    if (nblocking > setmin[i]) {
                        limb_diff = (n - nblocking) * limballoc;
                    } else {
                        limb_diff = (n - setmin[i]) * limballoc;
                        sum_min -= setmin[i];
                        setmin[i] = 0;
                    }
                    
                    // Dépiler la somme
                    sums_ptr[i] -= limb_diff;
                    mpn_zero(sums_ptr[i], limb_diff + limballoc);
                    sums_ptr[i] -= (!!cardinal) * limballoc;
                    
                    if (i) {
                        // Dépiler la co-somme inférieure
                        cosums_inf_ptr[i] -= limb_diff;
                        mpn_zero(cosums_inf_ptr[i], limb_diff + limballoc);
                        cosums_inf_ptr[i] -= (!!cardinal) * limballoc;
                    }
                    
                    // Mettre à jour les informations cardinales
                    cardinals[i] = cardinal;
                }
                
                while (0 < p && !setmin[p - 1]) {
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
            
            // Mettre à jour les sommes
            for (unsigned long j = 0; j < p; j++) {
                mp_limb_t *sumset = sums_ptr[j] + limballoc;
                mpn_copyi(sumset, sums_ptr[j], limballoc);
                sums_ptr[j] = sumset;
            }
            schur_number_translation(work1, partition[i], limballoc, n);
            mpn_ior_n(sums_ptr[i], sums_ptr[i], work1, limballoc);
            
            ADD_POINT(partition[i], n);
            
            // Recopier les co-sommes inférieures pour 0 < j <= i
            for (unsigned long j = 1; j <= i; j++) {
                mp_limb_t *cosumset = cosums_inf_ptr[j] + limballoc;
                mpn_copyi(cosumset, cosums_inf_ptr[j], limballoc);
                cosums_inf_ptr[j] = cosumset;
            }
            
            // Mettre à jour les co-sommes inférieures pour i < j < p
            for (unsigned long j = i + 1; j < p; j++) {
                mp_limb_t *cosumset = cosums_inf_ptr[j] + limballoc;
                mpn_and_n(cosumset, cosums_inf_ptr[j - 1], sums_ptr[j - 1], limballoc);
                cosums_inf_ptr[j] = cosumset;
            }
            
            if (p == pmax) {
                // Si p < pmax, alors toutes les co-sommes supérieures sont nulles
                
                // Recopier les co-sommes supérieures pour p > j >= i
                for (unsigned long j = i; j < p - 1; j++) {
                    mp_limb_t *cosumset = cosums_sup_ptr[j] + limballoc;
                    mpn_copyi(cosumset, cosums_sup_ptr[j], limballoc);
                    cosums_sup_ptr[j] = cosumset;
                }
                
                // Mettre à jour les co-sommes supérieures pour i > j >= 0
                for (unsigned long j = i; j > 0; j--) {
                    mp_limb_t *cosumset = cosums_sup_ptr[j - 1] + limballoc;
                    mpn_and_n(cosumset, cosums_sup_ptr[j], sums_ptr[j], limballoc);
                    cosums_sup_ptr[j - 1] = cosumset;
                }
            }
            
            // Mettre à jour les cardinaux
            cardinals[i]++;
            
            i = 0;
            if (n >= nsize) {
                limbsize++;
                nsize += GMP_NUMB_BITS;
            }
            is_new_branch = 1;
            nblocking = n0;
        }
        notSumFree = 1;
    }
    
    // Nettoyage
    free(work1);
    
    free(cardinals);
    free(setmin);
    
    for (unsigned long j = 0; j < pmax; j++) {
        free(sumpartition[j]);
        free(cosumpartition_inf[j]);
        free(cosumpartition_sup[j]);
    }
    free(sumpartition);
    free(sums_ptr);
    free(cosumpartition_inf);
    free(cosums_inf_ptr);
    
    action->iter_num = iter_num;
    
    return nbest;
}
