//
//  schurNumberSaveDistinctBegin.c
//  SchurNumber
//
//  Created by rubis on 26/06/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberIOAction.h"
#include "schurNumberPartitionStruc.h"

//  Définition des fonctions de recherche pour la relation d'ordre donnée
//  par le plus petit élément du dernier ensemble.
#define FIND_INF_INDEX_FUNC find_inf_index_setmin_order
#define FIND_SUP_INDEX_FUNC find_sup_index_setmin_order
#define FIND_EXTREMAL_INDICES_FUNC find_extremal_indices_setmin_order

#define INDEXING_FUNC(i) (sorted_indices[i])
#define CMP_FUNC(cmp_res, min1, min2) do{\
            cmp_res = (min1) - (min2);\
        } while(0)
#define FUNC_PARAMS size_t *sorted_indices
#define FUNC_ARGS sorted_indices

#include "schurNumberFindSort.h"
#undef schurNumberFindSort_h

#undef FIND_INF_INDEX_FUNC
#undef FIND_SUP_INDEX_FUNC
#undef FIND_EXTREMAL_INDICES_FUNC

#undef INDEXING_FUNC
#undef CMP_FUNC
#undef FUNC_PARAMS
#undef FUNC_ARGS

//  Définition des fonctions de recherche pour l'ordre lexicographique.
#define FIND_INF_INDEX_FUNC find_inf_index_lex_order
#define FIND_SUP_INDEX_FUNC find_sup_index_lex_order
#define FIND_EXTREMAL_INDICES_FUNC find_extremal_indices_lex_order

#define INDEXING_FUNC(i) (p * limbsize * sorted_indices[i])
// Attention, ce qui suit dépend de la manière dont sont organisés les ensembles au sein des partitions
#define CMP_FUNC(cmp_res, part1, part2) do {\
        unsigned long j = p;\
        do {\
            j--;\
            mp_limb_t *set1 = part1 + j * limbsize;\
            mp_limb_t *set2 = part2 + 2 * j * limballoc;\
            cmp_res = ((set1[minlimbsize - 1] & mask) - (set2[minlimbsize - 1] & mask));\
            if (!cmp_res && minlimbsize > 1) {\
                cmp_res = mpn_cmp(set1, set2, minlimbsize - 1);\
            }\
        } while (!cmp_res && j > 0);\
    } while(0)
#define FUNC_PARAMS size_t *sorted_indices, unsigned long p, mp_size_t limbsize, mp_size_t limballoc, mp_size_t minlimbsize, mp_limb_t mask
#define FUNC_ARGS sorted_indices, p, limbsize, limballoc, minlimbsize, mask

#include "schurNumberFindSort.h"


struct sum_additional_data_struc {
    char *sum_partition_buffer;                 // Tampon contenant les sommes du dernier ensemble des partitions
    size_t sum_partition_size;                  // Taille du tampon en octet
    FILE *sum_partition_stream;                 // Flux associé au tampon
    
    char *sorted_index_sum_partition_buffer;    // Tampon contenant les indices ordonnées des sommes du dernier ensemble des partitions
    size_t sorted_index_sum_partition_size;     // Taille du tampon en octet
    FILE *sorted_index_sum_partition_stream;    // Flux associé au tampon
};

unsigned long schur_number_save_distinct_begin_partition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /* Ajoute la partition si le début de ses p-1 premiers ensemble n'a pas déjà été trouvée, et si n >= nthreshold. */
    unsigned long  p = action->p;
    
    schur_number_intermediate_save_t *save = action->save;
    
    if (partition) {
        action->count_all++;
        
        if (save && !(action->count_all % SCHUR_NUMBER_SAVE_FREQUENCY)) {
            unsigned long nbest = schur_number_save_progression_update(save, n, partition);
            if (nbest > action->nbest) {
                schur_number_save_distinct_sum_partition(NULL, nbest, action);
            }
        }
    }
    
    // Regarder si nbest doit être modifié
    if (n > action->nbest) {
        action->nbest = n;
        if (save) {
            action->nbest = schur_number_save_best_upgrade(save, n, partition);
        }
        action->count_max = 0;
    }
    
    if (n >= action->nthreshold && partition) {
        if (action->partition_size < action->size_limit) {
            // Récupérer les données passées en argument supplémentaire
            unsigned long setmin = action->setmin;
            
            FILE *limbsize_stream = action->limbsize_stream;
            FILE *partition_stream = action->partition_stream;
            
            FILE *setmin_stream = action->auxiliary_stream;     // Flux vers les profondeurs de chaque partition
            FILE *sorted_index_sum_partition_stream = action->sorted_index_sum_partition_stream;
            
            /* Regarder si la somme de la partition a déjà été trouvée. */
            //mp_size_t limbsize = ((unsigned long)n>>6) + 1;
            mp_size_t limballoc = action->limballoc;
            mp_size_t limbsize = limballoc;
            //schur_number_sumset(work, partition[p - 1], partition[p - 1], 2 * limbsize, limbsize, 0, work + 2 * limbsize);
            
            fflush(partition_stream);
            fflush(setmin_stream);
            fflush(sorted_index_sum_partition_stream);
            
            size_t *sorted_indices = action->sorted_index_sum_partition_buffer;
            
            size_t i1 = 0;
            size_t i2 = action->count - 1;
            size_t i_card = 0;
            
            if (action->count > 0) {
                // Trouver le nombre i_card d'indices i tels que min(partition_buffer[limballoc * sorted_indices[i]]) == setmin
                
                /*cmp_block_t setmin_cmp = ^long(mp_limb_t *min1_p, mp_limb_t *min2_p) {
                    return *min1_p - setmin;
                };
                
                indexing_block_t sorted_index = ^size_t(size_t i) {
                    return sorted_indices[i];
                };*/
                
                i_card = find_extremal_indices_setmin_order(NULL, action->auxiliary_buffer, sorted_indices, &i1, &i2);
            }
            
            if (i_card) {
                // Regarder si des partitions coïncident dans l'intervalle [1, setmin-1]
                mp_size_t minlimbsize = INTEGER_2_LIMBSIZE(setmin + (setmin >> 1) - (setmin >> 2));
                mp_limb_t mask = GMP_NUMB_MAX >> ((GMP_NUMB_BITS - setmin - (setmin >> 1) + (setmin >> 2) - 1) % GMP_NUMB_BITS);
                
                //mp_size_t minlimbsize = INTEGER_2_LIMBSIZE(setmin - 1);
                //mp_limb_t mask = GMP_NUMB_MAX >> ((GMP_NUMB_BITS - setmin) % GMP_NUMB_BITS);
                /*cmp_block_t lex_partition_cmp = ^long(mp_limb_t *part1, mp_limb_t *part2) {
                    long setwise_cmp;
                    unsigned long j = p ;//- 1;
                    do {
                        j--;
                        // Attention, ce qui suit dépend de la manière dont sont organisés les ensembles au sein des partitions
                        mp_limb_t *set1 = part1 + j * limbsize;
                        mp_limb_t *set2 = part2 + 2 * j * limballoc;
                        setwise_cmp = ((set1[minlimbsize - 1] & mask) - (set2[minlimbsize - 1] & mask));
                        if (!setwise_cmp && minlimbsize > 1) {
                            setwise_cmp = mpn_cmp(set1, set2, minlimbsize - 1);
                        }
                    } while (!setwise_cmp && j > 0);
                    return setwise_cmp;
                };*/
                
                i_card = find_extremal_indices_lex_order(*partition, action->partition_buffer, sorted_indices, p, limbsize, limballoc, minlimbsize, mask, &i1, &i2);
            }
            
            if (i_card) {
                // Réaliser l'intersection de la partition déjà enregistrée avec la nouvelle
                mp_limb_t *partition_buffer = action->partition_buffer;
                size_t index = sorted_indices[i1];
                mp_limb_t *partition2 = partition_buffer + index * limbsize * p;
                for (unsigned long j = 0; j < p; j++) {
                    mpn_and_n(&partition2[j * limbsize], &partition2[j * limbsize], partition[j], limbsize);
                }
            } else {
                // Somme pas encore présente, donc ajouter la partition.
                fwrite(&limbsize, sizeof(mp_size_t), 1, limbsize_stream);
                
                for (unsigned long j = 0; j < p; j++) {
                    fwrite(partition[j], sizeof(mp_limb_t), limbsize, partition_stream);
                }
                
                // Ajouter aussi le plus petit élément de la partie p
                fwrite(&setmin, sizeof(unsigned long), 1, setmin_stream);
                
                // Placer son indice en MAX(i1 + 1, i2)
                size_t count = action->count;
                size_t last_index = count;
                if (count > i2) {
                    last_index = sorted_indices[count - 1];
                    mpn_copyd(sorted_indices + i2 + 1, sorted_indices + i2, count - i2 - 1);
                    sorted_indices[i2] = count;
                }
                /*for (size_t i = i1 + 1; i < action->count; i++) {
                 size_t old_index = sorted_indices[i];
                 sorted_indices[i] = index;
                 index = old_index;
                 }*/
                fwrite(&last_index, sizeof(size_t), 1, sorted_index_sum_partition_stream);
                
                action->count ++;
            }
        }
        if (n == action->nbest) {
            action->count_max ++;
        }
    }
    
    return action->nbest;
}

