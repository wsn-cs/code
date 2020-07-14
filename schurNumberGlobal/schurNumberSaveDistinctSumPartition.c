//
//  schurNumberSaveDistinctSumPartition.c
//  SchurNumber
//
//  Created by rubis on 25/06/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberIOAction.h"
#include "schurNumberPartitionStruc.h"

//  Définition des fonctions de recherche pour la relation d'ordre donnée
//  par les sommes.
#define FIND_INF_INDEX_FUNC find_inf_index_sumset_order
#define FIND_SUP_INDEX_FUNC find_sup_index_sumset_order
#define FIND_EXTREMAL_INDICES_FUNC find_extremal_indices_sumset_order

#define INDEXING_FUNC(i) (limballoc * sorted_indices[i])
#define CMP_FUNC(cmp_res, set1, set2) do{\
            cmp_res = mpn_cmp((set1), (set2), (limballoc));\
        } while(0)
#define FUNC_PARAMS size_t *sorted_indices, mp_size_t limballoc
#define FUNC_ARGS sorted_indices, limballoc
#include "schurNumberFindSort.h"

struct sum_additional_data_struc {
    char *sum_partition_buffer;                 // Tampon contenant les sommes du dernier ensemble des partitions
    size_t sum_partition_size;                  // Taille du tampon en octet
    FILE *sum_partition_stream;                 // Flux associé au tampon
    
    char *sorted_index_sum_partition_buffer;    // Tampon contenant les indices ordonnées des sommes du dernier ensemble des partitions
    size_t sorted_index_sum_partition_size;     // Taille du tampon en octet
    FILE *sorted_index_sum_partition_stream;    // Flux associé au tampon
};

void schur_number_sum_action_alloc(schur_number_action_t *action) {
    action->auxiliary_stream = open_memstream(&(action->auxiliary_buffer), &(action->auxiliary_size));
    action->sorted_index_sum_partition_stream = open_memstream(&(action->sorted_index_sum_partition_buffer), &(action->sorted_index_sum_partition_size));
    action->limbsize = PARTITION_2_LIMBSIZE(action->p);
}

unsigned long schur_number_save_distinct_sum_partition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /* Ajoute la partition si sa somme n'a pas déjà été trouvée, et si n == nbest. */
    unsigned long  p = action->p;
    FILE *limbsize_stream = action->limbsize_stream;
    FILE *partition_stream = action->partition_stream;
    
    FILE *sum_partition_stream = action->auxiliary_stream;
    FILE *sorted_index_sum_partition_stream = action->sorted_index_sum_partition_stream;
    
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
    /* else if (save) {
        // Assure que nbest = nbest_global dans tous les cas
        unsigned long nbest_global  = schur_number_save_get_best_global(save);
        if (nbest_global > action->nbest) {
            action->nbest = nbest_global;
            action->count_max = 0;
        }
    }*/
    
    if (n >= action->nthreshold && partition) {
        if (action->partition_size < action->size_limit) {
            mp_limb_t *work = action->work;
            mp_limb_t *sumset = action->sumset;
            
            /* Regarder si la somme de la partition a déjà été trouvée. */
            //mp_size_t limbsize = ((unsigned long)n>>6) + 1;
            mp_size_t limballoc = action->limballoc;
            mp_size_t limbsize = limballoc;
            //schur_number_sumset(work, partition[p - 1], partition[p - 1], 2 * limbsize, limbsize, 0, work + 2 * limbsize);
            
            fflush(sum_partition_stream);
            fflush(sorted_index_sum_partition_stream);
            
            size_t *sorted_indices = action->sorted_index_sum_partition_buffer;
            
            size_t i1 = 0;
            size_t i2 = action->count - 1;
            size_t i_card = 0;
            
            if (action->count > 0) {
                // Trouver le nombre i_card d'indices i tels que sum_partition_buffer[limballoc * sorted_indices[i]] == sumset
                
                /*cmp_block_t set_cmp = ^long(mp_limb_t *set1, mp_limb_t *set2) {
                    return mpn_cmp(set1, set2, limballoc);
                };
                
                indexing_block_t sorted_index = ^size_t(size_t i) {
                    return limballoc * sorted_indices[i];
                };*/
                
                i_card = find_extremal_indices_sumset_order(sumset, action->auxiliary_buffer, sorted_indices, limballoc, &i1, &i2);
            }
            /*
            if (i_card) {
                // Regarder si des partitions coïncident dans l'intervalle [1, setmin-1]
                fflush(action->partition_stream);
                unsigned long setmin = action->setmin;
                mp_size_t minlimbsize = INTEGER_2_LIMBSIZE(setmin - 1);
                mp_limb_t num_max = GMP_NUMB_MAX;
                mp_limb_t mask = GMP_NUMB_MAX >> ((GMP_NUMB_BITS - setmin) % GMP_NUMB_BITS);
                cmp_block_t lex_partition_cmp = ^long(mp_limb_t *part1, mp_limb_t *part2) {
                    long setwise_cmp;
                    unsigned long j = p - 1;
                    do {
                        j--;
                        // Attention, ce qui suit dépend de la manière dont sont organisés les ensembles au sein des partitions
                        mp_limb_t *set1 = part1 + j * limbsize;
                        mp_limb_t *set2 = part2 + 2 * j * limballoc;
                        setwise_cmp = ((set1[minlimbsize - 1] & mask) == (set2[minlimbsize - 1] & mask));
                        if (minlimbsize > 1) {
                            setwise_cmp = setwise_cmp && mpn_cmp(set1, set2, minlimbsize - 1);
                        }
                    } while (!setwise_cmp && j > 0);
                    return setwise_cmp;
                };
                i_card = find_extremal_indices2(*partition, action->partition_buffer, lex_partition_cmp, ^size_t(size_t i) {return p * limbsize * sorted_indices[i];}, &i1, &i2);
            }*/
            
            if (i_card) {
                // Réaliser l'intersection de la partition déjà enregistrée avec la nouvelle
                fflush(partition_stream);
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
                
                // Ajouter sa somme
                fwrite(work, sizeof(mp_limb_t), limballoc, sum_partition_stream);
                
                // Placer son indice en MIN(i1 + 1, i2)
                //i1++;
                size_t count = action->count;
                size_t last_index = count;
                if (count > i2) {
                    last_index = sorted_indices[count - 1];
                    mpn_copyi(sorted_indices + i2 + 1, sorted_indices + i2, count - i2);
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
