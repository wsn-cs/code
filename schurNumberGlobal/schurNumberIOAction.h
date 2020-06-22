//
//  schurNumberIOAction.h
//  SchurNumber
//
//  Created by rubis on 25/03/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#ifndef schurNumberIOAction_h
#define schurNumberIOAction_h

#include "schurNumberIO.h"
#include "schurNumberSaving.h"

enum schur_number_action_flag {SCHUR_NUMBER_DEFAULT, SCHUR_NUMBER_INTERMEDIATE_PRINT, SCHUR_NUMBER_PRINT_AND_STOP};

struct schurNumberIOAction {
    unsigned long p;        // Nombre d'ensembles par partition
    
    unsigned long nbest;     // Plus grande taille des partitions
    unsigned long iter_num; // Nombre d'itérations
    size_t count_max;       // Nombre de partitions de taille maximales
    size_t count_all;       // Nombre de partitions non prolongeables
    
    size_t count;           // Nombre de partitions contenues dans le tampon
    size_t count_limit;     // Nombre limite de partitions pouvant être contenues dans le tampon
    size_t size_limit;
    
    char *limbsize_buffer;      // Tampon contenant les nombres de limbes par ensembles de la partition
    size_t limbsize_size;       // Taille du tampon en octet
    FILE *limbsize_stream;      // Flux associé au tampon
    mp_size_t limbsize;
    
    char *partition_buffer;     // Tampon contenant les partitions
    size_t partition_size;      // Taille du tampon en octet
    FILE *partition_stream;     // Flux associé au tampon
    
    enum schur_number_action_flag action_flag;          // Drapeau indiquant si, au cours de l'exécution, les partitions doivent être imprimées
    
    char *sum_partition_buffer;                 // Tampon contenant les sommes du dernier ensemble des partitions
    size_t sum_partition_size;                  // Taille du tampon en octet
    FILE *sum_partition_stream;                 // Flux associé au tampon
    char *sorted_index_sum_partition_buffer;    // Tampon contenant les indices ordonnées des sommes du dernier ensemble des partitions
    size_t sorted_index_sum_partition_size;     // Taille du tampon en octet
    FILE *sorted_index_sum_partition_stream;    // Flux associé au tampon
    //unsigned long depth;                        // Profondeur de la partition
    
    mp_limb_t *work;
    
    unsigned long (*func)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
    
    schur_number_intermediate_save_t *save;
    
    size_t count_gathered_actions;                  // Nombre d'actions issues d'un éventuel regroupement
    struct schurNumberIOAction **gathered_actions;  // Tableau des actions issues d'un éventuel regroupement
};

typedef struct schurNumberIOAction schur_number_action_t;

typedef unsigned long (*schur_number_action_func_t)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);

void schur_number_action_alloc(schur_number_action_t *action, unsigned long p, unsigned long (*func)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action));
void schur_number_action_dealloc(schur_number_action_t *action);

size_t schur_number_action_print_partitions(schur_number_action_t *action);

unsigned long schur_number_action_total_iterations(const schur_number_action_t *action);
unsigned long schur_number_action_total_count_all(const schur_number_action_t *action);
unsigned long schur_number_action_total_count_max(const schur_number_action_t *action);
unsigned long schur_number_action_total_Nmax(const schur_number_action_t *action);

void schur_number_action_gather_copy(schur_number_action_t *r_action, schur_number_action_t **actions, size_t n_actions);
void schur_number_action_gather_nocopy(schur_number_action_t *r_action, schur_number_action_t **actions, size_t n_actions);

unsigned long schur_number_default_action(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
unsigned long schur_number_default_action_sync(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
unsigned long schur_number_sync_action(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
unsigned long schur_number_save_some_partition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
unsigned long schur_number_save_best_partition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
unsigned long schur_number_save_all_partition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);

unsigned long schur_number_save_distinct_sum_partition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
unsigned long schur_number_save_distinct_restrictedsum_partition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);

#endif /* schurNumberIOAction_h */
