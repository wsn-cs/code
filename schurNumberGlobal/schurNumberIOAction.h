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

struct schurNumberIOAction {
    unsigned long p;        // Nombre d'ensembles par partition
    
    unsigned long nmax;     // Plus grande taille des partitions
    unsigned long iter_num; // Nombre d'itérations
    size_t count_max;       // Nombre de partitions de taille maximales
    size_t count_all;       // Nombre de partitions non prolongeables
    
    size_t count;           // Nombre de partitions contenues dans le tampon
    size_t count_limit;     // Nombre limite de partitions pouvant être contenues dans le tampon
    
    char *limbsize_buffer;      // Tampon contenant les nombres de limbes par ensembles de la partition
    size_t limbsize_size;       // Taille du tampon en octet
    FILE *limbsize_stream;      // Flux associé au tampon
    
    char *partition_buffer;     // Tampon contenant les partitions
    size_t partition_size;      // Taille du tampon en octet
    FILE *partition_stream;     // Flux associé au tampon
    
    unsigned long (*func)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
    
    schur_number_intermediate_save_t *save;
    
    size_t n_buffers;               // Nombre de tampons issus d'un éventuel regroupement
    size_t *count_a;                // Tableau contenant le nombre de partitions contenues dans chaque tampon
    char **limbsize_buffer_a;
    size_t *limbsize_size_a;
    char **partition_buffer_a;      // Tableau de n_buffers tampons contenant les partitions
    size_t *partition_size_a;       // Tableau des tailles des n_buffers tampon en octet
};

typedef struct schurNumberIOAction schur_number_action_t;

typedef unsigned long (*schur_number_action_func_t)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);

void schurNumberActionAlloc(schur_number_action_t *action, unsigned long p, void (*func)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action));
void schurNumberActionDealloc(schur_number_action_t *action);

void schurNumberActionGatherCopy(schur_number_action_t *r_action, schur_number_action_t *actions, size_t n_actions);
void schurNumberActionGatherNoCopy(schur_number_action_t *r_action, schur_number_action_t *actions, size_t n_actions);

unsigned long schurNumberDefaultAction(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
unsigned long schurNumberSaveSomePartition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
unsigned long schurNumberSaveBestPartition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
unsigned long schurNumberSaveAllPartition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);

size_t schurNumberPrintPartitions(struct schurNumberIOAction *action);

#endif /* schurNumberIOAction_h */
