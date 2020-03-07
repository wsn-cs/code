//
//  schurNumberIO.h
//  schurNumberRecursive
//
//  Created by rubis on 28/01/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#ifndef schurNumberIO_h
#define schurNumberIO_h

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

struct schurNumberIOAction {
    unsigned long p;        // Nombre d'ensembles par partition
    
    unsigned long nmax;     // Plus grande taille des partitions
    unsigned long iter_num; // Nombre d'itérations
    size_t count_max;       // Nombre de partitions de taille maximales
    size_t count_all;       // Nombre de partitions non prolongeables
    
    size_t count;           // Nombre de partitions contenues dans le tampon
    
    char *limbsize_buffer;      // Tampon contenant les nombres de limbes par ensembles de la partition
    size_t limbsize_size;       // Taille du tampon en octet
    FILE *limbsize_stream;      // Flux associé au tampon
    
    char *partition_buffer;     // Tampon contenant les partitions
    size_t partition_size;      // Taille du tampon en octet
    FILE *partition_stream;     // Flux associé au tampon
    
    void (*func)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
    
    size_t n_buffers;               // Nombre de tampons issus d'un éventuel regroupement
    size_t *count_a;                // Tableau contenant le nombre de partitions contenues dans chaque tampon
    char **limbsize_buffer_a;
    size_t *limbsize_size_a;
    char **partition_buffer_a;      // Tableau de n_buffers tampons contenant les partitions
    size_t *partition_size_a;       // Tableau des tailles des n_buffers tampon en octet
};

typedef struct schurNumberIOAction schur_number_action_t;

typedef void (*schur_number_action_func_t)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);

unsigned long schurNumberGetSetMaximum(char *str);
void schurNumberGetSet(char *str, mp_limb_t *set, mp_limb_t *setinvert, mp_size_t limballoc);
unsigned long schurNumberGetPartition(unsigned long p, char **str, mp_limb_t **partition, mp_limb_t **partitioninvert, mp_size_t limballoc);

void schurNumberPrintSet(unsigned long n, mp_limb_t *set);
void schurNumberPrintPartition(unsigned long p, unsigned long n, mp_limb_t **partition);

void schurNumberActionAlloc(schur_number_action_t *action, unsigned long p, void (*func)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action));
void schurNumberActionDealloc(schur_number_action_t *action);

void schurNumberActionGatherCopy(schur_number_action_t *r_action, schur_number_action_t *actions, size_t n_actions);
void schurNumberActionGatherNoCopy(schur_number_action_t *r_action, schur_number_action_t *actions, size_t n_actions);

void schurNumberDefaultAction(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
void schurNumberSaveOnePartition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
void schurNumberSaveBestPartition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
void schurNumberSaveAllPartition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
size_t schurNumberPrintPartitions(struct schurNumberIOAction *action);

#endif /* schurNumberIO_h */
