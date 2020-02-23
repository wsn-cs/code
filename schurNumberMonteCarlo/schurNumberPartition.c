//
//  schurNumberPartition.c
//  schurNumberMonteCarlo
//
//  Created by wsn-cs on 09/05/2019.
//
// Ce fichier implémente toutes les fonctions gérant la mémoire associée à une partition.
//

#include "schurNumberNestedMonteCarloHeader.h"

void partition_init(unsigned int pmax, unsigned long nmax, partition_t *partitionstruc) {
    /*A partir de partitionstruc dont la mémoire correspondant à la structure partition_struc_t
     a déjà été allouée, crée les ensembles de la partition  et les remplit de 0.*/
    unsigned int i;
    mp_size_t limballoc;
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t **partition;
    mp_limb_t **partitioninvert;
    
    /*Allocation des grands entiers correspondant aux ensembles*/
    limballoc = (nmax / mp_bits_per_limb) + 1;
    work0 = calloc(limballoc, sizeof(mp_limb_t));
    work1 = calloc(limballoc, sizeof(mp_limb_t));
    partition = calloc(pmax, sizeof(mp_limb_t *));           //Tableau contenant la partition
    partitioninvert = calloc(pmax, sizeof(mp_limb_t *));     //Tableau contenant les ensembles "inverses" de la partition
    
    partitionstruc->work0 = work0;
    partitionstruc->work1 = work1;
    partitionstruc->partition = partition;
    partitionstruc->partitioninvert = partitioninvert;
    
    for (i=0; i<pmax; i++) {
        *partition = calloc(limballoc, sizeof(mp_limb_t));
        *partitioninvert = calloc(limballoc, sizeof(mp_limb_t));
        partition++;
        partitioninvert++;
    }
    
    /*Initialisation des variables*/
    partitionstruc->pmax = pmax;
    partitionstruc->p = 0;
    partitionstruc->n = 0;
    partitionstruc->limballoc = limballoc;
    partitionstruc->limbsize = 0;
}

#define REALLOC(ptr, reptr, size, resize) \
do { \
reptr = realloc(ptr, resize);\
while (reptr == NULL && resize > size) { \
resize -= resize >> 1;\
reptr = realloc(ptr, resize);\
} \
if (reptr != NULL) { \
ptr = reptr;\
} \
} while(0)

void partition_realloc(partition_t *partitionstruc, mp_limb_t **partitionbest) {
    /*Augmente la taille de la partition ainsi que l'ensemble partitionbest.
     Tente d'abord d'allouer une partition de taille double, puis diminue la taille
     au fur et à mesure des échecs. La fonction n'est pas terrible.*/
    unsigned int i, p;
    mp_size_t limballoc, limbrealloc;
    mp_limb_t *allocptr, *reallocptr;
    mp_limb_t **partition, **partitioninvert;
    
    limballoc = partitionstruc->limballoc;
    limbrealloc = 2*limballoc;
    p = partitionstruc->pmax;
    
    allocptr = partitionstruc->work0;
    REALLOC(allocptr, reallocptr, limballoc, limbrealloc);
    partitionstruc->work0 = allocptr;
    allocptr = partitionstruc->work1;
    REALLOC(allocptr, reallocptr, limballoc, limbrealloc);
    partitionstruc->work1 = allocptr;
    
    partition = partitionstruc->partition;
    partitioninvert = partitionstruc->partitioninvert;
    for (i=0; i<p; i++) {
        allocptr = *partition;
        REALLOC(allocptr, reallocptr, limballoc, limbrealloc);
        partition = &allocptr;
        
        allocptr = *partitioninvert;
        REALLOC(allocptr, reallocptr, limballoc, limbrealloc);
        partitioninvert = &allocptr;
        
        allocptr = *partitionbest;
        REALLOC(allocptr, reallocptr, limballoc, limbrealloc);
        partitionbest = &allocptr;
        
        partition++;
        partitioninvert++;
        partitionbest++;
    }
    for (i=0; i<p; i++) {
        partitioninvert--;
        mpn_lshift(*partitioninvert, *partitioninvert, limbrealloc, limbrealloc - limballoc);
    }
    partitionstruc->limballoc = limbrealloc;
}

void partition_unalloc(partition_t *partitionstruc) {
    /*Libère tous les grands entiers associés aux ensembles de partitionstruc.*/
    unsigned i, pmax;
    mp_limb_t **partition;
    mp_limb_t **partitioninvert;
    
    pmax = partitionstruc->pmax;
    partition = partitionstruc->partition;
    partitioninvert = partitionstruc->partitioninvert;
    for (i=0; i<pmax; i++) {
        free(*partition);
        free(*partitioninvert);
        partition++;
        partitioninvert++;
    }
    free(partitionstruc->partition);
    free(partitionstruc->partitioninvert);
    free(partitionstruc->work0);
    free(partitionstruc->work1);
}

void partition_copy(partition_t *partitionstrucd, partition_t *partitionstrucs) {
    /*Copier partitionstrucs dans partitionstrucs.*/
    unsigned int pmax;
    unsigned int p;
    unsigned int i;
    mp_size_t limballoc;
    mp_size_t limbsize;
    mp_size_t limballocs;
    mp_limb_t **partitions, **partitiond;
    mp_limb_t **partitioninverts, **partitioninvertd;
    
    pmax = partitionstrucd->pmax;
    limballoc = partitionstrucd->limballoc;
    
    p = partitionstrucs->p;
    limballocs = partitionstrucs->limballoc;
    limbsize = partitionstrucs->limbsize;
    
    if (p > pmax || limbsize > limballoc) {
        /*Augmenter la taile de partitionstrucd.
         A faire*/
    }
    
    partitions = partitionstrucs->partition;
    partitioninverts = partitionstrucs->partitioninvert;
    partitiond = partitionstrucd->partition;
    partitioninvertd = partitionstrucd->partitioninvert;
    for (i=0; i<p; i++) {
        mpn_zero(*partitiond, limballoc);
        mpn_zero(*partitioninvertd, limballoc);
        mpn_copyd(*partitiond, *partitions, limbsize);
        mpn_copyd(*partitioninvertd + limballoc - limbsize, *partitioninverts + limballocs - limbsize, limbsize);
        partitions++;
        partitiond++;
        partitioninverts++;
        partitioninvertd++;
    }
    
    partitionstrucd->p = p;
    partitionstrucd->n = partitionstrucs->n;
    partitionstrucd->limbsize = limbsize;
}
