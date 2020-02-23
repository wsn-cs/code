//
//  schurNumberNestedMonteCarloHeader.h
//  SchurNumber
//
//  Created by Gabriel Merlin on 22/04/2019.
//
//  Ce fichier est un header déclarant toutes les fonctions de calcul.
//

#ifndef schurNumberNestedMonteCarloHeader
#define schurNumberNestedMonteCarloHeader

#include <stdlib.h>
#include <gmp.h>

struct schur_partition_struct {
    unsigned int pmax;  // Nombre maximal de huches
    unsigned int p;     // Nombre de huches non vides
    
    unsigned long n;    // Entier courant
    mp_size_t limballoc;// Nombre de limbes alloués à chaque huche
    mp_size_t limbsize; // Nombre de limbes utilisés par chaque huche
    
    mp_limb_t **partition;
    mp_limb_t **partitioninvert;
    
    mp_limb_t *work0;
    mp_limb_t *work1;
};

typedef struct schur_partition_struct partition_t;

/*Fonctions d'initialisation, de réallocation et de libération*/
void partition_init(unsigned int pmax, unsigned long nmax, partition_t *partitionstruc);
void partition_realloc(partition_t *partitionstruc, mp_limb_t **partitionbest);
void partition_unalloc(partition_t *partitionstruc);
void partition_copy(partition_t *partitionstrucd, partition_t *partitionstrucs);

/*Fonctions d'affichage*/
void printPartition(unsigned int p, unsigned long n, mp_limb_t **partition);
unsigned int schurNumberScanPartitionFromFile(char *filename, partition_t *partitionstruc);
unsigned long schurNumberPrintPartitionToFile(char *filename, unsigned int p, unsigned long n, mp_limb_t **partition);

/*Fonction ajoutant les entiers un par un*/
unsigned long schurNumberSimpleMonteCarloLevelIteration(partition_t *sfpartitionstruc, unsigned int level, mp_limb_t **sfpartitionbest, unsigned int *pbestptr, unsigned int simulnum0);

unsigned long schurNumberWeakSimpleMonteCarloLevelIteration(partition_t *sfpartitionstruc, unsigned int level, mp_limb_t **sfpartitionbest, unsigned int *pbestptr, unsigned int simulnum0);

#endif
