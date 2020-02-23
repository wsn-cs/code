//
//  schurNumberSymmetricHeader.h
//  schurNumberSymmetric
//
//  Created by wsn-cs on 30/05/2019.
//

#ifndef schurNumberSymmetricHeader_h
#define schurNumberSymmetricHeader_h

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

void printSymmetricPartition(unsigned int p, unsigned long n, mp_limb_t **partition);

unsigned long schurNumberSymmetricBranch(unsigned long n, partition_t *partitionstruc, unsigned long depth,  char completesearch);

unsigned long schurNumberPermutationImproved(unsigned long n, partition_t *partitionstruc, unsigned long depth,  char completesearch);

unsigned long schurNumberSymmetricImposedPartition(unsigned long n, unsigned long p);

#endif
