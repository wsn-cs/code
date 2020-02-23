//
//  tests.c
//  schurNumberMonteCarlo
//
//  Created by rubis on 19/04/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include <stdlib.h>
#include <gmp.h>
#include <stdio.h>
#include "schurNumberNestedMonteCarloHeader.h"

#if GMP_NUMB_BITS == 64
#define GMP_2EXP_NUMB_BITS 6
#elif GMP_NUMB_BITS == 128
#define GMP_2EXP_NUMB_BITS 7
#elif GMP_NUMB_BITS == 32
#define GMP_2EXP_NUMB_BITS 5
#endif

void printPartition(unsigned long p, unsigned long n, mp_limb_t **partition) {
    /*Affiche une partition.*/
    unsigned long limbn;
    mp_size_t limbsize;
    unsigned long i;
    mp_bitcnt_t j;
    mp_limb_t *set;
    mp_limb_t limb;
    limbsize = (n>>6) + 1;
    printf("Partition:\n");
    for (i=0; i<p; i++) {
        printf("\t");
        set = partition[i];
        for (limbn=0; limbn<limbsize; limbn++) {
            limb = set[limbn];
            for (j=0; j<mp_bits_per_limb; j++) {
                if (limb & ((unsigned long)1<<j)) {
                    printf(" %lu", limbn*mp_bits_per_limb + j + 1);
                }
            }
        }
        printf("\n");
    }
}

void partitionunstack(partition_t *sfpartitionstruc, mp_size_t limbsize, unsigned int p) {
    /*Revient à la partition initiale.
     Après ajout, la partition a pour taille limbsize et compte p huches.*/
    unsigned int i;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    mp_size_t limballoc;
    mp_size_t wlimbsize;
    unsigned int shift;
    mp_limb_t mask1, mask2;
    
    /*Mise en place de la partition*/
    limballoc = sfpartitionstruc->limballoc;
    sfpartition = sfpartitionstruc->partition;
    sfpartitioninvert = sfpartitionstruc->partitioninvert;
    
    wlimbsize = sfpartitionstruc->limbsize;
    shift = mp_bits_per_limb - (sfpartitionstruc->n)%mp_bits_per_limb;
    mask1 = (mp_limb_t)GMP_NUMB_MAX>>shift;
    mask2 = (mp_limb_t)GMP_NUMB_MAX<<shift;
    for (i=0; i<p; i++) {
        mpn_zero(sfpartition[i] + wlimbsize, limbsize - wlimbsize);
        mpn_zero(sfpartitioninvert[i] + limballoc - limbsize -1, limbsize - wlimbsize);
        //mpn_zero(sfpartitioninvert[i] + limballoc - limbsize, limbsize - wlimbsize);
        sfpartition[i][wlimbsize-1] &= mask1;
        sfpartitioninvert[i][limballoc - wlimbsize] &= mask2;
    }
}

void setunstack(mp_limb_t *set, mp_limb_t *setinvert, unsigned long oldn, mp_limb_t curlimbsize, mp_limb_t limballoc) {
    mp_size_t wlimbsize;
    unsigned int shift;
    mp_limb_t mask1, mask2;
    
    wlimbsize = (oldn - 1)/mp_bits_per_limb + 1;
    shift = mp_bits_per_limb - oldn % mp_bits_per_limb;
    mask1 = (mp_limb_t)GMP_NUMB_MAX>>shift;
    mask2 = (mp_limb_t)GMP_NUMB_MAX<<shift;
    mpn_zero(set + wlimbsize, curlimbsize - wlimbsize);
    mpn_zero(setinvert + limballoc - curlimbsize, curlimbsize - wlimbsize);
    set[wlimbsize-1] &= mask1;
    setinvert[limballoc - wlimbsize] &= mask2;
}

void setappend(mp_limb_t *set, mp_limb_t *setinvert, unsigned long n, mp_size_t limballoc) {
    mp_size_t limbsize;
    unsigned long nmodbpl;
    unsigned int shift;
    mp_limb_t mask1, mask2;
    
    limbsize = n/mp_bits_per_limb + 1;
    nmodbpl = n%mp_bits_per_limb;
    shift = mp_bits_per_limb - nmodbpl;
    
    mask1 = (mp_limb_t)1<<nmodbpl;
    mask2 = (mp_limb_t)1<<(shift - 1);
    set[limbsize -1] |= mask1;
    setinvert[limballoc - limbsize] |= mask2;
}

void testPartitionunstack() {
    unsigned int i;
    mp_size_t limballoc;
    partition_t sfpartitionstruc;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    
    /*Initialisation*/
    limballoc = 4;
    sfpartition = calloc(4, sizeof(mp_limb_t *));
    sfpartitioninvert = calloc(4, sizeof(mp_limb_t *));
    for (i=0; i<4; i++) {
        sfpartition[i] = calloc(limballoc, sizeof(mp_limb_t));
        sfpartitioninvert[i] = calloc(limballoc, sizeof(mp_limb_t));
    }
    sfpartitionstruc.limballoc = limballoc;
    sfpartitionstruc.limbsize = 1;
    sfpartitionstruc.partition = sfpartition;
    sfpartitionstruc.partitioninvert = sfpartitioninvert;
    sfpartitionstruc.pmax = 4;
    
    /*Création de la partition initiale,
     égale à {}, {}, ø, ø*/
}

char testSetunstack() {
    char b;
    mp_limb_t *set;
    mp_limb_t *setinvert;
    
    /*Initialisation*/
    set = calloc(3, sizeof(mp_limb_t));
    setinvert = calloc(3, sizeof(mp_limb_t));
    
    /*Valeur initiale {1,…,65}.
     Valeur finale {1,…,68}*/
    *set = GMP_NUMB_MAX;
    set[1] = (mp_limb_t)15;
    setinvert[2] = GMP_NUMB_MAX;
    setinvert[1] = GMP_NUMB_MAX ^ (GMP_NUMB_MAX >> 4);
    printf("Ensemble final 1:\n");
    printf("\t%lx %lx %lx\n", set[2], set[1], *set);
    printf("\t%lx %lx %lx\n", setinvert[2], setinvert[1], *setinvert);
    setunstack(set, setinvert, 65, 2, 3);
    printf("Ensemble dépilé 1:\n");
    printf("\t%lx %lx %lx\n", set[2], set[1], *set);
    printf("\t%lx %lx %lx\n", setinvert[2], setinvert[1], *setinvert);
    b = (*set == GMP_NUMB_MAX && set[1] == (mp_limb_t)1 && setinvert[2] == GMP_NUMB_MAX && setinvert[1] == (GMP_NUMB_MAX ^ (GMP_NUMB_MAX >> 1)));
    printf("%i\n", b);
    
    /*Valeur initiale {1,…,64}.
     Valeur finale {1,…,68}*/
    *set = GMP_NUMB_MAX;
    set[1] = (mp_limb_t)15;
    setinvert[2] = GMP_NUMB_MAX;
    setinvert[1] = GMP_NUMB_MAX ^ (GMP_NUMB_MAX >> 4);
    printf("Ensemble final 2:\n");
    printf("\t%lx %lx %lx\n", set[2], set[1], *set);
    printf("\t%lx %lx %lx\n", setinvert[2], setinvert[1], *setinvert);
    setunstack(set, setinvert, 64, 2, 3);
    printf("Ensemble dépilé 2:\n");
    printf("\t%lx %lx %lx\n", set[2], set[1], *set);
    printf("\t%lx %lx %lx\n", setinvert[2], setinvert[1], *setinvert);
    b = (*set == GMP_NUMB_MAX && set[1] == (mp_limb_t)0 && setinvert[2] == GMP_NUMB_MAX && setinvert[1] == (mp_limb_t)0);
    printf("%i\n", b);
    
    /*Valeur initiale {1,…,62}.
     Valeur finale {1,…,68}*/
    *set = GMP_NUMB_MAX;
    set[1] = (mp_limb_t)15;
    setinvert[2] = GMP_NUMB_MAX;
    setinvert[1] = GMP_NUMB_MAX ^ (GMP_NUMB_MAX >> 4);
    printf("Ensemble final 3:\n");
    printf("\t%lx %lx %lx\n", set[2], set[1], *set);
    printf("\t%lx %lx %lx\n", setinvert[2], setinvert[1], *setinvert);
    setunstack(set, setinvert, 62, 2, 3);
    printf("Ensemble dépilé 3:\n");
    printf("\t%lx %lx %lx\n", set[2], set[1], *set);
    printf("\t%lx %lx %lx\n", setinvert[2], setinvert[1], *setinvert);
    b = (*set == (GMP_NUMB_MAX>>2) && set[1] == (mp_limb_t)0 && setinvert[2] == (GMP_NUMB_MAX<<2) && setinvert[1] == (mp_limb_t)0);
    printf("%i\n", b);
    
    /*Valeur initiale {1,…,60}.
     Valeur finale {1,…,62}*/
    *set = GMP_NUMB_MAX;
    set[1] = (mp_limb_t)0;
    setinvert[2] = GMP_NUMB_MAX;
    setinvert[1] = (mp_limb_t)0;
    printf("Ensemble final 4:\n");
    printf("\t%lx %lx %lx\n", set[2], set[1], *set);
    printf("\t%lx %lx %lx\n", setinvert[2], setinvert[1], *setinvert);
    setunstack(set, setinvert, 60, 1, 3);
    printf("Ensemble dépilé 4:\n");
    printf("\t%lx %lx %lx\n", set[2], set[1], *set);
    printf("\t%lx %lx %lx\n", setinvert[2], setinvert[1], *setinvert);
    b = (*set == (GMP_NUMB_MAX>>4) && setinvert[2] == (GMP_NUMB_MAX<<4));
    printf("%i\n", b);
    
    return b;
}

char testSetappend() {
    char b;
    mp_limb_t *set;
    mp_limb_t *setinvert;
    
    /*Initialisation*/
    set = calloc(2, sizeof(mp_limb_t));
    setinvert = calloc(2, sizeof(mp_limb_t));
    
    /*Ensemble {1}*/
    *set = (mp_limb_t)1;
    setinvert[1] = GMP_NUMB_MAX ^ (GMP_NUMB_MAX >> 1);
    setappend(set, setinvert, 1, 2);
    printf("Ensemble 1:\n");
    printf("\t%lx %lx\n", set[1], *set);
    printf("\t%lx %lx\n", setinvert[1], *setinvert);
    b = (*set == 3 && setinvert[1] == (GMP_NUMB_MAX ^ (GMP_NUMB_MAX >> 2)));
    printf("%i\n", b);
    
    /*Ensemble {1,…,63}*/
    *set = GMP_NUMB_MAX>>1;
    setinvert[1] = GMP_NUMB_MAX<<1;
    setappend(set, setinvert, 63, 2);
    printf("Ensemble 2:\n");
    printf("\t%lx %lx\n", set[1], *set);
    printf("\t%lx %lx\n", setinvert[1], *setinvert);
    b = (*set == GMP_NUMB_MAX && set[1] == 0 && setinvert[1] == GMP_NUMB_MAX && *setinvert == 0);
    printf("%i\n", b);
    
    /*Ensemble {64}*/
    *set = GMP_NUMB_MAX ^ (GMP_NUMB_MAX >> 1);
    setinvert[1] = 1;
    set[1] = 0;
    *setinvert = 0;
    setappend(set, setinvert, 64, 2);
    printf("Ensemble 3:\n");
    printf("\t%lx %lx\n", set[1], *set);
    printf("\t%lx %lx\n", setinvert[1], *setinvert);
    b = (*set == (GMP_NUMB_MAX ^ (GMP_NUMB_MAX >> 1)) && set[1] == 1 && setinvert[1] == 1 && *setinvert == (GMP_NUMB_MAX ^ (GMP_NUMB_MAX >> 1)));
    printf("%i\n", b);
    
    return b;
}

void testScan() {
    partition_t partitionstruc;
    //schurNumberScanPartitionFromFile("/Users/rubis/Documents/Partitons_Schur/Partitions_sans_sommes_vrac/Partition_sans_sommes1.txt", &partitionstruc);
    schurNumberScanPartitionFromFile("/Users/rubis/Documents/Partitons_Schur/Partitions_sans_sommes_vrac/test_scan.txt", &partitionstruc);
    printPartition(partitionstruc.p, partitionstruc.n, partitionstruc.partition);
    printPartition(partitionstruc.p, 64*partitionstruc.limballoc, partitionstruc.partitioninvert);
}

int main(int argc, const char * argv[]) {
    testScan();
    return 0;
}
