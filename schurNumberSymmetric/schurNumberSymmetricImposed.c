//
//  schurNumberSymmetricImposed.c
//  schurNumberSymmetric
//
//  Created by rubis on 16/11/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include <stdio.h>
#include "schurNumberSymmetricHeader.h"
#include "schurNumberSymmetricInitialPartition.h"


#define ADD_POINT(set, x) set[x / mp_bits_per_limb] |= ((unsigned long)1 << (x % mp_bits_per_limb))

#define DELETE_POINT(set, x) set[x / mp_bits_per_limb] ^= ((unsigned long)1 << (x % mp_bits_per_limb))

#define GET_POINT(set, x) set[x / mp_bits_per_limb] & ((unsigned long)1 << (x % mp_bits_per_limb))

unsigned long schurNumberSymmetricInitialPartition4[4][12] = SCHUR_INITIAL_PARTITION_4;
unsigned long schurNumberSymmetricLengthsInitialPartition4[4] = SCHUR_LENGTHS_INITIAL_PARTITION_4;

unsigned long schurNumberSymmetricInitialPartition5[5][40] = SCHUR_INITIAL_PARTITION_5;
unsigned long schurNumberSymmetricLengthsInitialPartition5[5] = SCHUR_LENGTHS_INITIAL_PARTITION_5;

unsigned long schurNumberSymmetricBuildInitialPartition(unsigned long n, unsigned long p, mp_limb_t **partition) {
    /* Cette fonction attribue aux ensembles de partition leur valeur initiale, qui forment une partition sans-somme à p ensembles. Elle renvoie la profondeur de la (p+1)-ième partie.*/
    
    unsigned int i, j, k;
    unsigned long *lengths;
    unsigned long *values;
    unsigned int array_size;
    unsigned long depth;
    
    switch (p-1) {
        case 4:
            lengths = schurNumberSymmetricLengthsInitialPartition4;
            values = schurNumberSymmetricInitialPartition4;
            array_size = 12;
            break;
            
        case 5:
            lengths = schurNumberSymmetricLengthsInitialPartition5;
            values = schurNumberSymmetricInitialPartition5;
            array_size = 40;
            break;
            
        default:
            ADD_POINT(partition[1], 1);
            ADD_POINT(partition[1], (n - 1));
            return 1;
    }
    
    k = 0;
    
    for (i=0; i<p-1; i++) {
        for (j=0; j < lengths[i]; j++) {
            k = array_size * i + j;
            //printf("%u : %lu\n", k, values[k]);
            ADD_POINT(partition[i], values[k]);
            ADD_POINT(partition[i], (n - values[k]));
        }
    }
    
    depth = values[*lengths - 1] + 1;
    ADD_POINT(partition[i], depth); // désormais i = p-1
    ADD_POINT(partition[i], (n - depth));
    
    for (j = n - 2 * depth + 1; j < 2 * depth; j++) {
        ADD_POINT(partition[i], j);
    }
    
    return depth;
}

void buildInitialPartition5(unsigned long n, unsigned long x, mp_limb_t **partition) {
    /* La fonction initialise les ensembles de partition et leur attribue leur valeur initiale*/
    
    unsigned int i, j;
    mp_size_t size;
    /*unsigned long values[5][34] = {{1, 5, 8, 12, 15, 19, 25, 29, 32, 36, 39, 43, 118, 122, 125, 129, 132, 136, 142, 146, 149, 153, 156, 160},
        {2, 6, 7, 17, 18, 26, 27, 37, 38, 42, 119, 123, 124, 134, 135, 143, 144, 154, 155, 159},
        {3, 4, 9, 10, 11, 16, 28, 33, 34, 35, 40, 41, 120, 121, 126, 127, 128, 133, 145, 150, 151, 152, 157, 158},
        {13, 14, 20, 21, 22, 23, 24, 30, 31, 130, 131, 137, 138, 139, 140, 141, 147, 148},
        {44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117}};*/
    //unsigned long lengths[5] = {24, 20, 24, 18, 34};
    
    /*unsigned long values[5][44] = {{1, 5, 8, 12, 15, 19, 25, 29, 32, 36, 39, 43, 118, 122, 125, 129, 132, 136, 142, 146, 149, 153, 156, 160},
        {2, 6, 7, 17, 18, 26, 27, 37, 38, 42, 119, 123, 124, 134, 135, 143, 144, 154, 155, 159},
        {3, 4, 9, 10, 11, 16, 28, 33, 34, 35, 40, 41, 120, 121, 126, 127, 128, 133, 145, 150, 151, 152, 157, 158},
        {13, 14, 20, 21, 22, 23, 24, 30, 31, 130, 131, 137, 138, 139, 140, 141, 147, 148},
        {74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87}
    };*/
    //unsigned long lengths[5] = {24, 20, 24, 18, 132 - 2*x - 1};
    
    unsigned long values[5][24] = {{1, 5, 8, 12, 15, 19, 25, 29, 32, 36, 39, 43, 118, 122, 125, 129, 132, 136, 142, 146, 149, 153, 156, 160},
        {2, 6, 7, 17, 18, 26, 27, 37, 38, 42, 119, 123, 124, 134, 135, 143, 144, 154, 155, 159},
        {3, 4, 9, 10, 11, 16, 28, 33, 34, 35, 40, 41, 120, 121, 126, 127, 128, 133, 145, 150, 151, 152, 157, 158},
        {13, 14, 20, 21, 22, 23, 24, 30, 31, 130, 131, 137, 138, 139, 140, 141, 147, 148},
        {44, 117}
    };
    unsigned long lengths[5] = {24, 20, 24, 18, 2};
    
    /*for (i=44; i <= x; i++) {
        values[4][i-44 + 15] = i;
        values[4][i-44 + x-44 + 16] = n - i;
    }
    
    for (i=59; i <= 160 - 2*x; i++) {
        values[4][i-59 + 2*(x-44) + 17] = i;
        values[4][i-59 + (160 - 2*x - 59) + 2*(x-44) + 18] = n - i;
    }*/
    
    size = (n / mp_bits_per_limb) + 1;
    
    for (i=0; i<5; i++) {
        //partition[i] = calloc(size, sizeof(mp_limb_t));
        
        for (j=0; j < lengths[i]; j++) {
//            printf("%lu = %lu * %i + %lu\n", values[i][j], values[i][j]/mp_bits_per_limb, mp_bits_per_limb, values[i][j]%mp_bits_per_limb);
            ADD_POINT(partition[i], values[i][j]);
        }
    }
}

void buildInitialPartition6(unsigned long n, unsigned long x, mp_limb_t **partition) {
    /* La fonction initialise les ensembles de partition et leur attribue leur valeur initiale*/
    
    unsigned int i, j;
    mp_size_t size;
    
    unsigned long values[5][40] = {{1, 5, 8, 12, 15, 19, 25, 29, 32, 36, 39, 43, 45, 49, 52, 56, 105,
        109, 112, 116, 118, 122, 125, 129, 132, 136, 142, 146, 149, 153,
        156, 160},
        {2, 6, 7, 17, 18, 26, 27, 37, 38, 42, 46, 47, 50, 62, 66,
            95, 99, 111, 114, 115, 119, 123, 124, 134, 135, 143, 144, 154, 155,
            159},
        {3, 4, 9, 10, 11, 16, 28, 33, 34, 35, 40, 41, 53, 58, 59, 60,
                65, 96, 101, 102, 103, 108, 120, 121, 126, 127, 128, 133, 145, 150,
                151, 152, 157, 158},
        {13, 14, 20, 21, 22, 23, 24, 30, 31, 63, 69,
                    73, 88, 92, 98, 130, 131, 137, 138, 139, 140, 141, 147, 148},
        {44, 48, 51, 54, 55, 57, 61, 64, 67, 68, 70, 71, 72, 74, 75, 76, 77, 78,
                        79, 80, 81, 82, 83, 84, 85, 86, 87, 89, 90, 91, 93, 94, 97, 100,
                        104, 106, 107, 110, 113, 117}};
    unsigned long lengths[5] = {32, 30, 34, 24, 40};
    
    size = (n / mp_bits_per_limb) + 1;
    
    for (i=0; i<6; i++) {
        //partition[i] = calloc(size, sizeof(mp_limb_t));
        
        for (j=0; j < lengths[i]; j++) {
            //            printf("%lu = %lu * %i + %lu\n", values[i][j], values[i][j]/mp_bits_per_limb, mp_bits_per_limb, values[i][j]%mp_bits_per_limb);
            ADD_POINT(partition[i], values[i][j]);
            ADD_POINT(partition[i], (n - values[i][j]));
        }
    }
}

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
                    printf(" %lu", limbn*mp_bits_per_limb + j);
                }
            }
        }
        printf("\n");
    }
}

void printSet(unsigned long n, mp_limb_t *set) {
    /*Affiche un ensemble.*/
    unsigned long limbn;
    mp_size_t limbsize;
    mp_bitcnt_t j;
    mp_limb_t limb;
    limbsize = (n>>6) + 1;
    printf("Ensemble:\n");
    for (limbn=0; limbn<limbsize; limbn++) {
        limb = set[limbn];
        for (j=0; j<mp_bits_per_limb; j++) {
            if (limb & ((unsigned long)1<<j)) {
                printf(" %lu", limbn*mp_bits_per_limb + j);
            }
        }
    }
    printf("\n");
}

unsigned long schurNumberSymmetricImposedPartition(unsigned long n, unsigned long p) {
    
    mp_limb_t **partition;
    mp_limb_t *set, *set1, *set2, *set3;
    mp_limb_t overflow;
    mp_size_t size;
    unsigned long x, y, xmax;
    unsigned long depth, surface;
    unsigned long xrem;
    unsigned int shift;
    unsigned int i;
    char notSumFree;
    
    size = (n / mp_bits_per_limb) + 1;
    set1 = calloc(size, sizeof(mp_limb_t));
    set2 = calloc(size, sizeof(mp_limb_t));
    set3 = calloc(size, sizeof(mp_limb_t));
    partition = calloc(p, sizeof(mp_limb_t *));
    for (i=0; i<p; i++) {
        partition[i] = calloc(size, sizeof(mp_limb_t));
    }
    depth = schurNumberSymmetricBuildInitialPartition(n, p, partition);
    surface = n - 2 * depth + 1;
    printPartition(p, n-1, partition);
    
    x = depth + 1;
    xmax = x;
    i = 0;
    
    while ((depth < x) && (x < surface)) {
        /*Essayer de placer x dans un ensemble*/
        notSumFree = 1;
        while (notSumFree && i < p) {
            /*Regarder si x et n - x peuvent être placés dans l'ensemble i.
             Le cas dès qu'il existe y dans partition[i] tel que (y + x) ou (y - x) y appartienne aussi.*/
            
            /*while (GET_POINT(partition[4], x)) {
                x++;
            }*/
            xrem = x;
            set = partition[i];
            mpn_copyi(set1, set, size); // set1 = set - x
            mpn_copyi(set2, set, size); // set2 = set + x
            
            while (xrem > 0) {
                shift = xrem % mp_bits_per_limb;
                if (!shift) {
                    shift = mp_bits_per_limb - 1;
                }
                //printSet(n, set1);
                //overflow = set1[0] & (GMP_NUMB_MASK >> (mp_bits_per_limb - shift));
                overflow = mpn_rshift(set1, set1, size, shift);
                for (y=0; y < shift; y++) {
                    if (GET_POINT((&overflow), y)) {
                        ADD_POINT(set1, ((y - mp_bits_per_limb) + n));
                        //printf("%lu\n", (y - mp_bits_per_limb) + n);
                        //printSet(n, set1);
                    }
                }
                //printSet(n, set1);
//                set1[0] |= overflow;
                //printSet(n, set2);
                
                overflow = mpn_lshift(set2, set2, size, shift);
                for (y=0; y < shift; y++) {
                    if (GET_POINT((&overflow), y)) {
                        ADD_POINT(set2, ((y + size * mp_bits_per_limb) - n));
                        //printf("%lu %lu\n", y, (y + size * mp_bits_per_limb) - n);
                    }
                }
                for (y=n; y < n + mp_bits_per_limb; y++) {
                    if (GET_POINT(set2, y)) {
                        DELETE_POINT(set2, y);
                        ADD_POINT(set2, (y - n));
                        //printf("%lu\n", y);
                    }
                }
                //printSet(n, set2);
                //set2[size - 1] |= overflow;
                xrem -= shift;
            }
            
            mpn_and_n(set3, set, set1, size);
            if (mpn_zero_p(set3, size)) {
                mpn_and_n(set3, set, set2, size);
                notSumFree = !mpn_zero_p(set3, size);
            }
            i += notSumFree;
        }
        
        if (!notSumFree && x <= surface) {
            /*x et n - x sont placés dans l'ensemble i*/
            ADD_POINT(partition[i], x);
            ADD_POINT(partition[i], (n - x));
            x++;
            if (x > xmax) {
                xmax = x;
                printPartition(p, n, partition);
            }
            /*if (x > n / 2) {
                printPartition(5, n, partition);
            }*/
            i = 0;
        } else {
            /*x ne peut être placé nul part*/
            x--;
            /*while (GET_POINT(partition[4], x)) {
                x--;
            }*/
            i = 0;
            while (i < p && !(GET_POINT(partition[i], x))) {
                i++;
            }
            DELETE_POINT(partition[i], x);
            DELETE_POINT(partition[i], (n - x));
            i++;
        }
        
    }
    
    printPartition(p, n-1, partition);
    printf("%lu\n", xmax);
    
    /*Nettoyage*/
    for (i=0; i<p; i++) {
        free(partition[i]);
    }
    free(partition);
    
    free(set1);
    free(set2);
    free(set3);
    
    return x;
}
