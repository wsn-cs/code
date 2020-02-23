//
//  schurNumberRecursiveIteration.c
//  schurNumberRecursive
//
//  Created by rubis on 01/12/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include "schurNumberRecursiveBuild.h"
#include <stdio.h>


#define ADD_POINT(set, x) set[(x) / mp_bits_per_limb] |= ((unsigned long)1 << ((x) % mp_bits_per_limb))

#define DELETE_POINT(set, x) set[(x) / mp_bits_per_limb] ^= ((unsigned long)1 << ((x) % mp_bits_per_limb))

#define GET_POINT(set, x) (set[(x) / mp_bits_per_limb] & ((unsigned long)1 << ((x) % mp_bits_per_limb)))

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

char schurNumberRecursivePartitionComplete(unsigned long n, unsigned long l, unsigned int p, mp_limb_t **partition, mp_limb_t **reverse_partition) {
    /*Complète une partition sans-somme de [1, l] U {(i + 1)l + i, 1 ≤ i ≤ [(n-l)/(l+1)]} en une partition sans-somme de [1, n].*/
    unsigned int j;
    unsigned long x, xrem;
    unsigned int shift;
    char notSumFree;
    mp_size_t size;
    mp_limb_t *set, *set1, *set2, *set3;
    
    size = (n / mp_bits_per_limb) + 1;
    set1 = calloc(size, sizeof(mp_limb_t));
    set2 = calloc(size, sizeof(mp_limb_t));
    set3 = calloc(size, sizeof(mp_limb_t));
    
    x = l + 1;
    j = 0;
    
    while (x > l) {
        
        // Vérifier que x n'a pas déjà été placé.
        if (l && ((x / l) - (x % l) == 1)) {
            x++;
        }
        
        /*Essayer de placer x dans un ensemble*/
        notSumFree = 1;
        while (notSumFree && j < p) {
            // Regarder si x peut être placé dans l'ensemble j, ce qui n'est pas le cas uniquement lorsqu'il existe y dans partition[j] tel que (y + x) ou (y - x) y appartienne aussi.
            xrem = x;
            set = partition[j];
            
            while (xrem > 0) {
                shift = xrem % mp_bits_per_limb;
                if (!shift) {
                    shift = mp_bits_per_limb - 1;
                }
                
                if (xrem == x) {
                    mpn_rshift(set1, set, size, shift);     // set1 = set - x
                    mpn_lshift(set2, set, size, shift);    // set2 = set + x
                } else {
                    mpn_rshift(set1, set1, size, shift);     // set1 = set - x
                    mpn_lshift(set2, set2, size, shift);    // set2 = set + x
                }
                
                xrem -= shift;
            }
            
            printSet(n, set);
            printSet(n, set1);
            printSet(n, set2);
            mpn_and_n(set3, set, set1, size);
            if (mpn_zero_p(set3, size)) {
                mpn_and_n(set3, set, set2, size);
                notSumFree = !mpn_zero_p(set3, size);
            }
            j += notSumFree;
        }
        
        if (!notSumFree) {
            /*x est placé dans l'ensemble j*/
            ADD_POINT(partition[j], x);
            x++;
            if (x > n) {
                free(set1);
                free(set2);
                free(set3);
                return 1;
            }
            j = 0;
        } else {
            /*x ne peut être placé nul part*/
            x--;
            if ((x / l) - (x % l) == 1) {
                x--;
            }
            j = 0;
            while (j < p && !GET_POINT(partition[j], x)) {
                j++;
            }
            DELETE_POINT(partition[j], x);
            j++;
        }
        
    }

    free(set1);
    free(set2);
    free(set3);
    
    return 0;
}

char schurNumberRecursiveIteration(unsigned long *integer_stack, size_t stack_size, unsigned int p, mp_limb_t **partition, mp_limb_t **reverse_partition) {
    /*Cette fonction exécute une itération du processus récursif : elle parcourt, à partir d'une partition de l = integer_stack[0], l'intégralité des partitions sans-sommes de n = stack_size[1], et pour chacune d'elle, s'appelle récursivement pour construire les partitions sans-sommes subséquentes de stack_size[2].*/
    unsigned long l, n;
    unsigned long i, imax;
    unsigned int j, k;
    
    l = integer_stack[0];
    n = integer_stack[1];
    imax = (n - l) / (l + 1);
    
    mp_size_t size;
    mp_limb_t *mask, *reverse_mask;
    size = (n / mp_bits_per_limb) + 1;
    mask = calloc(size, sizeof(mp_limb_t));
    reverse_mask = calloc(size, sizeof(mp_limb_t));
    for (i=1; i<=l; i++) {
        ADD_POINT(mask, i);
    }
    for (i=1; i<=imax; i++) {
        ADD_POINT(mask, (i + 1) * l + i);
    }
    
    // Répartir d'abord les entiers (i + 1)l + i, 1 ≤ i ≤ [(n-l)/(l+1)] = imax parmi les p ensembles
    i = 1;
    j = 0;
    while (i > 0) {
        
        if (i > imax) {
            // Tous les entiers à placer en premier ayant été répartis, il convient à présent de construire une partition sans-somme de [1, n] tout entier grâce à la fonction schurNumberRecursivePartitionComplete.
            
            if (schurNumberRecursivePartitionComplete(n, l, p, partition, reverse_partition)) {
                // Une partition sans-somme de [1, n] a été construite.
                if ((stack_size == 2) || schurNumberRecursiveIteration(integer_stack + 1, stack_size - 1, p, partition, reverse_partition)) {
                    free(mask);
                    return 1;
                }
                // Revenir à la partition initiale.
                for (k=0; k<p; k++) {
                    mpn_and_n(partition[k], partition[k], mask, size);
                }
            }
            
            // Tenter maintenant de placer imax dans l'ensemble suivant
            i--;
            DELETE_POINT(partition[j], (i + 1) * l + i);
            j++;
            while (j == p && i > 0) {
                i--;
                // Trouver l'ensemble contenant (i + 1)l + i
                j = 0;
                while (!GET_POINT(partition[j], (i + 1) * l + i)) {
                    j++;
                }
                // Retirer (i + 1)l + i de cet ensemble
                DELETE_POINT(partition[j], (i + 1) * l + i);
                j++;
            }
            
        } else {
            ADD_POINT(partition[j], (i + 1) * l + i);
            // Incrémenter i.
            i++;
            j = 0;
        }
        
    }
    
    free(mask);
    
    return 0;
}
