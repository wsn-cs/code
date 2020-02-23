//
//  schurNumberRecursiveBuild.c
//  schurNumberRecursive
//
//  Created by rubis on 01/12/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include "schurNumberRecursiveBuild.h"
#include <math.h>
#include <stdio.h>

#define ADD_POINT(set, x) set[(x) / mp_bits_per_limb] |= ((unsigned long)1 << ((x) % mp_bits_per_limb))

#define DELETE_POINT(set, x) set[(x) / mp_bits_per_limb] ^= ((unsigned long)1 << ((x) % mp_bits_per_limb))

#define GET_POINT(set, x) (set[(x) / mp_bits_per_limb] & ((unsigned long)1 << ((x) % mp_bits_per_limb)))

unsigned long * indermediateIntegerStack(unsigned long n, size_t *stack_size_ptr, int recursion_depth) {
    /*Cette fonction construit la pile d'entiers l1 < l2 < … utile pour la construction récursive de la partition sans-somme, et la renvoie.*/
    
    unsigned long *stack;
    
    if (!n) {
        stack = calloc(recursion_depth + 1, sizeof(unsigned long));
        stack[0] = n;
        *stack_size_ptr = recursion_depth + 1;
        return stack;
    }
    
    unsigned long l;
    
    l = lround(sqrt(n + 1) - 1);
    if (l + ((n - l) / (l+1)) > l - 1 + ((n - l + 1) / (l+1))) {
        l--;
    }
    /*if (!l) {
        l = 1;
    }*/
    
    stack = indermediateIntegerStack(l, stack_size_ptr, recursion_depth + 1);
    stack[*stack_size_ptr - recursion_depth - 1] = n;
    
    return stack;
}

char schurNumberRecursiveBuild(unsigned long n, unsigned int p) {
    /*Cette fonction essaye de construire une partition sans-somme de [1, n] à p ensembles. Elle renvoie 1 si une telle partition existe et 0 sinon.
     Elle fait appel à la fonction schurNumberRecursiveIteration.*/
    
    unsigned int i;
    char sum_free_partition_found;
    
    //Détermination de la taille des ensembles en limbs
    mp_size_t limb_size = (n / mp_bits_per_limb) + 1;
    
    // Création de la pile d'entiers intermédiaires
    unsigned long *integer_stack;
    size_t stack_size;
    integer_stack = indermediateIntegerStack(n, &stack_size, 0);
    
    for (i=0; i<stack_size; i++) {
        printf("%lu\n", integer_stack[i]);
    }
    
    // Allocation de la partition
    mp_limb_t **partition;
    mp_limb_t **reverse_partition;
    partition = calloc(p, sizeof(mp_limb_t *));
    reverse_partition = calloc(p, sizeof(mp_limb_t *));
    for (i=0; i<p; i++) {
        partition[i] = calloc(limb_size, sizeof(mp_limb_t));
        reverse_partition[i] = calloc(limb_size, sizeof(mp_limb_t));
    }
    
    sum_free_partition_found = schurNumberRecursiveIteration(integer_stack, stack_size, p, partition, reverse_partition);
    printPartition(p, n, partition);
    
    // Libération de la partition
    for (i=0; i<p; i++) {
        free(partition[i]);
        free(reverse_partition[i]);
    }
    free(partition);
    free(reverse_partition);
    free(integer_stack);
    
    return sum_free_partition_found;
}
