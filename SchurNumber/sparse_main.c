//
//  sparse_main.c
//  SchurNumber
//
//  Created by rubis on 24/02/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include </Users/rubis/Downloads/gmp-6.1.2/gmp.h>

typedef struct mpz_structlist {
    unsigned long size;
    unsigned long count;
    mpz_t *array;
} mpz_list;

void mpz_list_alloc(mpz_list *l, unsigned long size) {
    l = calloc(1, sizeof(mpz_list));
    l->size = size;
    l->count = 0;
    l->array = calloc(size, sizeof(mpz_t));
}

void mpz_list_realloc(mpz_list *l, unsigned long size) {
    l->size = size;
    l->array = realloc(l->array, size);
}

void mpz_list_free(mpz_list *l) {
    unsigned long i;
    for (i=0; i<l->count; i++) {
        mpz_clear(l->array[i]);
    }
    free(l);
}

void mpz_list_append(mpz_list *l, mpz_t a) {
    unsigned long i;
    i = l->count;
    if (i >= l->size) {
        mpz_list_realloc(l, 2*l->size);
    }
    mpz_init_set(l->array[i], a);
    l->count = i+1;
}

char mpz_list_contain(mpz_list *l, mpz_t a) {
    /* Renvoie 1 si a appartient à l, 0 sinon.
     La liste l est supposée ordonnée. */
    unsigned long i0;
    unsigned long i1;
    unsigned long i2;
    char notFound;
    mpz_t *array;
    i1 = 0;
    i2 = l->count;
    i0 = i2;
    notFound = 1;
    array = l->array;
    while (notFound || i1 != i0) {
        i0 = (i2 - i1)>>1;
        notFound = mpz_cmp(a, array[i0]);
        if (notFound > 0) {
            // a > array[i]
            i1 = i0;
        }
        if (notFound < 0) {
            // a < array[i]
            i2 = i0;
        }
    }
    return !notFound;
}

void sparseSumFreeStep(mpz_list *sfi, unsigned long n, mp_limb_t *work) {
    unsigned long i;
    unsigned long c;
    unsigned long limbnum;
    mpz_ptr a;
    mpz_t part;
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t *work2;
    c = sfi->count;
    limbnum = (n>>7)+1;
    work0 = calloc(limbnum, sizeof(mp_limb_t));
    work1 = calloc(limbnum, sizeof(mp_limb_t));
    work2 = calloc(limbnum, sizeof(mp_limb_t));
    for (i=0; i<c; i++) {
        a = sfi->array[i];
        mpn_copyd(work1, a->_mp_d, limbnum);
        mpn_copyd(work2, &(a->_mp_d[limbnum]), limbnum);
        mpn_and_n(work0, work1, work2, limbnum);
        if (mpn_zero_p(work0, limbnum)) {
            mpz_init2(part, n);
            mpz_set(part, a);
            mpz_setbit(part, n);
            mpz_list_append(sfi, part);
        }
    }
    free(work0);
    free(work1);
    free(work2);
}

unsigned long schurNumber(unsigned long pmax) {
    unsigned long n;
    unsigned long i;
    unsigned long p;
    unsigned long nbest;
    char shouldContinue;
    char isSumFree;
    mpz_t a;
    mpz_list *sf = NULL;
    mpz_t *sfpartition;
    n = pmax;
    nbest = p;
    mpz_list_alloc(sf, 32);
    sfpartition = calloc(p, sizeof(mpz_t));
    for (i=0; i<p; i++) {
        mpz_init_set_ui(a, 1<<i);
        mpz_init_set_ui(sfpartition[i], 1<<i);
        mpz_list_append(sf, a);
    }
    
    /*Itération jusqu'à trouver S(p)*/
    shouldContinue = 1;
    i = 0;
    p = 1;
    while (shouldContinue) { //Arrêt si sfpartition[p-1]=={p,p+1,…,n}
        // Placer n+1 dans une des huches en conservant l'inadditivité
        while (i < p) {
            mpz_setbit(sfpartition[i], n); // Ajoute n à la huche i
            isSumFree = mpz_list_contain(sf, sfpartition[i]);
            if (isSumFree) {
                n++;
                if (n+p > nbest) {
                    sparseSumFreeStep(sf, n, NULL);
                    nbest = n+p;
                }
                i = 0;
            } else {
                mpz_clrbit(sfpartition[i], n); // Retire n de la huche i
                i++;
            }
        }
        // Dépiler
        n --;
        i = 0;
        while (i < p && !mpz_tstbit(sfpartition[i], n)) {
            i++;
        }
        mpz_clrbit(sfpartition[i], n);
        i++;
    }
    
    mpz_list_free(sf);
    return nbest;
}

int main(int argc, const char * argv[]) {
    mpz_t b;
    mpz_init_set_ui(b, 7);
    printf("%i ", b->_mp_size);
    mpz_clrbit(b, 2); // b = 3
    gmp_printf("%Zd ", b);
    printf("%i\n", b->_mp_size);
    mpz_clear(b);
    return 0;
}

