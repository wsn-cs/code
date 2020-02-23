//
//  homogenous_sparse_polynomial.c
//  schurNumberGeometric
//
//  Created by rubis on 11/03/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include "homogenous_sparse_polynomial.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc/malloc.h>

sparse_polynomial * sparse_polynomial_alloc(size_t coeffalloc, mp_size_t indlimbsize, size_t degree) {
    //unsigned long i;
    //long *coeff;
    sparse_polynomial *poly;
    
    poly = malloc(sizeof(sparse_polynomial));
    poly->coeffalloc = coeffalloc;
    poly->coeffsize = 0;
    poly->degree = degree;
    poly->n = 0;
    
    poly->indlimbsize = indlimbsize;
    printf("%lu %lu %lu %lu\n", indlimbsize, coeffalloc, sizeof(mp_limb_t), coeffalloc * indlimbsize * sizeof(mp_limb_t));
    poly->indexes = calloc(coeffalloc, indlimbsize * sizeof(mp_limb_t));
    
    printf("%p %p %p\n", poly->indexes, poly->indexes + (coeffalloc - 1) * indlimbsize * sizeof(mp_limb_t), poly->indexes + coeffalloc * indlimbsize * sizeof(mp_limb_t));
    printf("%lu\n", malloc_size(poly->indexes));
    
    //coeff = calloc(coeffalloc, sizeof(long));
    //poly->coefficients = coeff;
    poly->coefficients = calloc(coeffalloc, sizeof(long));
    printf("%p\n", poly->coefficients);
    /*for (i=0; i<coeffalloc; i++) {
        mpq_init(*coeff);
        mpq_set_si(*coeff, 0, 1);
        coeff++;
    }*/
    
    return poly;
}


void sparse_polynomial_copy(sparse_polynomial *polyr, sparse_polynomial *polys) {
    size_t i;
    size_t coeffsize, coeffalloc;
    mp_size_t indlimbsize;
    mp_limb_t *inds, *indr;
    long *coeffs, *coeffr;
    
    coeffsize = polys->coeffsize;
    indlimbsize = polys->indlimbsize;
    coeffalloc = polyr->coeffalloc;
    inds = polys->indexes;
    coeffs = polys->coefficients;
    indr = polyr->indexes;
    coeffr = polyr->coefficients;
    if (coeffsize > coeffalloc) {
        coeffalloc *= 2;
        indr = realloc(indr, coeffalloc * indlimbsize * sizeof(mp_limb_t));
        coeffr = realloc(coeffr, coeffalloc * sizeof(long));
        polyr->indexes = indr;
        polyr->coefficients = coeffr;
    } else if (indlimbsize > polyr->indlimbsize) {
        indr = realloc(indr, coeffalloc * indlimbsize * sizeof(mp_limb_t));
    }
    
    for (i=0; i<coeffsize; i++) {
        mpn_copyi(indr, inds, indlimbsize);
        //mpq_set(*coeffr, *coeffs);
        *coeffr = *coeffs;
        inds += indlimbsize;
        indr += indlimbsize;
        coeffs++;
        coeffr++;
    }
    
    polyr->coeffsize = polys->coeffsize;
    polyr->degree = polys->degree;
    polyr->n = polys->n;
}

char sparse_polynomial_realloc(sparse_polynomial *poly, size_t coeffalloc, mp_size_t indlimbsize) {
    size_t i;
    size_t coeffsize;
    mp_size_t oldindlimbsize;
    mp_size_t diffindlimbsize;
    mp_limb_t *indexes;
    long *coefficients;
    mp_limb_t *oldindptr;
    mp_limb_t *newindptr;
    
    //printf("Réallocation:\n");
    //sparse_polynomial_print(poly, 2);
//    for (i=0; i< poly->coeffsize * indlimbsize; i++) {
//        printf("%lu\n", *(poly->indexes + i * sizeof(mp_limb_t)));
//    }
    indexes = realloc(poly->indexes, coeffalloc * indlimbsize * sizeof(mp_limb_t));
//    for (i=0; i< coeffalloc * indlimbsize; i++) {
//        printf("%lu\n", *(indexes + i * sizeof(mp_limb_t)));
//    }
    coefficients = realloc(poly->coefficients, coeffalloc * sizeof(long));
    if (indexes && coefficients) {
        coeffsize = poly->coeffsize;
        oldindlimbsize = poly->indlimbsize;
        diffindlimbsize = indlimbsize - oldindlimbsize;
        oldindptr = indexes + oldindlimbsize*(coeffsize - 1);
        newindptr = indexes + indlimbsize*(coeffsize - 1);
        /*Parcourir les coefficients dans le sens des indices décroissants*/
        for (i=0; i<coeffsize; i++) {
            mpn_zero(newindptr + oldindlimbsize, diffindlimbsize);
            mpn_copyi(newindptr, oldindptr, oldindlimbsize);
            oldindptr -= oldindlimbsize;
            newindptr -= indlimbsize;
        }
        //for (i=oldcoeffalloc; i<coeffalloc; i++) {
        //    mpq_init(coefficients[i]);
        //}
        poly->indexes = indexes;
        poly->coefficients = coefficients;
        poly->coeffalloc = coeffalloc;
        poly->indlimbsize = indlimbsize;
        //sparse_polynomial_print(poly, 3);
        return 1;
    } else if (indexes) {
        poly->indexes = indexes;
    } else if (coefficients) {
        poly->coefficients = coefficients;
    }
    printf("Reallocation failed\n");
    return 0;
}

void sparse_polynomial_free(sparse_polynomial *poly) {
    //size_t i;
    //size_t coeffalloc;
    //mpq_t *coeff;
    
    free(poly->indexes);
    
    //coeffalloc = poly->coeffalloc;
    //coeff = poly->coefficients;
    //for (i=0; i<coeffalloc; i++) {
    //    mpq_clear(*coeff);
    //    coeff++;
    //}
    free(poly->coefficients);
    
    free(poly);
}
