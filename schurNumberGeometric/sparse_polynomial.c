//
//  sparse_polynomial.c
//  schurNumberGeometric
//
//  Created by rubis on 11/03/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include <stdlib.h>
#include <gmp.h>

struct sparse_polynomial_struc {
    size_t coeffalloc;          // Nombre maximal de coefficients
    size_t coeffsize;           // Nombre de coefficients utilisés
    size_t degree;              // Degré total
    
    mp_size_t indlimballoc;     // Nombre de limbes alloués par indice
    mp_limb_t *indexes;         // Tableau contenant les indices ordonnés lexicographiquement des coefficients ≠ 0
    
    mp_size_t coefflimballoc;   // Nombre de limbes alloués par coefficient
    mp_size_t coefflimbsize;    // Nombre de limbes utilisés par coefficient
    mp_limb_t *coefficients;    // Tableau des coefficients
};

typedef struct sparse_polynomial_struc sparse_polynomial;

sparse_polynomial * sparse_polynomial_alloc(size_t coeffalloc, mp_size_t indlimballoc, mp_size_t coefflimballoc) {
    sparse_polynomial *poly;
    poly = malloc(sizeof(sparse_polynomial));
    poly->coeffalloc = coeffalloc;
    poly->degree = 0;
    
    poly->indlimballoc = indlimballoc;
    poly->coefflimballoc = coefflimballoc;
    poly->coefflimbsize = 1;
    
    poly->indexes = calloc(coeffalloc, indlimballoc * sizeof(mp_limb_t));
    poly->coefficients = calloc(coeffalloc, coefflimballoc * sizeof(mp_limb_t));
    
    return poly;
}

char sparse_polynomial_realloc(sparse_polynomial *poly, size_t coeffalloc, mp_size_t indlimballoc, mp_size_t coefflimballoc) {
    size_t i;
    size_t coeffsize;
    mp_size_t oldindlimballoc;
    mp_size_t oldcoefflimballoc;
    mp_limb_t *indexes;
    mp_limb_t *coefficients;
    mp_limb_t *oldindptr;
    mp_limb_t *oldcoeffptr;
    mp_limb_t *newindptr;
    mp_limb_t *newcoeffptr;
    
    indexes = realloc(poly->indexes, coeffalloc * indlimballoc * sizeof(mp_limb_t));
    if (indexes) {
        coefficients = realloc(poly->coefficients, coeffalloc * coefflimballoc * sizeof(mp_limb_t));
        if (coefficients) {
            coeffsize = poly->coeffsize;
            oldindlimballoc = poly->indlimballoc;
            oldcoefflimballoc = poly->coefflimballoc;
            oldindptr = poly->indexes + oldindlimballoc*(coeffsize - 1);
            oldcoeffptr = poly->coefficients + oldcoefflimballoc*(coeffsize - 1);
            newindptr = indexes + indlimballoc*(coeffsize - 1);
            newcoeffptr = coefficients + coefflimballoc*(coeffsize - 1);
            for (i=0; i<coeffsize; i++) {
                mpn_copyi(newindptr, oldindptr, oldindlimballoc);
                mpn_copyi(newcoeffptr, oldcoeffptr, oldcoefflimballoc);
                oldindptr -= oldindlimballoc;
                oldcoeffptr -= oldcoefflimballoc;
                newindptr -= indlimballoc;
                newcoeffptr -= coefflimballoc;
            }
            poly->indexes = indexes;
            poly->coefficients = coefficients;
            poly->coeffalloc = coeffalloc;
            poly->indlimballoc = indlimballoc;
            poly->coefflimballoc = coefflimballoc;
            return 1;
        }
    }
    return 0;
}

void sparse_polynomial_free(sparse_polynomial *poly) {
    free(poly->indexes);
    free(poly->coefficients);
    free(poly);
}

void sparse_polynomial_add(sparse_polynomial *polydst, sparse_polynomial *polysrc1, sparse_polynomial *polysrc2) {
    char b;
    mp_size_t i1;
    mp_size_t i2;
    mp_size_t coeffsize0;
    mp_size_t coeffsize1;
    mp_size_t coeffsize2;
    mp_limb_t *indst;
    mp_limb_t *coeffdst;
    mp_limb_t *inddsthigh;
    mp_limb_t * coeffdsthigh;
    mp_size_t indlimbsizemin;
    mp_size_t indlimbsizediff;
    size_t degdst;
    size_t coeffallocdst;
    mp_limb_t *indsrc1;
    mp_limb_t *coeffsrc1;
    size_t deg1;
    mp_limb_t *indsrc2;
    mp_limb_t *coeffsrc2;
    size_t deg2;
    mp_limb_t *indsrc2diff;
    
    coeffsize1 = polysrc1->coeffsize;
    indlimbsizemin = polysrc1->indlimbsize;
    coeffsize2 = polysrc2->coeffsize;
    indlimbsizediff = polysrc2->indlimbsize;
    if (indlimbsizemin > indlimbsizediff) {
        return sparse_polynomial_add(polydst, polysrc2, polysrc1);
    }
    indlimbsizediff -= indlimbsizemin;
    if (coeffsize1 > coeffsize2) {
        coeffsize0 = coeffsize1;
    } else {
        coeffsize0 = coeffsize2;
    }
    
    indst = polydst->indexes;
    coeffdst = polydst->coefficients;
    coeffallocdst = polydst->coeffalloc;
    degdst = 0;
    
    indsrc1 = polysrc1->indexes;
    coeffsrc1 = polysrc1->coefficients;
    deg1 = polysrc1->degree;
    i1 = 0;
    
    indsrc2 = polysrc2->indexes;
    coeffsrc2 = polysrc2->coefficients;
    deg2 = polysrc2->degree;
    i2 = 0;
    indsrc2diff = indsrc2 + indlimbsizediff;
    
    while (i1 < deg1 || i2 < deg2) {
        if (mpn_zero_p(indsrc2diff, indlimbsizediff)) {
            b = mpn_cmp(indsrc1, indsrc2, indlimbsizemin);
            if (b == 0) {
                /*Les deux indices sont identiques: il faut additionner*/
                if (degdst + 1 >= coeffallocdst) {
                    sparse_polynomial_realloc(polydst, 2*coeffallocdst, polydst->indlimballoc, polydst->coefflimballoc);
                }
                mpn_add_n(indst, coeffsrc1, coeffsrc2, coeffsize0);
                mpn_copyi(inddsthigh, coeffsrc2, coeffsizediff);
            } else if (b > 0) {
                
            } else {
                
            }
            
        } else {
            
        }
        
    }
}

void sparse_polynomial_reduce(sparse_polynomial *poly1, sparse_polynomial *poly2) {
    /*Réduis poly1 modulo poly2 en suivant l'ordre lexicographique.*/
}
