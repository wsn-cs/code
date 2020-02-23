//
//  homogenous_sparse_polynomial_mul.c
//  schurNumberGeometric
//
//  Created by rubis on 16/03/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include "homogenous_sparse_polynomial.h"
#include </Users/rubis/Documents/Partitons_Schur/SchurNumber/SchurNumber/rscan1.c>

void sparse_polynomial_mul_monom(sparse_polynomial *polyr, sparse_polynomial *polys, mp_limb_t *monomindex, mp_size_t monomlimbsize, long monomcoeff) {
    /*Effectue polyr = monomcoeff.X^(monomindex).polys
     Suppose que polyr->indlimbsize = polys->indlimbsize = monomlimbsize*/
    size_t i, coeffsize, coeffalloc;
    mp_limb_t *inds, *indr;
    long *coeffs, coeffr;
    
    /*Vérifie si polyr peut contenir le résultat*/
    coeffsize = polys->coeffsize;
    coeffalloc = polyr->coeffalloc;
    if (coeffsize > coeffalloc) {
        coeffalloc *= 2;
        sparse_polynomial_realloc(polyr, coeffalloc, monomlimbsize);
    }
    
    /*Initialisation pour polys*/
    inds = polys->indexes;
    coeffs = polys->coefficients;
    
    /*Initialisation pour polyr*/
    indr = calloc(monomlimbsize, sizeof(mp_limb_t));
    polyr->coeffsize = 0;
    
    /*Détermination du degré du monôme
    if (mpn_zero_p(monomindex, monomlimbsize)) {
        polyr->degree = polys->degree;
    } else {
        polyr->degree = polys->degree + mpn_rscan1(monomindex, monomlimbsize);
    }*/
    
    for (i=0; i<coeffsize; i++) {
        mpn_ior_n(indr, inds, monomindex, monomlimbsize);
        coeffr = *coeffs * monomcoeff;
        sparse_polynomial_add_monom(polyr, polyr, indr, monomlimbsize, coeffr);
        inds += monomlimbsize;
        coeffs++;
    }
    
    /*Déallocation*/
    free(indr);
}

void sparse_polynomial_add_mul_monom(sparse_polynomial *polyr, sparse_polynomial *polys1, sparse_polynomial *polys2, mp_limb_t *monomindex, mp_size_t monomlimbsize, long monomcoeff) {
    /*Effectue polyr = polys1 + monomcoeff.X^(monomindex).polys2
     Suppose que polyr->indlimbsize = polys1->indlimbsize = polys2->indlimbsize = monomlimbsize.
     Le pointeur polys2 doit être distinct de polyr.*/
    size_t i, coeffsize;
    mp_limb_t *inds2, *indm;
    long *coeffs2, coeffm;
    sparse_polynomial wpolyr;
    
    /*Initialisation pour polys2*/
    coeffsize = polys2->coeffsize;
    inds2 = polys2->indexes;
    coeffs2 = polys2->coefficients;
    
    /*Initialisation pour les monômes à ajouter successivement*/
    indm = calloc(monomlimbsize, sizeof(mp_limb_t));
    //mpq_init(coeffm);
    
    /*Initialisation de wpolyr, copie de polyr*/
    wpolyr.coeffalloc = polyr->coeffalloc;
    wpolyr.coeffsize = 0;
    wpolyr.degree = polys1->degree;
    wpolyr.n = 0;
    wpolyr.indlimbsize = monomlimbsize;
    wpolyr.indexes = polyr->indexes;
    wpolyr.coefficients = polyr->coefficients;
    
    /*Première étape pour initialiser polyr*/
    mpn_ior_n(indm, inds2, monomindex, monomlimbsize);
    coeffm = *coeffs2 * monomcoeff;
    sparse_polynomial_add_monom(&wpolyr, polys1, indm, monomlimbsize, coeffm);
    inds2 += monomlimbsize;
    coeffs2++;
    
    for (i=1; i<coeffsize; i++) {
        mpn_ior_n(indm, inds2, monomindex, monomlimbsize);
        coeffm = *coeffs2 * monomcoeff;
        sparse_polynomial_add_monom(&wpolyr, &wpolyr, indm, monomlimbsize, coeffm);
        inds2 += monomlimbsize;
        coeffs2++;
    }
    
    /*Propriétés remises à polyr*/
    polyr->coeffalloc = wpolyr.coeffalloc;
    polyr->coeffsize = wpolyr.coeffsize;
    polyr->degree = wpolyr.degree;
    polyr->n = wpolyr.n;
    polyr->indexes = wpolyr.indexes;
    polyr->coefficients = wpolyr.coefficients;
    
    /*Déallocation*/
    free(indm);
    //mpq_clear(coeffm);
}

void sparse_polynomial_mul(sparse_polynomial *polyr, sparse_polynomial *polys1, sparse_polynomial *polys2) {
    /*Effectue polyr = polys1.polys2
     Suppose que polyr->indlimbsize = polys1->indlimbsize = polys2->indlimbsize.
     Les pointeurs polys1 et polys2 doivent être distincts de polyr.*/
    size_t coeffsize, i;
    mp_size_t indlimbsize;
    mp_limb_t *ind2;
    long *coeff2;
    sparse_polynomial *swapoly;
    
    coeffsize = polys2->coeffsize;
    if (coeffsize > polys1->coeffsize) {
        /*Intervertit polys1 et polys2*/
        swapoly = polys2;
        polys2 = polys1;
        polys1 = swapoly;
        coeffsize = polys2->coeffsize;
    }
    indlimbsize = polys2->indlimbsize;
    ind2 = polys2->indexes;
    coeff2 = polys2->coefficients;
    
    /*Initialisation de polyr*/
    polyr->coeffsize = 0;
    
    for (i=0; i<coeffsize; i++) {
        sparse_polynomial_add_mul_monom(polyr, polyr, polys1, ind2, indlimbsize, *coeff2);
        ind2 += indlimbsize;
        coeff2++;
    }
}
