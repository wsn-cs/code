//
//  homogenous_sparse_polynomial_inversible.c
//  schurNumberGeometric
//
//  Created by rubis on 16/03/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include "homogenous_sparse_polynomial.h"
#include <stdio.h>

void mpq_homogenous_sparse_polynomial_div_monom(mpq_homogenous_sparse_polynomial *polyr, mpq_homogenous_sparse_polynomial *polys, mp_limb_t *monomindex, mp_size_t monomlimbsize) {
    /*Effectue polyr = polys(Xn = X0)
     Suppose que polyr->indlimbsize = polys->indlimbsize = monomlimbsize.
     Le pointeur polys doit être distinct de polyr.*/
    size_t i, coeffsize;
    mp_limb_t *inds, *indm;
    mpq_t *coeffs;
    mpq_homogenous_sparse_polynomial wpolyr;
    
    /*Initialisation pour polys*/
    coeffsize = polys->coeffsize;
    inds = polys->indexes;
    coeffs = polys->coefficients;
    
    /*Initialisation pour les monômes à ajouter successivement*/
    indm = calloc(monomlimbsize, sizeof(mp_limb_t));
    
    /*Initialisation de wpolyr, copie de polyr*/
    wpolyr.coeffalloc = polyr->coeffalloc;
    wpolyr.coeffsize = 0;
    wpolyr.degree = polys->degree;
    wpolyr.n = 0;
    wpolyr.indlimbsize = monomlimbsize;
    wpolyr.indexes = polyr->indexes;
    wpolyr.coefficients = polyr->coefficients;
    
    for (i=0; i<coeffsize; i++) {
        mpn_andn_n(indm, inds, monomindex, monomlimbsize);
        mpq_homogenous_sparse_polynomial_add_monom(&wpolyr, &wpolyr, indm, monomlimbsize, coeffs);
        inds += monomlimbsize;
        coeffs++;
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
}

void mpq_homogenous_sparse_polynomial_add_div_monom(mpq_homogenous_sparse_polynomial *polyr, mpq_homogenous_sparse_polynomial *polys1, mpq_homogenous_sparse_polynomial *polys2, mp_limb_t *monomindex, mp_size_t monomlimbsize) {
    /*Effectue polyr = X0.polys1 + polys2(Xn = X0)
     Suppose que polyr->indlimbsize = polys1->indlimbsize = polys2->indlimbsize = monomlimbsize.
     Le pointeur polys2 doit être distinct de polyr.*/
    size_t i, coeffsize;
    mp_limb_t *inds2, *indm;
    mpq_t *coeffs2;
    mpq_homogenous_sparse_polynomial wpolyr;
    
    /*Initialisation pour polys2*/
    coeffsize = polys2->coeffsize;
    inds2 = polys2->indexes;
    coeffs2 = polys2->coefficients;
    
    /*Initialisation pour les monômes à ajouter successivement*/
    indm = calloc(monomlimbsize, sizeof(mp_limb_t));
    
    /*Initialisation de wpolyr, copie de polyr*/
    wpolyr.coeffalloc = polyr->coeffalloc;
    wpolyr.coeffsize = 0;
    wpolyr.degree = polys2->degree;
    wpolyr.n = 0;
    wpolyr.indlimbsize = monomlimbsize;
    wpolyr.indexes = polyr->indexes;
    wpolyr.coefficients = polyr->coefficients;
    
    /*Première étape pour initialiser polyr*/
    mpn_andn_n(indm, inds2, monomindex, monomlimbsize);
    mpq_homogenous_sparse_polynomial_add_monom(&wpolyr, polys1, indm, monomlimbsize, coeffs2);
    inds2 += monomlimbsize;
    coeffs2++;
    
    for (i=1; i<coeffsize; i++) {
        mpn_andn_n(indm, inds2, monomindex, monomlimbsize);
        mpq_homogenous_sparse_polynomial_add_monom(&wpolyr, &wpolyr, indm, monomlimbsize, coeffs2);
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
}

char mpq_homogenous_sparse_polynomial_inversible(mpq_homogenous_sparse_polynomial *poly, mpq_homogenous_sparse_polynomial *polyinv, mpq_homogenous_sparse_polynomial *wpolyrem, mpq_homogenous_sparse_polynomial *wpolyquot, mpq_homogenous_sparse_polynomial *wpoly3, mpq_homogenous_sparse_polynomial *wpoly4, mpq_homogenous_sparse_polynomial *wpoly5) {
    /*Renvoie 1 si il existe polyinv tel que poly.polyinv = 1*/
    char isinv;
    size_t coeffsize, coeffsize2;
    size_t degree;
    size_t nmax, n, ndiv64, nmod64;
    mp_size_t indlimbsize;
    mpq_homogenous_sparse_polynomial poly1, poly2;
    mp_limb_t *nind, *ind2;
    mp_limb_t mask;
    mpq_t one;
    
    /*Initialisation des variables*/
    coeffsize = poly->coeffsize;
    degree = poly->degree;
    nmax = poly->n;
    indlimbsize = poly->indlimbsize;
    
    /*Inversion du coefficient constant*/
    if (!coeffsize || !mpn_zero_p(poly->indexes, indlimbsize)) {
        /*Coefficient constant nul*/
        return 0;
    }
    polyinv->coeffsize = 1;
    polyinv->degree = 1;
    mpn_zero(polyinv->indexes, indlimbsize);
    mpq_inv(*polyinv->coefficients, *poly->coefficients);
    
    /*Allocation*/
    mpq_init(one);
    mpq_set_si(one, -1, 1);
    nind = calloc(indlimbsize, sizeof(mp_limb_t));
    
    /*Initialisation de poly1, égal successivement à poly(X0,…,X(n-1),0,…)*/
    poly1.coeffalloc = poly->coeffalloc;
    poly1.coeffsize = 1;
    poly1.degree = degree;
    poly1.n = 0;
    poly1.indlimbsize = indlimbsize;
    poly1.indexes = poly->indexes;
    poly1.coefficients = poly->coefficients;
    
    /*Initialisation de poly2, égal à poly - poly1*/
    poly2.coeffalloc = poly->coeffalloc;
    poly2.coeffsize = 1;
    poly2.degree = degree;
    poly2.n = 1;
    poly2.indlimbsize = indlimbsize;
    poly2.indexes = poly->indexes + indlimbsize;
    poly2.coefficients = poly->coefficients + 1;
    
    n = 1;
    *nind = (mp_limb_t)1;
    isinv = 1;
    while (isinv && n <= nmax) {
        /*Détermination de la taille de poly2*/
        coeffsize2 = 0;
        ind2 = poly2.indexes;
        ndiv64 = (n-1)>>6;
        nmod64 = (n-1)%64;
        mask = 1<<nmod64;
        while (ind2[ndiv64] & mask) {
            ind2 += indlimbsize;
            coeffsize2++;
        }
        poly2.coeffsize = coeffsize2;
        //wpoly3 = poly2(Xn = X0) + X0.poly1
        printf("\npoly1:\n");
        mpq_homogenous_sparse_polynomial_print(&poly1, 2);
        printf("\npoly2:\n");
        mpq_homogenous_sparse_polynomial_print(&poly2, 2);
        //mpq_homogenous_sparse_polynomial_add_mul_monom(wpoly3, &poly2, &poly1, nind, indlimbsize, &one);
        mpq_homogenous_sparse_polynomial_add_div_monom(wpoly3, &poly1, &poly2, nind, indlimbsize);
        printf("\npoly2(Xn = X0) + X0.poly1:\n");
        mpq_homogenous_sparse_polynomial_print(wpoly3, 2);
        
        //wpoly4 = poly2(Xn = X0).polyinv
        printf("\npolyinv:\n");
        mpq_homogenous_sparse_polynomial_print(polyinv, 2);
        mpq_homogenous_sparse_polynomial_mul(wpoly5, &poly2, polyinv);
        mpq_homogenous_sparse_polynomial_div_monom(wpoly4, wpoly5, nind, indlimbsize);
        printf("\npoly2.polyinv:\n");
        mpq_homogenous_sparse_polynomial_print(wpoly4, 2);
        isinv = mpq_homogenous_sparse_polynomial_reduce(wpolyrem, wpolyquot, wpoly4, wpoly3, wpoly5);
        
        printf("\nwpolyrem + wpolyquot.wpoly3:\n");
        mpq_homogenous_sparse_polynomial_mul(wpoly4, wpolyquot, wpoly3);
        mpq_homogenous_sparse_polynomial_add(wpoly5, wpolyrem, wpoly4);
        mpq_homogenous_sparse_polynomial_print(wpoly5, 2);
        //Incrémenter polyinv
        mpq_homogenous_sparse_polynomial_mul_monom(wpolyrem, wpolyquot, nind, indlimbsize, &one);
        mpq_homogenous_sparse_polynomial_addincr(polyinv, wpolyrem);
        poly1.n = n;
        poly1.coeffsize += coeffsize2;
        printf("\npoly1.polyinv:\n");
        mpq_homogenous_sparse_polynomial_mul(wpolyrem, &poly1, polyinv);
        mpq_homogenous_sparse_polynomial_print(wpolyrem, 2);
        n++;
        mpn_lshift(nind, nind, indlimbsize, 1);
        poly2.n = n;
        poly2.indexes += indlimbsize * coeffsize2;
        poly2.coefficients += coeffsize2;
    }
    
    /*Déallocation*/
    mpq_clear(one);
    free(nind);
    
    return isinv;
}

char mpq_homogenous_sparse_polynomial_inversible2(mpq_homogenous_sparse_polynomial *poly, mpq_homogenous_sparse_polynomial *polyinv, mpq_homogenous_sparse_polynomial *wpoly) {
    /*Renvoie 1 si il existe polyinv tel que poly.polyinv = 1*/
    size_t coeffsize, coeffsize3;
    size_t degree;
    size_t nmax;
    size_t k, kmax, i, j;
    mp_size_t indlimbsize;
    mp_limb_t *ind3, *wind3, *wind;
    mpq_t coeff;
    
    /*Initialisation des variables*/
    coeffsize = poly->coeffsize;
    degree = poly->degree;
    nmax = poly->n;
    kmax = 1<<nmax;
    indlimbsize = poly->indlimbsize;
    mpq_init(coeff);
    //ind3 = calloc(indlimbsize, sizeof(mp_limb_t));
    wind = calloc(indlimbsize, sizeof(mp_limb_t));
    
    /*Inversion du coefficient constant*/
    if (!coeffsize || !mpn_zero_p(poly->indexes, indlimbsize)) {
        /*Coefficient constant nul*/
        return 0;
    }
    polyinv->coeffsize = 1;
    polyinv->degree = 1;
    mpn_zero(polyinv->indexes, indlimbsize);
    mpq_inv(*polyinv->coefficients, *poly->coefficients);
    
    /*Initialisation de wpoly = polyinv.poly*/
    mpq_homogenous_sparse_polynomial_mul_monom(wpoly, poly, polyinv->indexes, indlimbsize, polyinv->coefficients);
    k = 0;
    
    while (wpoly->coeffsize != 1 && k<kmax) {
        /*Ajouter ∑ a_I X^I à polyinv, où I parcoure l'ensemble des parties à k éléments de [1, nmax]*/
        k++;
        coeffsize3 = wpoly->coeffsize;
        ind3 = wpoly->indexes;
        for (i=1; i<coeffsize3; i++) {
            /*Parcourir les indices à la recherche d'un correspondant à une partie à k éléments*/
            ind3 += indlimbsize;
            if (mpn_popcount(ind3, indlimbsize) == k) {
                //mpn_copyi(ind3, *(&wpoly->indexes + i*indlimbsize), indlimbsize);
                printf("%lu %lu %lu\n", i, *ind3, mpn_popcount(ind3, indlimbsize));
                mpq_set_si(coeff, 0, 1);
                wind3 = wpoly->indexes;
                for (j=0; j<coeffsize3; j++) {
                    mpn_andn_n(wind, ind3, wind3, indlimbsize);
                    if (mpn_cmp(wind, ind3, indlimbsize) == 0) {
                        mpq_add(coeff, coeff, wpoly->coefficients[j]);
                    }
                    wind3 += indlimbsize;
                }
                if (mpq_sgn(coeff) != 0) {
                    /*coeff != 0*/
                    mpq_div(coeff, wpoly->coefficients[i], coeff);
                    mpq_neg(coeff, coeff);
                    //mpq_homogenous_sparse_polynomial_add_monom(polyinv, polyinv, ind3, indlimbsize, &coeff);
                    mpn_copyi(wind, ind3, indlimbsize);
                    mpq_homogenous_sparse_polynomial_add_mul_monom(wpoly, wpoly, poly, wind, indlimbsize, &coeff);
                }
            }
        }
    }
    
    /*Déallocation*/
    mpq_clear(coeff);
    //free(ind3);
    free(wind);
    
    return (wpoly->coeffsize == 1);
}
