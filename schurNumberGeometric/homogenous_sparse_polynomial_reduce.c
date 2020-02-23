//
//  homogenous_sparse_polynomial_reduce.c
//  schurNumberGeometric
//
//  Created by rubis on 16/03/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include "homogenous_sparse_polynomial.h"
#include <stdio.h>

char mpq_homogenous_sparse_polynomial_reduce(mpq_homogenous_sparse_polynomial *polyrem, mpq_homogenous_sparse_polynomial *polyquot, mpq_homogenous_sparse_polynomial *poly1, mpq_homogenous_sparse_polynomial *poly2, mpq_homogenous_sparse_polynomial *wpoly) {
    /*Calcule les polynômes de degrés minimaux tels que poly1 = polyquot*poly2 + polyrem.
     Renvoie 1 si polyrem = 0, et 0 sinon.*/
    
    char poly2cst;
    size_t wcoeffsz, coeffsz2;
    size_t deg1, deg2, degq, degm;
    mp_size_t indlimbsize;
    mp_limb_t *wleadingind, *leadingind2, *monomind;
    mp_size_t setbitind2;
    mpq_t *wleadingcoeff, *leadingcoeff2, monomcoeff;
    
    /*Aucun test de compatibilité dimensionnelle*/
    wcoeffsz = poly1->coeffsize;
    coeffsz2 = poly2->coeffsize;
    deg1 = poly1->degree;
    deg2 = poly2->degree;
    indlimbsize = poly1->indlimbsize;
    
    /*Définition initiale à 0 de polyrem et de polyquot*/
    polyrem->coeffsize = 0;
    polyrem->degree = deg1;
    polyquot->coeffsize = 0;
    /*degq = deg1 - deg2;
    if (deg1 < deg2) {
        mpq_homogenous_sparse_polynomial_copy(polyrem, poly1);
        return !(polyrem->coeffsize);
    }
    polyquot->degree = degq;*/
    
    /*Définition de wpoly*/
    mpq_homogenous_sparse_polynomial_copy(wpoly, poly1);    //Polynôme restant à réduire
    //printf("wpoly:\n");
    //mpq_homogenous_sparse_polynomial_print(wpoly, 2);
    //printf("\npoly1:\n");
    //mpq_homogenous_sparse_polynomial_print(poly1, 2);
    
    /*Terme dominant de poly2*/
    leadingind2 = &poly2->indexes[indlimbsize * (coeffsz2 - 1)];
    leadingcoeff2 = &poly2->coefficients[coeffsz2 - 1];
    //setbitind2 = mpn_popcount(leadingind2, indlimbsize);    //Nombre de bits égaux à 1 dans leadingind2
    //degm = degq + setbitind2;   //Degré maximal d'un monôme à multiplier par poly2
    poly2cst = mpn_zero_p(leadingind2, indlimbsize);
    
    /*Allocation du monôme*/
    monomind = calloc(indlimbsize, sizeof(mp_limb_t));
    mpq_init(monomcoeff);
    
    //printf("\npoly1:\n");
    //mpq_homogenous_sparse_polynomial_print(poly1, 2);
    //printf("\npoly2:\n");
    //mpq_homogenous_sparse_polynomial_print(poly2, 2);
    
    while (wcoeffsz > 0) {
        /*Tant qu'il reste des coefficients potentiellement réductibles*/
        
        /*Terme dominant de wpoly*/
        wleadingind = &wpoly->indexes[indlimbsize * (wcoeffsz - 1)];
        wleadingcoeff = &wpoly->coefficients[wcoeffsz - 1];
        
        mpn_and_n(monomind, wleadingind, leadingind2, indlimbsize);
        while (mpn_cmp(monomind, leadingind2, indlimbsize) == 0 || poly2cst) {
        //while ((mpn_cmp(monomind, leadingind2, indlimbsize) == 0 || poly2cst) && mpn_popcount(wleadingind, indlimbsize) <= degm) {
            /*Il existe un monôme tel que lc1 = lc2 * monôme*/
            mpn_xor_n(monomind, wleadingind, leadingind2, indlimbsize);
            mpq_div(monomcoeff, *wleadingcoeff, *leadingcoeff2);
            mpq_homogenous_sparse_polynomial_add_monom(polyquot, polyquot, monomind, indlimbsize, &monomcoeff);
            mpq_neg(monomcoeff, monomcoeff);
            mpq_homogenous_sparse_polynomial_add_mul_monom(wpoly, wpoly, poly2, monomind, indlimbsize, &monomcoeff);
            /*Nouveau terme dominant de wpoly*/
            wcoeffsz = wpoly->coeffsize;
            if (!wcoeffsz) {
                break;
            }
            wleadingind = &wpoly->indexes[indlimbsize * (wcoeffsz - 1)];
            wleadingcoeff = &wpoly->coefficients[wcoeffsz - 1];
            mpn_and_n(monomind, wleadingind, leadingind2, indlimbsize);
        }
        
        if (!wcoeffsz) {
            break;
        }
        mpn_copyi(monomind, wleadingind, indlimbsize);
        mpq_set(monomcoeff, *wleadingcoeff);
        mpq_homogenous_sparse_polynomial_add_monom(polyrem, polyrem, monomind, indlimbsize, &monomcoeff);
        //mpq_neg(monomcoeff, monomcoeff);
        //mpq_neg(monomcoeff, *wleadingcoeff);
        //mpq_homogenous_sparse_polynomial_add_monom(wpoly, wpoly, monomind, indlimbsize, &monomcoeff);
        wpoly->coeffsize--;
        wcoeffsz = wpoly->coeffsize;
    }
    
    /*Déallocation du monôme*/
    free(monomind);
    mpq_clear(monomcoeff);
    
    printf("\npolyrem:\n");
    mpq_homogenous_sparse_polynomial_print(polyrem, 2);
    printf("\npolyquot:\n");
    mpq_homogenous_sparse_polynomial_print(polyquot, 2);
    
    return !(polyrem->coeffsize);
}
