//
//  homogenous_sparse_polynomial_add.c
//  schurNumberGeometric
//
//  Created by rubis on 16/03/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include "homogenous_sparse_polynomial.h"

void sparse_polynomial_add(sparse_polynomial *polyr, sparse_polynomial *polys1, sparse_polynomial *polys2) {
    /*Effectue polyr = polys1 + polys2
     Suppose que polyr->indlimbsize = polys1->indlimbsize = polys2->indlimbsize.
     Les pointeurs polys1 et polys2 doivent être distincts de polyr.*/
    char indcmp;
    size_t n, nmod64;
    size_t i1, i2;
    mp_size_t indlimbsize;
    mp_limb_t *inds1, *inds2, *indr;
    mpq_t *coeffs1, *coeffs2, *coeffr;
    size_t coeffsize1, coeffsize2, coeffsizer;
    size_t coeffalloc;
    
    indlimbsize = polys1->indlimbsize;
    coeffsize1 = polys1->coeffsize;
    coeffsize2 = polys2->coeffsize;
    coeffsizer = 0;
    
    /*Vérifie si polyr peut contenir le résultat*/
    coeffalloc = polyr->coeffalloc;
    if (coeffsize1 + coeffsize2 > coeffalloc) {
        coeffalloc *= 2;
        sparse_polynomial_realloc(polyr, coeffalloc, indlimbsize);
    }
    
    /*Initialisation pour polys1*/
    inds1 = polys1->indexes;
    coeffs1 = polys1->coefficients;
    
    /*Initialisation pour polys2*/
    inds2 = polys2->indexes;
    coeffs2 = polys2->coefficients;
    
    /*Initialisation pour polyr*/
    indr = polyr->indexes;
    coeffr = polyr->coefficients;
    
    i1 = 0;
    i2 = 0;
    while (i1 < coeffsize1 && i2 < coeffsize2) {
        /*Comparer les indices inds1 et inds2*/
        indcmp = mpn_cmp(inds1, inds2, indlimbsize);
        if (indcmp == 0) {
            /*Les indices sont identiques: additionner.*/
            mpn_copyi(indr, inds1, indlimbsize);
            mpq_add(*coeffr, *coeffs1, *coeffs2);
            if (!mpq_sgn(*coeffr)) {
                /*coeffr == 0*/
                coeffsizer--;
                indr--;
                coeffr--;
            }
            i1++;
            i2++;
            inds1 += indlimbsize;
            inds2 += indlimbsize;
            coeffs1++;
            coeffs2++;
        } else if (indcmp > 0) {
            /*L'indice 1 est plus grand que l'indice 2.*/
            mpn_copyi(indr, inds2, indlimbsize);
            mpq_set(*coeffr, *coeffs2);
            i2++;
            inds2 += indlimbsize;
            coeffs2++;
        } else {
            /*L'indice 1 est plus petit que l'indice 2.*/
            mpn_copyi(indr, inds1, indlimbsize);
            mpq_set(*coeffr, *coeffs1);
            i1++;
            inds1 += indlimbsize;
            coeffs1++;
        }
        
        coeffsizer++;
        indr += indlimbsize;
        coeffr++;
    }
    
    if (i1 == coeffsize1) {
        /*Tous les coefficients de polys1 ont été parcourus*/
        while (i2 < coeffsize2) {
            /*Parcourir les coefficients restant de polys2*/
            mpn_copyi(indr, inds2, indlimbsize);
            mpq_set(*coeffr, *coeffs2);
            i2++;
            inds2 += indlimbsize;
            coeffs2++;
            
            coeffsizer++;
            indr += indlimbsize;
            coeffr++;
        }
    } else {
        /*Il reste des coefficients de polys1 à parcourir*/
        while (i1 < coeffsize1) {
            /*Parcourir les coefficients restant de polys1*/
            mpn_copyi(indr, inds1, indlimbsize);
            mpq_set(*coeffr, *coeffs1);
            i1++;
            inds1 += indlimbsize;
            coeffs1++;
            
            coeffsizer++;
            indr += indlimbsize;
            coeffr++;
        }
    }
    polyr->coeffsize = coeffsizer;
    polyr->degree = polys1->degree;
    
    /*Nombre d'indéterminées*/
    n = polys1->n;
    if (n < polys2->n) {
        n = polys2->n;
    }
    nmod64 = n%64;
    indr--;
    while (n > 0 && !(*indr & (mp_limb_t)1<<(nmod64-1))) {
        n--;
        //ndiv64 = n>>6;
        nmod64 = n%64;
        if (nmod64 == 0) {
            indr--;
        }
    }
    polyr->n = n;
}

void sparse_polynomial_add_monom(sparse_polynomial *polyr, sparse_polynomial *polys, mp_limb_t *monomindex, mp_size_t monomlimbsize, long monomcoeff) {
    /*Ajoute monomcoeff.X^(monomind) à polys et met le résultat dans polyr.
     Suppose que polyr->indlimbsize = polys->indlimbsize = monomlimbsize*/
    char b;
    char monomhasbeenplaced;
    size_t n, nmod64;
    size_t i, coeffsize, coeffalloc;
    mp_limb_t *inds, *indr, *indw1, *indw2;
    long *coeffs, *coeffr, coeffw1, coeffw2;
    
    /*Vérifie si polyr peut contenir le résultat*/
    coeffsize = polys->coeffsize;
    coeffalloc = polyr->coeffalloc;
    if (coeffsize >= coeffalloc) {
        coeffalloc *= 2;
        sparse_polynomial_realloc(polyr, coeffalloc, monomlimbsize);
    }
    
    /*Initialisation pour polys*/
    inds = polys->indexes;
    coeffs = polys->coefficients;
    
    /*Initialisation pour polyr*/
    indr = polyr->indexes;
    coeffr = polyr->coefficients;
    polyr->coeffsize = coeffsize;
    polyr->degree = polys->degree;
    
    i = 0;
    monomhasbeenplaced = 0;
    while (i < coeffsize) {
        b = mpn_cmp(inds, monomindex, monomlimbsize);
        if (b > 0) {
            /*inds > monomindex*/
            monomhasbeenplaced = 1;
            polyr->coeffsize ++;    //Le monôme crée un nouveau coefficient.
            /*Conserver les coefficients inds et coeffs de polys*/
            indw1 = calloc(monomlimbsize, sizeof(mp_limb_t));
            indw2 = calloc(monomlimbsize, sizeof(mp_limb_t));
            mpn_copyi(indw1, inds, monomlimbsize);
            coeffw1 = *coeffs;
            /*Placer le monôme*/
            mpn_copyi(indr, monomindex, monomlimbsize);
            *coeffr = monomcoeff;
            indr += monomlimbsize;
            coeffr++;
            //i--;
            /*Parcourir le reste de polys*/
            while (i < coeffsize) {
                /*Conserver les coefficients inds et coeffs*/
                inds += monomlimbsize;
                coeffs++;
                mpn_copyi(indw2, inds, monomlimbsize);
                coeffw2 = *coeffs;
                /*Remplacer le coefficient*/
                mpn_copyi(indr, indw1, monomlimbsize);
                *coeffr = coeffw1;
                /*Recopier*/
                mpn_copyi(indw1, indw2, monomlimbsize);
                coeffw1 = coeffw2;
                i++;
                indr += monomlimbsize;
                coeffr++;
            }
            free(indw1);
            free(indw2);
            break;
        } else if (b == 0) {
            /*inds == monomindex*/
            monomhasbeenplaced = 1;
            mpn_copyi(indr, inds, monomlimbsize);
            *coeffr = *coeffs + monomcoeff;
            i++;
            inds += monomlimbsize;
            coeffs++;
            if (!(*coeffr)) {
                /*coeffr == 0*/
                polyr->coeffsize --;
            } else {
                /*coeffr != 0*/
                indr += monomlimbsize;
                coeffr++;
            }
            break;
        }
        /*inds < monomindex*/
        mpn_copyi(indr, inds, monomlimbsize);
        *coeffr = *coeffs;
        i++;
        inds += monomlimbsize;
        coeffs++;
        indr += monomlimbsize;
        coeffr++;
    }
    
    while (i < coeffsize) {
        /*Si il reste des coefficients à parcourir*/
        mpn_copyi(indr, inds, monomlimbsize);
        *coeffr = *coeffs;
        i++;
        inds += monomlimbsize;
        coeffs++;
        indr += monomlimbsize;
        coeffr++;
    }
    
    if (!monomhasbeenplaced) {
        polyr->coeffsize ++;    //Le monôme crée un nouveau coefficient.
        //printf("%lu %lu\n", *monomindex, monomindex[1]);
        mpn_copyi(indr, monomindex, monomlimbsize);
        //printf("%p %lu %p %lu\n", indr, *indr, indr + indlimbincr - sizeof(mp_limb_t), *(indr + indlimbincr - sizeof(mp_limb_t)));
        *coeffr = monomcoeff;
        //printf("%p\n", coeffr);
        indr += monomlimbsize;
        //printf("%p\n", indr);
    }

    /*Nombre d'indéterminées*/
    n = 64*monomlimbsize;
    nmod64 = n%64;
    //printf("%p %lu %p %lu\n", indr - 2*sizeof(mp_limb_t), *(indr - 2*sizeof(mp_limb_t)), indr - sizeof(mp_limb_t), *(indr - sizeof(mp_limb_t)));
    indr--;
    //printf("%p %lu %p %lu\n", indr - sizeof(mp_limb_t), *(indr - sizeof(mp_limb_t)), indr, *indr);
    while (n > 0 && !(*indr & (mp_limb_t)1<<(nmod64-1))) {
        n--;
        //ndiv64 = n>>6;
        nmod64 = n%64;
        if (nmod64 == 1) {
            indr--;
            //printf("%lu %p %lu\n", n, indr, *indr);
        }
    }
    polyr->n = n;
}

void sparse_polynomial_addincr(sparse_polynomial *polyr, sparse_polynomial *polys) {
    /*Effectue polyr += polys.
     Suppose que polyr->indlimbsize = polys1->indlimbsize = polys2->indlimbsize.*/
    size_t coeffsize, i;
    mp_size_t indlimbsize;
    mp_limb_t *inds;
    long *coeffs;
    
    coeffsize = polys->coeffsize;
    indlimbsize = polys->indlimbsize;
    inds = polys->indexes;
    coeffs = polys->coefficients;
    
    for (i=0; i<coeffsize; i++) {
        sparse_polynomial_add_monom(polyr, polyr, inds, indlimbsize, *coeffs);
        inds += indlimbsize;
        coeffs++;
    }
}
