//
//  schur_number_poly.c
//  schurNumberGeometric
//
//  Created by rubis on 08/03/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include "homogenous_sparse_polynomial.h"

void mpq_homogenous_sparse_polynomial_print(mpq_homogenous_sparse_polynomial *poly, unsigned int p) {
    size_t coeffsize, i;
    mp_size_t indlimbsize;
    mp_limb_t *ind;
    mpq_t *coeff;
    
    coeffsize = poly->coeffsize;
    indlimbsize = poly->indlimbsize;
    ind = poly->indexes;
    coeff = poly->coefficients;
    printf("n: %lu\n", poly->n);
    //printf("deg: %lu\n", poly->degree);
    for (i=0; i<coeffsize; i++) {
        gmp_printf("%lu %Qd\n", *ind, coeff);
        ind += indlimbsize;
        coeff++;
    }
}

unsigned long schurNumber0(unsigned int p) {
    unsigned long n, ndiv2, i, nmax;
    unsigned int j;
    mpq_homogenous_sparse_polynomial *sfpartpoly;
    mpq_homogenous_sparse_polynomial *sfpartpolyinv;
    mpq_homogenous_sparse_polynomial *wpoly1;
    mpq_homogenous_sparse_polynomial *wpoly2;
    mpq_homogenous_sparse_polynomial *wpoly3;
    mpq_homogenous_sparse_polynomial *wpoly4;
    mpq_homogenous_sparse_polynomial *wpoly5;
    size_t partcoeffsize, sfcoeffsize, coeffalloc;
    mp_size_t indlimbsize;
    mp_limb_t *sfindex, *partindex, *index0;
    mp_limb_t *sfindex0, *sfindex1, *sfindex2;
    mpq_t *coeff0, *partcoeff, *sfcoeff;
    mpz_t pexp, sfpexp;
    
    indlimbsize = p;
    nmax = mp_bits_per_limb * p;
    if (p >= 64) {
        indlimbsize *= p;
        nmax *= p;
    }
    coeffalloc = 1024*p*p;
    /*Allocation des polynômes*/
    sfpartpoly = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    sfpartpolyinv = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    wpoly1 = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    wpoly2 = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    wpoly3 = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    wpoly4 = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    wpoly5 = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    
    /*Allocation des variables*/
    sfindex0 = calloc(indlimbsize, sizeof(mp_limb_t));
    sfindex1 = calloc(indlimbsize, sizeof(mp_limb_t));
    sfindex2 = calloc(indlimbsize, sizeof(mp_limb_t));
    mpz_init_set_si(pexp, 1);
    mpz_init_set_si(sfpexp, 1);
    
    index0 = sfpartpoly->indexes;
    mpn_zero(index0, indlimbsize);
    coeff0 = sfpartpoly->coefficients;
    mpq_set_si(*coeff0, 0, 1);
    sfpartpoly->coeffsize = 1;
    
    /*Remplir wpoly4 avec partpoly et wpoly5 avec sfpoly*/
    partindex = wpoly4->indexes;
    partcoeff = wpoly4->coefficients;
    coeff0 = partcoeff;
    mpn_zero(partindex, indlimbsize);
    mpq_set_si(*partcoeff, 0, 1);
    partindex += indlimbsize;
    partcoeff++;
    partcoeffsize = 1;
    
    sfpartpoly->coeffsize = 1;
    sfindex = wpoly5->indexes;
    sfcoeff = wpoly5->coefficients;
    sfcoeffsize = 0;
    
    for (n=1; n<=p; n++) {
        ndiv2 = (n+1)>>1;
        if (mpz_cmp_si(sfpexp, ndiv2*(ndiv2 + 1)) <= 0) {
            mpz_mul_si(sfpexp, sfpexp, p+1);
        }
        mpz_mul_si(pexp, pexp, p+1);
        mpz_sub(&((*coeff0)->_mp_num), &((*coeff0)->_mp_num), pexp);
        
        for (j=0; j<p; j++) {
            /*Ajout de (p+1)^n X_{n,j} à partpoly*/
            if (n==1 && !j) {
                *partindex = (mp_limb_t)1;
            } else {
                mpn_lshift(partindex, partindex - indlimbsize, indlimbsize, 1);
            }
            mpz_set(&(*partcoeff)->_mp_num, pexp);
            mpz_set_si(&(*partcoeff)->_mp_den, 1);
            partcoeff++;
            partcoeffsize++;
            
            /*Ajout de pexp X{n,j} ∑_{1≤i<ndiv2} X{i,j} X{n-i+1,j} à sfpoly*/
            mpn_zero(sfindex1, indlimbsize);
            *sfindex1 = (mp_limb_t)1<<j;
            mpn_zero(sfindex2, indlimbsize);
            sfindex2[n>>6] = (mp_limb_t)1<<(n%64 - j);
            for (i=0; i<ndiv2; i++) {
                mpn_ior_n(sfindex0, sfindex1, sfindex2, indlimbsize);
                mpn_ior_n(sfindex, partindex, sfindex0, indlimbsize);
                mpz_set(&(*sfcoeff)->_mp_num, sfpexp);
                mpz_set_ui(&(*sfcoeff)->_mp_den, 1);
                sfindex += indlimbsize;
                sfcoeff++;
                sfcoeffsize++;
                mpn_lshift(sfindex1, sfindex1, indlimbsize, p);
                mpn_rshift(sfindex2, sfindex2, indlimbsize, p);
            }
            partindex += indlimbsize;
        }
    }
    wpoly4->coeffsize = partcoeffsize;
    wpoly4->n = n*p;
    printf("Partpoly :\n");
    mpq_homogenous_sparse_polynomial_print(wpoly4, p);
    wpoly5->coeffsize = sfcoeffsize;
    wpoly5->n = n*p;
    printf("SFpoly :\n");
    mpq_homogenous_sparse_polynomial_print(wpoly5, p);
    mpq_homogenous_sparse_polynomial_add(sfpartpoly, wpoly4, wpoly5);
    printf("SFPartpoly :\n");
    mpq_homogenous_sparse_polynomial_print(sfpartpoly, p);
    
    while (!mpq_homogenous_sparse_polynomial_inversible(sfpartpoly, sfpartpolyinv, wpoly1, wpoly2, wpoly3, wpoly4, wpoly5) && n<4) {
        ndiv2 = (n+1)>>1;
        if (mpz_cmp_si(sfpexp, ndiv2*(ndiv2 + 1)) <= 0) {
            mpz_mul_si(sfpexp, sfpexp, p+1);
        }
        mpz_mul_si(pexp, pexp, p+1);
        mpz_sub(&((*coeff0)->_mp_num), &((*coeff0)->_mp_num), pexp);
        
        if (n >= nmax) {
            indlimbsize *= 2;
            nmax *= 2;
            mpq_homogenous_sparse_polynomial_realloc(sfpartpoly, coeffalloc, indlimbsize);
        }
        
        for (j=0; j<p; j++) {
            /*Ajout de (p+1)^n X_{n,j} à partpoly*/
            if (n==1 && !j) {
                *partindex = (mp_limb_t)1;
            } else {
                mpn_lshift(partindex, partindex - indlimbsize, indlimbsize, 1);
            }
            mpz_set(&(*partcoeff)->_mp_num, pexp);
            mpz_set_si(&(*partcoeff)->_mp_den, 1);
            partcoeff++;
            partcoeffsize++;
            
            /*Ajout de pexp X{n,j} ∑_{1≤i<ndiv2} X{i,j} X{n-i+1,j} à sfpoly*/
            mpn_zero(sfindex1, indlimbsize);
            *sfindex1 = (mp_limb_t)1<<j;
            mpn_zero(sfindex2, indlimbsize);
            sfindex2[n>>6] = (mp_limb_t)1<<(n%64 - j);
            for (i=0; i<ndiv2; i++) {
                mpn_ior_n(sfindex0, sfindex1, sfindex2, indlimbsize);
                mpn_ior_n(sfindex, partindex, sfindex0, indlimbsize);
                mpz_set(&(*sfcoeff)->_mp_num, sfpexp);
                mpz_set_ui(&(*sfcoeff)->_mp_den, 1);
                sfindex += indlimbsize;
                sfcoeff++;
                sfcoeffsize++;
                mpn_lshift(sfindex1, sfindex1, indlimbsize, p);
                mpn_rshift(sfindex2, sfindex2, indlimbsize, p);
            }
            partindex += indlimbsize;
        }

        n++;
    }
    
    /*Déallocation des polynômes*/
    mpq_homogenous_sparse_polynomial_free(sfpartpoly);
    mpq_homogenous_sparse_polynomial_free(sfpartpolyinv);
    mpq_homogenous_sparse_polynomial_free(wpoly1);
    mpq_homogenous_sparse_polynomial_free(wpoly2);
    mpq_homogenous_sparse_polynomial_free(wpoly3);
    mpq_homogenous_sparse_polynomial_free(wpoly4);
    mpq_homogenous_sparse_polynomial_free(wpoly5);
    
    /*Déallocation des autres variables*/
    free(sfindex0);
    free(sfindex1);
    free(sfindex2);
    mpz_clear(pexp);
    mpz_clear(sfpexp);
    
    return n;
}

unsigned long schurNumber(unsigned int p) {
    unsigned long n, ndiv2, nmod2, i, nmax;
    unsigned int j;
    mpq_homogenous_sparse_polynomial *sfpartpoly;
    mpq_homogenous_sparse_polynomial *sfpartpolyinv;
    mpq_homogenous_sparse_polynomial *wpoly1;
    mpq_homogenous_sparse_polynomial *wpoly2;
    mpq_homogenous_sparse_polynomial *wpoly3;
    mpq_homogenous_sparse_polynomial *wpoly4;
    mpq_homogenous_sparse_polynomial *wpoly5;
    size_t coeffalloc;
    mp_size_t indlimbsize;
    mp_limb_t *index, *partindex;
    mp_limb_t *sfindex0, *sfindex1, *sfindex2;
    mpq_t *coeff0, coeff;
    mpz_t pexp, sfpexp, sfcoeffmax;
    
    indlimbsize = p;
    nmax = mp_bits_per_limb * p;
    if (p >= 64) {
        indlimbsize *= p;
        nmax *= p;
    }
    coeffalloc = 1024*p*p;
    /*Allocation des polynômes*/
    sfpartpoly = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    sfpartpolyinv = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    wpoly1 = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    wpoly2 = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    wpoly3 = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    wpoly4 = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    wpoly5 = mpq_homogenous_sparse_polynomial_alloc(coeffalloc, indlimbsize, 3);
    
    /*Allocation des variables*/
    sfindex0 = calloc(indlimbsize, sizeof(mp_limb_t));
    sfindex1 = calloc(indlimbsize, sizeof(mp_limb_t));
    sfindex2 = calloc(indlimbsize, sizeof(mp_limb_t));
    partindex = calloc(indlimbsize, sizeof(mp_limb_t));
    mpq_init(coeff);
    mpz_init_set_si(pexp, p+1);
    mpz_init_set_si(sfpexp, 1);
    mpz_init_set_ui(sfcoeffmax, 0);
    
    index = sfpartpoly->indexes;
    mpn_zero(index, indlimbsize);
    coeff0 = sfpartpoly->coefficients;
    mpq_set_si(*coeff0, 0, 1);
    sfpartpoly->coeffsize = 1;
    
    /*Initialement, sfpartpoly = (p+1)(∑ X_{1,j} - 1)*/
    *partindex = (mp_limb_t)1;
    mpz_sub(&((*coeff0)->_mp_num), &((*coeff0)->_mp_num), pexp);
    for (j=0; j<p; j++) {
        /*Ajout de (p+1) X_{1,j} à partpoly*/
        if (j) {
            mpn_lshift(partindex, partindex, indlimbsize, 1);
        }
        mpz_set(&coeff->_mp_num, pexp);
        mpz_set_ui(&coeff->_mp_den, 1);
        mpq_homogenous_sparse_polynomial_add_monom(sfpartpoly, sfpartpoly, partindex, indlimbsize, &coeff);
    }
    
    for (n=2; n<=p; n++) {
        ndiv2 = n>>1;
        nmod2 = n%2;
        /*if (mpz_cmp_si(sfpexp, ndiv2*(ndiv2 + 1)) <= 0) {
            mpz_mul_si(sfpexp, sfpexp, p+1);
        }*/
        mpz_addmul_ui(sfcoeffmax, sfpexp, p*ndiv2*(ndiv2 + nmod2));
        if (mpz_cmp(sfpexp, sfcoeffmax) <= 0) {
            mpz_addmul_ui(sfcoeffmax, sfpexp, (p+1)*p*ndiv2*(ndiv2 + nmod2));
            mpz_mul_si(sfpexp, sfpexp, p+2);
        }
        mpz_mul_si(pexp, pexp, p+1);
        mpz_sub(&((*coeff0)->_mp_num), &((*coeff0)->_mp_num), pexp);
        
        for (j=0; j<p; j++) {
            /*Ajout de (p+1)^n X_{n,j} à partpoly*/
            mpn_lshift(partindex, partindex, indlimbsize, 1);
            mpz_set(&coeff->_mp_num, pexp);
            mpz_set_ui(&coeff->_mp_den, 1);
            mpq_homogenous_sparse_polynomial_add_monom(sfpartpoly, sfpartpoly, partindex, indlimbsize, &coeff);
            //printf("SFPartpoly %u:\n", j);
            //mpq_homogenous_sparse_polynomial_print(sfpartpoly, p);
            
            /*Ajout de pexp X{n,j} ∑_{1≤i<ndiv2} X{i,j} X{n-i+1,j} à sfpoly*/
            mpn_zero(sfindex1, indlimbsize);
            *sfindex1 = (mp_limb_t)1<<j;
            mpn_zero(sfindex2, indlimbsize);
            sfindex2[(p*(n-1) + j)>>6] = (mp_limb_t)1<<((p*(n-1) + j)%64);
            mpz_set(&coeff->_mp_num, sfpexp);
            for (i=0; i<ndiv2; i++) {
                mpn_ior_n(sfindex0, sfindex1, sfindex2, indlimbsize);
                mpn_ior_n(sfindex0, partindex, sfindex0, indlimbsize);
                mpq_homogenous_sparse_polynomial_add_monom(sfpartpoly, sfpartpoly, sfindex0, indlimbsize, &coeff);
                mpn_lshift(sfindex1, sfindex1, indlimbsize, p);
                mpn_rshift(sfindex2, sfindex2, indlimbsize, p);
                //printf("SFPartpoly %lu %u:\n", i, j);
                //mpq_homogenous_sparse_polynomial_print(sfpartpoly, p);
            }
        }
    }
    printf("SFPartpoly :\n");
    mpq_homogenous_sparse_polynomial_print(sfpartpoly, p);
    
    while (!mpq_homogenous_sparse_polynomial_inversible2(sfpartpoly, sfpartpolyinv, wpoly1) && n<8) {
    /*while (!mpq_homogenous_sparse_polynomial_inversible(sfpartpoly, sfpartpolyinv, wpoly1, wpoly2, wpoly3, wpoly4, wpoly5) && n<8) {*/
        ndiv2 = n>>1;
        nmod2 = n%2;
        mpz_addmul_ui(sfcoeffmax, sfpexp, p*ndiv2*(ndiv2 + nmod2));
        if (mpz_cmp(sfpexp, sfcoeffmax) <= 0) {
            mpz_addmul_ui(sfcoeffmax, sfpexp, (p+1)*p*ndiv2*(ndiv2 + nmod2));
            mpz_mul_si(sfpexp, sfpexp, p+2);
        }
        mpz_mul_si(pexp, pexp, p+1);
        coeff0 = sfpartpoly->coefficients;
        mpz_sub(&((*coeff0)->_mp_num), &((*coeff0)->_mp_num), pexp);
        
        if (n >= nmax) {
            indlimbsize *= 2;
            nmax *= 2;
            mpq_homogenous_sparse_polynomial_realloc(sfpartpoly, coeffalloc, indlimbsize);
            mpq_homogenous_sparse_polynomial_realloc(sfpartpoly, coeffalloc, indlimbsize);
            mpq_homogenous_sparse_polynomial_realloc(wpoly1, coeffalloc, indlimbsize);
            mpq_homogenous_sparse_polynomial_realloc(wpoly2, coeffalloc, indlimbsize);
            mpq_homogenous_sparse_polynomial_realloc(wpoly3, coeffalloc, indlimbsize);
            mpq_homogenous_sparse_polynomial_realloc(wpoly4, coeffalloc, indlimbsize);
            mpq_homogenous_sparse_polynomial_realloc(wpoly5, coeffalloc, indlimbsize);
            sfindex0 = realloc(sfindex0, indlimbsize * sizeof(mp_limb_t));
            sfindex1 = realloc(sfindex1, indlimbsize * sizeof(mp_limb_t));
            sfindex2 = realloc(sfindex2, indlimbsize * sizeof(mp_limb_t));
            partindex = realloc(partindex, indlimbsize * sizeof(mp_limb_t));
        }
        
        for (j=0; j<p; j++) {
            /*Ajout de (p+1)^n X_{n,j} à partpoly*/
            mpn_lshift(partindex, partindex, indlimbsize, 1);
            mpz_set(&coeff->_mp_num, pexp);
            mpz_set_ui(&coeff->_mp_den, 1);
            mpq_homogenous_sparse_polynomial_add_monom(sfpartpoly, sfpartpoly, partindex, indlimbsize, &coeff);
            //printf("SFPartpoly %u:\n", j);
            //mpq_homogenous_sparse_polynomial_print(sfpartpoly, p);
            
            /*Ajout de pexp X{n,j} ∑_{1≤i<ndiv2} X{i,j} X{n-i+1,j} à sfpoly*/
            mpn_zero(sfindex1, indlimbsize);
            *sfindex1 = (mp_limb_t)1<<j;
            mpn_zero(sfindex2, indlimbsize);
            sfindex2[(p*(n-2) + j)>>6] = (mp_limb_t)1<<((p*(n-2) + j)%64);
            mpz_set(&coeff->_mp_num, sfpexp);
            for (i=0; i<ndiv2; i++) {
                mpn_ior_n(sfindex0, sfindex1, sfindex2, indlimbsize);
                mpn_ior_n(sfindex0, partindex, sfindex0, indlimbsize);
                mpq_homogenous_sparse_polynomial_add_monom(sfpartpoly, sfpartpoly, sfindex0, indlimbsize, &coeff);
                mpn_lshift(sfindex1, sfindex1, indlimbsize, p);
                mpn_rshift(sfindex2, sfindex2, indlimbsize, p);
                //printf("SFPartpoly %lu %u:\n", i, j);
                //mpq_homogenous_sparse_polynomial_print(sfpartpoly, p);
            }
        }
        printf("SFPartpoly :\n");
        mpq_homogenous_sparse_polynomial_print(sfpartpoly, p);
        
        n++;
    }
    
    /*Déallocation des polynômes*/
    mpq_homogenous_sparse_polynomial_free(sfpartpoly);
    mpq_homogenous_sparse_polynomial_free(sfpartpolyinv);
    mpq_homogenous_sparse_polynomial_free(wpoly1);
    mpq_homogenous_sparse_polynomial_free(wpoly2);
    mpq_homogenous_sparse_polynomial_free(wpoly3);
    mpq_homogenous_sparse_polynomial_free(wpoly4);
    mpq_homogenous_sparse_polynomial_free(wpoly5);
    
    /*Déallocation des autres variables*/
    free(sfindex0);
    free(sfindex1);
    free(sfindex2);
    free(partindex);
    mpq_clear(coeff);
    mpz_clear(pexp);
    mpz_clear(sfpexp);
    mpz_clear(sfcoeffmax);
    
    return n;
}

