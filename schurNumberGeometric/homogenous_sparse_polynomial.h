//
//  homogenous_sparse_polynomial.h
//  SchurNumber
//
//  Created by rubis on 16/03/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#ifndef homogenous_sparse_polynomial_h
#define homogenous_sparse_polynomial_h

#include <stdlib.h>
#include <gmp.h>

#if GMP_NUMB_BITS == 64
#define GMP_2EXP_NUMB_BITS 6
#elif GMP_NUMB_BITS == 128
#define GMP_2EXP_NUMB_BITS 7
#elif GMP_NUMB_BITS == 32
#define GMP_2EXP_NUMB_BITS 5
#endif

struct sparse_polynomial_struc {
    size_t coeffalloc;          // Nombre maximal de coefficients
    size_t coeffsize;           // Nombre de coefficients utilisés
    size_t degree;              // Degré homogène
    size_t n;                   // Nombre d'indéterminées, i.e. poly appartient à Q[X0,…,Xn]
    
    mp_size_t indlimbsize;      // Nombre de limbes utilisés par indice
    mp_limb_t *indexes;         // Tableau contenant les indices ordonnés lexicographiquement des coefficients ≠ 0
    
    //mpq_t *coefficients;        // Tableau des coefficients
    long *coefficients;         // Tableau des coefficients
};

typedef struct sparse_polynomial_struc sparse_polynomial;

sparse_polynomial * sparse_polynomial_alloc(size_t coeffalloc, mp_size_t indlimbsize, size_t degree);

char sparse_polynomial_realloc(sparse_polynomial *poly, size_t coeffalloc, mp_size_t indlimbsize);

void sparse_polynomial_copy(sparse_polynomial *polyr, sparse_polynomial *polys);

void sparse_polynomial_free(sparse_polynomial *poly);

void sparse_polynomial_add(sparse_polynomial *polyr, sparse_polynomial *polys1, sparse_polynomial *polys2);

void sparse_polynomial_add_monom(sparse_polynomial *polyr, sparse_polynomial *polys, mp_limb_t *monomindex, mp_size_t monomlimbsize, long monomcoeff);

void sparse_polynomial_addincr(sparse_polynomial *polyr, sparse_polynomial *polys);

//void sparse_polynomial_sub_monom(sparse_polynomial *polyr, sparse_polynomial *polys, mp_limb_t *monomindex, mp_size_t monomlimbsize, mpq_t *monomcoeff);

void sparse_polynomial_mul_monom(sparse_polynomial *polyr, sparse_polynomial *polys, mp_limb_t *monomindex, mp_size_t monomlimbsize, long monomcoeff);

void sparse_polynomial_add_mul_monom(sparse_polynomial *polyr, sparse_polynomial *polys1, sparse_polynomial *polys2, mp_limb_t *monomindex, mp_size_t monomlimbsize, long monomcoeff);

void sparse_polynomial_mul(sparse_polynomial *polyr, sparse_polynomial *polys1, sparse_polynomial *polys2);

//void sparse_polynomial_add_div_monom(sparse_polynomial *polyr, sparse_polynomial *polys1, sparse_polynomial *polys2, mp_limb_t *monomindex, mp_size_t monomlimbsize);

//void sparse_polynomial_sub_mul_monom(sparse_polynomial *polyr, sparse_polynomial *polys1, sparse_polynomial *polys2, mp_limb_t *monomindex, mp_size_t monomlimbsize, mpq_t *monomcoeff);

//char sparse_polynomial_reduce(sparse_polynomial *polyrem, sparse_polynomial *polyquot, sparse_polynomial *poly1, sparse_polynomial *poly2, sparse_polynomial *wpoly);

//char sparse_polynomial_inversible(sparse_polynomial *poly, sparse_polynomial *polyinv, sparse_polynomial *wpolyrem, sparse_polynomial *wpolyquot, sparse_polynomial *wpoly3, sparse_polynomial *wpoly4, sparse_polynomial *wpoly5);

void sparse_polynomial_print(sparse_polynomial *poly, unsigned int p);

#endif /* homogenous_sparse_polynomial_h */
