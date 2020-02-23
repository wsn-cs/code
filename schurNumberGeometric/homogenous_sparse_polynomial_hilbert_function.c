//
//  homogenous_sparse_polynomial_hilbert_function.c
//  schurNumberGeometric
//
//  Created by rubis on 19/03/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include "homogenous_sparse_polynomial.h"

mpz_t * mpq_homogenous_sparse_polynomial_Hilbert_function(mpq_homogenous_sparse_polynomial *poly) {
    size_t n,max;
    mpz_t * hfvalues;
    mp_size_t indlimbsize;
    mp_limb_t *indm;
    
    nmax = poly->n;
    indlimbsize = poly->indlimbsize;
    
    indm = calloc(indlimbsize, sizeof(mp_limb_t));
    
    for (n=0; n<nmax; n++) {
        <#statements#>
    }
}
