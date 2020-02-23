//
//  main.c
//  schurNumberGeometric
//
//  Created by rubis on 28/02/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include <stdio.h>
//#include "homogenous_sparse_polynomial.h"
//#include "schur_number_poly.c"
#include "schur_number_order_reduction.c"

void sparse_polynomial_print(sparse_polynomial *poly, unsigned int p) {
    size_t coeffsize, i;
    mp_size_t indlimbsize;
    mp_limb_t *ind;
    long *coeff;
    
    coeffsize = poly->coeffsize;
    indlimbsize = poly->indlimbsize;
    ind = poly->indexes;
    coeff = poly->coefficients;
    printf("n: %lu\n", poly->n);
    //printf("deg: %lu\n", poly->degree);
    for (i=0; i<coeffsize; i++) {
        printf("%lu %ld\n", *ind, *coeff);
        ind += indlimbsize;
        coeff++;
    }
}

int main(int argc, const char * argv[]) {
    //printf("Schur Number S(2) = %lu\n", schurNumber(2));
    printf("Schur Number S(2) = %lu\n", schurNumberWithOrderReduction(2));
    printf("Schur Number S(3) = %lu\n", schurNumberWithOrderReduction(3));
    printf("Schur Number S(4) = %lu\n", schurNumberWithOrderReduction(4));
    return 0;
}
