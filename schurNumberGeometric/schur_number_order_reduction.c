//
//  schur_number_order_reduction.c
//  schurNumberGeometric
//
//  Created by rubis on 30/03/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include "homogenous_sparse_polynomial.h"
#include <string.h>

unsigned long schurNumberWithOrderReduction(unsigned long p) {
    unsigned long i, imodbpl, idivbpl, imax;
    unsigned long n, nmax;
    unsigned long sfnp, sfn, sfp;
    unsigned long sfnpmodbpl, sfnpdivbpl, sfndiv2, sfnmodp;
    mp_size_t indlimbsize;
    mp_limb_t *indexes;
    long *coefficients;
    sparse_polynomial onepoly;
    sparse_polynomial partpoly;
    sparse_polynomial wpoly;
    
    /*Initialisation des variables*/
    n = 1;
    nmax = p * mp_bits_per_limb;
    indlimbsize = p;
    
    /*Initialisation de partpoly*/
    partpoly.coeffalloc = p+1;
    partpoly.coeffsize = p+1;
    partpoly.degree = 1;
    partpoly.n = p;
    partpoly.indlimbsize = indlimbsize;
    indexes = calloc(p+1, indlimbsize * sizeof(mp_limb_t));
    partpoly.indexes = indexes;
    indexes += indlimbsize;
    coefficients = calloc(p+1, sizeof(long));
    partpoly.coefficients = coefficients;
    *coefficients = (long)-1;
    coefficients++;
    for (i=0; i<p; i++) {
        *indexes = (mp_limb_t)1<<i;
        *coefficients = (long)1;
        indexes += indlimbsize;
        coefficients++;
    }
    
    /*Initialisation de wpoly*/
    wpoly.coeffalloc = p+1;
    wpoly.indlimbsize = indlimbsize;
    indexes = calloc(p+1, indlimbsize * sizeof(mp_limb_t));
    wpoly.indexes = indexes;
    coefficients = calloc(p+1, sizeof(long));
    wpoly.coefficients = coefficients;
    
    /*Initialisation de onepoly*/
    onepoly.coeffalloc = 1<<(p+1);
    onepoly.coeffsize = 1;
    onepoly.degree = 0;
    onepoly.n = 0;
    onepoly.indlimbsize = indlimbsize;
    indexes = calloc(1<<(p+1), indlimbsize * sizeof(mp_limb_t));
    onepoly.indexes = indexes;
    coefficients = calloc(1<<(p+1), sizeof(long));
    onepoly.coefficients = coefficients;
    *coefficients = (long)1;
    
    while (onepoly.coeffsize) {
        /*Tant que onepoly ≠ 0*/
        //sfn = n - p;
        sfnp = p*(n - p);
        if (sfnp < p) {
            sfnp = p;
        }
        while (sfnp < p*n && onepoly.coeffsize) {
            /*sfnp = sfp + p*sfn*/
            sfn = sfnp / p;
            sfp = sfnp % p;
            sfnpmodbpl = sfnp % mp_bits_per_limb;
            sfnpdivbpl = sfnp / mp_bits_per_limb;
            
            //printf("Possibilité %lu %lu %lu\n", sfn, (mp_limb_t)1<<sfnmodbpl, sfnmodbpl);
            if (indexes[sfnpdivbpl] & ((mp_limb_t)1<<sfnpmodbpl)) {
                /*Un monôme X_{sfn+1} X_{i+1} X_{sfn -i} pourrait réduire onepoly*/
                i = sfp;
                imax = i + (sfn>>1)*p;
                while (i <= imax) {
                    imodbpl = i % mp_bits_per_limb;
                    idivbpl = i / mp_bits_per_limb;
                    i += p;
                    sfnpmodbpl = (sfnp - i + sfp) % mp_bits_per_limb;
                    sfnpdivbpl = (sfnp - i + sfp) / mp_bits_per_limb;
                    //printf("Monôme %lu %lu %lu %lu %lu\n", sfnp, sfn, sfp, i-p, ((mp_limb_t)1<<imodbpl) | ((mp_limb_t)1 << sfnpmodbpl) | ((mp_limb_t)1<<sfnp));
                    //i += p;
                    if (indexes[idivbpl] & ((mp_limb_t)1<<imodbpl) && indexes[sfnpdivbpl] & ((mp_limb_t)1 << sfnpmodbpl)) {
                        /*Le monôme X_{sfn+1} X_{i+1} X_{sfn -i} réduit effectivement onepoly.
                         Il faut supprimer le dernier coefficient de onepoly.*/
                        //printf("Réduction par %lu:\n", ((mp_limb_t)1<<imodbpl) | ((mp_limb_t)1 << sfnpmodbpl) | ((mp_limb_t)1<<sfnp));
                        onepoly.coeffsize--;
                        mpn_copyi(indexes, indexes + indlimbsize, onepoly.coeffsize * indlimbsize);
                        memmove(coefficients, coefficients + 1, onepoly.coeffsize * sizeof(long));
                        //mpn_copyi(coefficients, coefficients + sizeof(long), onepoly.coeffsize);
                        //sparse_polynomial_print(&onepoly, p);
                        //if (!onepoly.coeffsize) {
                        //    break;
                        //}
                        sfnp--;
                        //i = sfndiv2;
                        break;
                    }
                }
            }
            sfnp++;
        }
        
        /*Tester si partpoly réduit onepoly.*/
        sparse_polynomial_mul_monom(&wpoly, &partpoly, indexes, indlimbsize, *coefficients);
        if (!mpn_cmp(indexes, wpoly.indexes, indlimbsize)) {
            /*Si le monôme le plus bas de wpoly est identique à celui de onepoly, effectuer onepoly += wpoly.*/
            //printf("Réduction %lu:\n", n);
            //sparse_polynomial_print(&wpoly, p);
            //sparse_polynomial_print(&onepoly, p);
            sparse_polynomial_addincr(&onepoly, &wpoly);
            //sparse_polynomial_print(&onepoly, p);
            indexes = onepoly.indexes;
            coefficients = onepoly.coefficients;
        } else {
            /*Incrémenter n*/
            n++;
            if (n >= nmax) {
                printf("Reallocation failed.");
                break;
            }
            /*Modifier partpoly*/
            partpoly.n += p;
            indexes = partpoly.indexes + indlimbsize;
            for (i=0; i<p; i++) {
                mpn_lshift(indexes, indexes, indlimbsize, p);
                indexes += indlimbsize;
            }
            //sparse_polynomial_print(&partpoly, p);
            indexes = onepoly.indexes;
        }
    }
    
    /*Libération de mémoire*/
    free(partpoly.indexes);
    free(partpoly.coefficients);
    free(wpoly.indexes);
    free(wpoly.coefficients);
    free(onepoly.indexes);
    free(onepoly.coefficients);
    
    return n-2;
}

unsigned long weakSchurNumberWithOrderReduction(unsigned long p) {
    unsigned long i, imodbpl, idivbpl;
    unsigned long n, nmax;
    unsigned long sfn, sfnmodbpl, sfndivbpl, sfndiv2;
    mp_size_t indlimbsize;
    mp_limb_t *indexes;
    long *coefficients;
    sparse_polynomial onepoly;
    sparse_polynomial partpoly;
    sparse_polynomial wpoly;
    
    /*Initialisation des variables*/
    n = 1;
    nmax = p * mp_bits_per_limb;
    indlimbsize = p;
    
    /*Initialisation de partpoly*/
    partpoly.coeffalloc = p+1;
    partpoly.coeffsize = p+1;
    partpoly.degree = 1;
    partpoly.n = p;
    partpoly.indlimbsize = indlimbsize;
    indexes = calloc(p+1, indlimbsize * sizeof(mp_limb_t));
    partpoly.indexes = indexes;
    indexes += indlimbsize;
    coefficients = calloc(p+1, sizeof(long));
    partpoly.coefficients = coefficients;
    *coefficients = (long)-1;
    coefficients++;
    for (i=0; i<p; i++) {
        *indexes = (mp_limb_t)1<<i;
        *coefficients = (long)1;
        indexes += indlimbsize;
        coefficients++;
    }
    
    /*Initialisation de wpoly*/
    wpoly.coeffalloc = p+1;
    wpoly.indlimbsize = indlimbsize;
    indexes = calloc(p+1, indlimbsize * sizeof(mp_limb_t));
    wpoly.indexes = indexes;
    coefficients = calloc(p+1, sizeof(long));
    wpoly.coefficients = coefficients;
    
    /*Initialisation de onepoly*/
    onepoly.coeffalloc = 1<<(p+1);
    onepoly.coeffsize = 1;
    onepoly.degree = 0;
    onepoly.n = 0;
    onepoly.indlimbsize = indlimbsize;
    indexes = calloc(1<<(p+1), indlimbsize * sizeof(mp_limb_t));
    onepoly.indexes = indexes;
    coefficients = calloc(1<<(p+1), sizeof(long));
    onepoly.coefficients = coefficients;
    *coefficients = (long)1;
    
    while (onepoly.coeffsize) {
        /*Tant que onepoly ≠ 0*/
        //sfn = n - p;
        sfn = 2;
        sfnmodbpl = sfn % mp_bits_per_limb;
        sfndivbpl = sfn / mp_bits_per_limb ;
        sfndiv2 = sfn>>1;
        while (sfn < n && onepoly.coeffsize) {
            if (indexes[sfndivbpl] & ((mp_limb_t)1<<sfnmodbpl)) {
                /*Un monôme X_{sfn+1} X_{i+1} X_{sfn -i} pourrait réduire onepoly*/
                i = 0;
                imodbpl = 0;
                idivbpl = 0;
                sfnmodbpl = (sfn - i) % mp_bits_per_limb;
                sfndivbpl = (sfn - i) / mp_bits_per_limb;
                while (i <= sfndiv2) {
                    if (indexes[idivbpl] & ((mp_limb_t)1<<imodbpl) && indexes[sfndivbpl] & ((mp_limb_t)1 << sfnmodbpl)) {
                        /*Le monôme X_{sfn+1} X_{i+1} X_{sfn -i} réduit effectivement onepoly.
                         Il faut supprimer le dernier coefficient de onepoly.*/
                        printf("Réduction par %lu:\n", ((mp_limb_t)1<<imodbpl) | ((mp_limb_t)1 << sfnmodbpl) | (1<<sfn));
                        onepoly.coeffsize--;
                        mpn_copyi(indexes, indexes + indlimbsize, onepoly.coeffsize);
                        memmove(coefficients, coefficients + 1, onepoly.coeffsize * sizeof(long));
                        //mpn_copyi(coefficients, coefficients + sizeof(long), onepoly.coeffsize);
                        sparse_polynomial_print(&onepoly, p);
                        if (!onepoly.coeffsize) {
                            break;
                        }
                    }
                    i++;
                    imodbpl = i % mp_bits_per_limb;
                    idivbpl = i / mp_bits_per_limb;
                    sfnmodbpl = (sfn - i) % mp_bits_per_limb;
                    sfndivbpl = (sfn - i) / mp_bits_per_limb;
                }
            }
            sfn++;
            sfnmodbpl = sfn % mp_bits_per_limb;
            sfndivbpl = sfn / mp_bits_per_limb ;
            sfndiv2 = sfn>>1;
        }
        
        /*Tester si partpoly réduit onepoly.*/
        sparse_polynomial_mul_monom(&wpoly, &partpoly, indexes, indlimbsize, *coefficients);
        if (!mpn_cmp(indexes, wpoly.indexes, indlimbsize)) {
            /*Si le monôme le plus bas de wpoly est identique à celui de onepoly, effectuer onepoly += wpoly.*/
            //printf("Réduction %lu:\n", n);
            //sparse_polynomial_print(&wpoly, p);
            //sparse_polynomial_print(&onepoly, p);
            sparse_polynomial_addincr(&onepoly, &wpoly);
            //sparse_polynomial_print(&onepoly, p);
        } else {
            /*Incrémenter n*/
            n++;
            if (n >= nmax) {
                printf("Reallocation failed.");
                break;
            }
            /*Modifier partpoly*/
            partpoly.n += p;
            indexes = partpoly.indexes + indlimbsize;
            for (i=0; i<p; i++) {
                mpn_lshift(indexes, indexes, indlimbsize, p);
                indexes += indlimbsize;
            }
            indexes = onepoly.indexes;
        }
    }
    
    /*Libération de mémoire*/
    free(partpoly.indexes);
    free(partpoly.coefficients);
    free(wpoly.indexes);
    free(wpoly.coefficients);
    free(onepoly.indexes);
    free(onepoly.coefficients);
    
    return n-1;
}

unsigned long schurNumberWithOrderReduction0(unsigned int p) {
    char notreduced;
    unsigned long n, ndiv2, nmax;
    size_t I, sfIbegin, partIbegin;
    //size_t Iminimal;
    size_t *sfszarray;
    unsigned int i;
    size_t indlimbsize;
    mp_limb_t *ind, *polyindex, *oneindex;
    mp_limb_t *partindex;
    mp_limb_t *sfindex0, *sfindex1, *sfindex2;
    sparse_polynomial *onereduce;
    //long coeff;
    long *onecoeff;
    size_t sfsize, sfalloc;
    sparse_polynomial **sfpolys;  //Tableau des polynômes traduisant le caractère sans-somme
    size_t partsize, partalloc;
    sparse_polynomial **partpolys;//Tableau des polynômes traduisant le caractère partitionnant
    sparse_polynomial *poly, *sfpoly, *wpoly;
    
    /*Allocation*/
    nmax = p * mp_bits_per_limb;
    sfalloc = (nmax>>1) * (nmax>>1);
    sfszarray = calloc(p, sizeof(size_t));
    sfpolys = calloc(sfalloc, sizeof(sparse_polynomial *));
    
    partalloc = p * nmax;
    partpolys = calloc(partalloc, sizeof(sparse_polynomial *));
    
    indlimbsize = p;
    ind = calloc(p, sizeof(mp_limb_t));
    partindex = calloc(p, sizeof(mp_limb_t));
    sfindex0 = calloc(p, sizeof(mp_limb_t));
    sfindex1 = calloc(p, sizeof(mp_limb_t));
    sfindex2 = calloc(p, sizeof(mp_limb_t));
    wpoly = sparse_polynomial_alloc(1<<p, indlimbsize, 1);
    //mpq_init(coeff);
    
    /*Assignation de 1 à onereduce*/
    onereduce = sparse_polynomial_alloc(1<<p, indlimbsize, 0);
    printf("%lu %lu\n", *ind, ind[1]);
    sparse_polynomial_add_monom(onereduce, onereduce, ind, indlimbsize, (long)1);
    sparse_polynomial_print(onereduce, p);
    oneindex = onereduce->indexes;
    onecoeff = onereduce->coefficients;
    
    /*Remplissage initiale de sfpolys*/
    n = 2;
    for (i=0; i<p; i++) {
        sfpoly = NULL;
        sfpoly = sparse_polynomial_alloc(1, indlimbsize, 0);
        *ind = (mp_limb_t)((mp_limb_t)1<<i) | ((mp_limb_t)1<<(p+i));
        printf("%lu %lu\n", *ind, ind[1]);
        sparse_polynomial_add_monom(sfpoly, sfpoly, ind, indlimbsize, (long)1);
        sparse_polynomial_print(sfpoly, p);
        sfpolys[i] = sfpoly;
    }
    sfsize = p;
    sfIbegin = 0;
    sfszarray[2 % p] = p;
    //Iminimal = p;
    
    /*Remplissage initiale de partpolys*/
    *partpolys = sparse_polynomial_alloc(p+1, indlimbsize, 0);
    partpolys[1] = sparse_polynomial_alloc(p+1, indlimbsize, 0);
    *partindex = (mp_limb_t)0;
    sparse_polynomial_add_monom(*partpolys, *partpolys, partindex, indlimbsize, -1);
    sparse_polynomial_print(*partpolys, p);
    sparse_polynomial_add_monom(partpolys[1], partpolys[1], partindex, indlimbsize, -1);
    sparse_polynomial_print(partpolys[1], p);
    *partindex = (mp_limb_t)1;
    for (i=0; i<p; i++) {
        sparse_polynomial_add_monom(*partpolys, *partpolys, partindex, indlimbsize, (long)1);
        mpn_lshift(partindex, partindex, indlimbsize, 1);
    }
    for (i=0; i<p; i++) {
        sparse_polynomial_add_monom(partpolys[1], partpolys[1], partindex, indlimbsize, (long)1);
        mpn_lshift(partindex, partindex, indlimbsize, 1);
    }
    partsize = 2;
    partIbegin = 0;
    
    while (onereduce->coeffsize) {
        /*Tant que onereduce ≠ 0*/
        notreduced = 1;
        
        /*Réduire onereduce par sfpolys*/
        I = sfIbegin;
        while (I < sfsize && notreduced) {
            /*Tester si sfpolys[I] réduit onereduce*/
            poly = sfpolys[I];
            polyindex = poly->indexes;
            mpn_and_n(ind, polyindex, oneindex, indlimbsize);
            if (!mpn_cmp(ind, polyindex, indlimbsize)) {
                /*Retrancher à onereduce*/
                //mpq_neg(coeff, *onecoeff);
                mpn_copyd(ind, oneindex, indlimbsize);
                printf("Réduction par le monôme d'indice %lu au rang %lu.\n", I, n);
                //if (Iminimal > I) {
                //    Iminimal = I;
                //}
                sparse_polynomial_print(onereduce, p);
                sparse_polynomial_add_monom(onereduce, onereduce, ind, indlimbsize, - *onecoeff);
                sparse_polynomial_print(onereduce, p);
                notreduced = 0;
                oneindex = onereduce->indexes;
                onecoeff = onereduce->coefficients;
            }
            I++;
        }
        
        /*Réduire onereduce par partpolys*/
        I = partIbegin;
        while (I < partsize && notreduced) {
            /*Tester si partpolys[I] réduit onereduce*/
            poly = partpolys[I];
            sparse_polynomial_mul_monom(wpoly, poly, oneindex, indlimbsize, *onecoeff);
            if (!mpn_cmp(wpoly->indexes, oneindex, indlimbsize)) {
                printf("Réduction selon une partition:\n");
                sparse_polynomial_print(onereduce, p);
                sparse_polynomial_print(wpoly, p);
                sparse_polynomial_addincr(onereduce, wpoly);
                sparse_polynomial_print(onereduce, p);
                notreduced = 0;
                oneindex = onereduce->indexes;
                onecoeff = onereduce->coefficients;
            }
            I++;
        }
        
        if (notreduced) {
            /*Incrémenter n*/
            n++;
            ndiv2 = n>>1;
            
            if (n >= nmax) {
                printf("Reallocation failed.\n");
                break;
            }
            
            //printf("n: %lu \t Iminimal: %lu \t sfsize: %lu\n", n-1, Iminimal, sfsize);
            
            /*Augmenter partpolys en ajoutant Xn1 + … + Xnp - 1 et sfpolys*/
            sfIbegin = sfszarray[ n % p];
            partIbegin = partsize;
            poly = sparse_polynomial_alloc(p+1, indlimbsize, (long)1);
            mpn_zero(ind, indlimbsize);
            sparse_polynomial_add_monom(poly, poly, ind, indlimbsize, (long)-1);
            for (i=0; i<p; i++) {
                sparse_polynomial_add_monom(poly, poly, partindex, indlimbsize, (long)1);
                
                /*Augmenter sfpolys selon la huche i*/
                *sfindex1 = (mp_limb_t)1<<i;
                mpn_rshift(sfindex2, partindex, indlimbsize, p);
                for (I=0; I<ndiv2; I++) {
                    mpn_ior_n(ind, sfindex1, sfindex2, indlimbsize);
                    mpn_ior_n(sfindex0, ind, partindex, indlimbsize);
                    sfpoly = sparse_polynomial_alloc(1, indlimbsize, 3);
                    sparse_polynomial_add_monom(sfpoly, sfpoly, sfindex0, indlimbsize, (long)1);
                    sfpolys[sfsize] = sfpoly;
                    printf("%lu\n", sfsize);
                    sfsize++;
                    mpn_lshift(sfindex1, sfindex1, indlimbsize, p);
                    mpn_rshift(sfindex2, sfindex2, indlimbsize, p);
                    sparse_polynomial_print(sfpoly, p);
                }
                
                mpn_lshift(partindex, partindex, indlimbsize, (long)1);
            }
            partpolys[partsize] = poly;
            partsize++;
            //Iminimal = sfsize;
            sfszarray[n % p] = sfsize;
        }
    }
    
    //printf("n: %lu \t Iminimal: %lu \t sfsize: %lu\n", n, Iminimal, sfsize);
    
    /*Déallocation*/
    for (I=0; I<sfsize; I++) {
        sparse_polynomial_free(sfpolys[I]);
    }
    free(sfpolys);
    
    for (I=0; I<partsize; I++) {
        sparse_polynomial_free(partpolys[I]);
    }
    free(partpolys);
    
    free(ind);
    free(partindex);
    free(sfindex0);
    free(sfindex1);
    free(sfindex2);
    
    free(sfszarray);
    
    sparse_polynomial_free(onereduce);
    sparse_polynomial_free(wpoly);
    //mpq_clear(coeff);
    
    return n-1;
}

unsigned long schurNumberWithOrderReduction1(unsigned int p) {
    char notreduced;
    unsigned long n, ndiv2, nmodbpl, ndivbpl, nmax;
    size_t I;
    unsigned int i;
    size_t indlimbsize;
    mp_limb_t *ind, *polyindex, *oneindex;
    mp_limb_t *partindex;
    mp_limb_t *sfindex0, *sfindex1, *sfindex2;
    sparse_polynomial *onereduce;
    mpq_t one;
    mpq_t negone;
    mpq_t coeff;
    mpq_t *onecoeff;
    size_t sfsize, sfalloc;
    sparse_polynomial **sfpolys;  //Tableau des polynômes traduisant le caractère sans-somme
    size_t partsize, partalloc;
    sparse_polynomial **partpolys;//Tableau des polynômes traduisant le caractère partitionnant
    sparse_polynomial *poly, *sfpoly, *wpoly;
    
    /*Allocation*/
    nmax = p * mp_bits_per_limb;
    sfalloc = (nmax>>1) * (nmax>>1);
    sfsize = 0;
    sfpolys = calloc(sfalloc, sizeof(sparse_polynomial *));
    
    partalloc = p * nmax;
    partsize = 0;
    partpolys = calloc(partalloc, sizeof(sparse_polynomial *));
    
    indlimbsize = p;
    ind = calloc(p, sizeof(mp_limb_t));
    partindex = calloc(p, sizeof(mp_limb_t));
    sfindex0 = calloc(p, sizeof(mp_limb_t));
    sfindex1 = calloc(p, sizeof(mp_limb_t));
    sfindex2 = calloc(p, sizeof(mp_limb_t));
    wpoly = sparse_polynomial_alloc(1<<(2*p), indlimbsize, 1);
    mpq_init(coeff);
    
    /*Assignation de 1 à onereduce*/
    mpq_init(one);
    mpq_set_si(one, 1, 1);
    mpq_init(negone);
    mpq_set_si(negone, -1, 1);
    onereduce = sparse_polynomial_alloc(1<<(2*p), indlimbsize, 0);
    sparse_polynomial_add_monom(onereduce, onereduce, ind, indlimbsize, &one);
    oneindex = onereduce->indexes;
    onecoeff = onereduce->coefficients;
    
    /*Remplissage initiale de sfpolys*/
    n = 2;
    ndiv2 = 1;
    for (i=0; i<p; i++) {
        sfpoly = sparse_polynomial_alloc(1, indlimbsize, 0);
        *ind = (1<<i) | (1<<(p+i));
        sparse_polynomial_add_monom(sfpoly, sfpoly, ind, indlimbsize, &one);
        //sparse_polynomial_print(sfpoly, p);
        sfpolys[i] = sfpoly;
    }
    sfsize = p;
    
    /*Remplissage initiale de partpolys*/
    *partpolys = sparse_polynomial_alloc(p+1, indlimbsize, 0);
    partpolys[1] = sparse_polynomial_alloc(p+1, indlimbsize, 0);
    *partindex = (mp_limb_t)0;
    sparse_polynomial_add_monom(*partpolys, *partpolys, partindex, indlimbsize, &negone);
    sparse_polynomial_add_monom(partpolys[1], partpolys[1], partindex, indlimbsize, &negone);
    *partindex = (mp_limb_t)1;
    for (i=0; i<p; i++) {
        sparse_polynomial_add_monom(*partpolys, *partpolys, partindex, indlimbsize, &one);
        mpn_lshift(partindex, partindex, indlimbsize, 1);
    }
    for (i=0; i<p; i++) {
        sparse_polynomial_add_monom(partpolys[1], partpolys[1], partindex, indlimbsize, &one);
        mpn_lshift(partindex, partindex, indlimbsize, 1);
    }
    partsize = 2;
    
    while (onereduce->coeffsize) {
        /*Tant que onereduce ≠ 0*/
        notreduced = 1;
        
        /*Réduire onereduce par sfpolys*/
        I = 0;
        poly = *sfpolys;
        while (I < sfsize && notreduced) {
            /*Tester si sfpolys[I] réduit onereduce*/
            polyindex = poly->indexes;
            mpn_and_n(ind, polyindex, oneindex, indlimbsize);
            if (!mpn_cmp(ind, polyindex, indlimbsize)) {
                /*Retrancher à onereduce*/
                mpq_neg(coeff, *onecoeff);
                mpn_copyd(ind, oneindex, indlimbsize);
                //printf("Réduction par %lu:\n", *ind);
                //sparse_polynomial_print(onereduce, p);
                sparse_polynomial_add_monom(onereduce, onereduce, ind, indlimbsize, &coeff);
                //sparse_polynomial_print(onereduce, p);
                notreduced = 0;
                oneindex = onereduce->indexes;
                onecoeff = onereduce->coefficients;
            }
            I++;
            poly = sfpolys[I];
        }
        
        /*Réduire onereduce par partpolys*/
        I = 0;
        poly = *partpolys;
        while (I < partsize && notreduced) {
            /*Tester si partpolys[I] réduit onereduce*/
            sparse_polynomial_mul_monom(wpoly, poly, oneindex, indlimbsize, onecoeff);
            if (!mpn_cmp(wpoly->indexes, oneindex, indlimbsize)) {
                //printf("Réduction selon une partition:\n");
                //sparse_polynomial_print(onereduce, p);
                //sparse_polynomial_print(wpoly, p);
                sparse_polynomial_addincr(onereduce, wpoly);
                //sparse_polynomial_print(onereduce, p);
                notreduced = 0;
                oneindex = onereduce->indexes;
                onecoeff = onereduce->coefficients;
            }
            I++;
            poly = partpolys[I];
        }
        
        if (notreduced) {
            /*Incrémenter n*/
            nmodbpl = n%mp_bits_per_limb;
            ndivbpl = n>>GMP_2EXP_NUMB_BITS;
            n++;
            ndiv2 = n>>1;
            
            if (n >= nmax) {
                printf("Reallocation failed.\n");
                break;
            }
            
            /*Augmenter partpolys en ajoutant Xn1 + … + Xnp - 1 et sfpolys*/
            poly = sparse_polynomial_alloc(p+1, indlimbsize, 1);
            mpn_zero(ind, indlimbsize);
            sparse_polynomial_add_monom(poly, poly, ind, indlimbsize, &negone);
            for (i=0; i<p; i++) {
                sparse_polynomial_add_monom(poly, poly, partindex, indlimbsize, &one);
                
                /*Augmenter sfpolys selon la huche i*/
                *sfindex1 = (mp_limb_t)1<<i;
                mpn_rshift(sfindex2, partindex, indlimbsize, p);
                for (I=0; I<ndiv2; I++) {
                    mpn_ior_n(ind, sfindex1, sfindex2, indlimbsize);
                    mpn_ior_n(sfindex0, ind, partindex, indlimbsize);
                    sfpoly = sparse_polynomial_alloc(1, indlimbsize, 3);
                    sparse_polynomial_add_monom(sfpoly, sfpoly, sfindex0, indlimbsize, &one);
                    sfpolys[sfsize] = sfpoly;
                    sfsize++;
                    mpn_lshift(sfindex1, sfindex1, indlimbsize, p);
                    mpn_rshift(sfindex2, sfindex2, indlimbsize, p);
                    //sparse_polynomial_print(sfpoly, p);
                }
                
                mpn_lshift(partindex, partindex, indlimbsize, 1);
            }
            partpolys[partsize] = poly;
            partsize++;
        }
    }
    
    /*Déallocation*/
    for (I=0; I<sfsize; I++) {
        sparse_polynomial_free(sfpolys[I]);
    }
    free(sfpolys);
    
    for (I=0; I<partsize; I++) {
        sparse_polynomial_free(partpolys[I]);
    }
    free(partpolys);
    
    free(ind);
    free(partindex);
    free(sfindex0);
    free(sfindex1);
    free(sfindex2);
    
    sparse_polynomial_free(onereduce);
    sparse_polynomial_free(wpoly);
    mpq_clear(one);
    mpq_clear(negone);
    mpq_clear(coeff);
    
    return n-1;
}

