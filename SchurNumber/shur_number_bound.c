//
//  shur_number_bound.c
//  SchurNumber
//
//  Created by rubis on 04/03/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include <stdlib.h>
#include <gmp.h>

unsigned long schurNumberIterBound(unsigned long pmax, unsigned long *nbests, mpz_t iternum) {
    unsigned long n;
    unsigned long i;
    unsigned long p;
    char isSumFree;
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t *set;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    mp_size_t limballoc;    // Nombre de limbes alloué à chaque ensemble de sfpartition et à work0 et work1
    mp_size_t limbsize;     // Nombre de limbes utilisés par les ensembles de sfpartition
    mp_size_t wlimbsize;    // Nombre de limbes utilisés par work0 et work1
    mp_size_t wblimbsize;   // Idem mais pour tester nbests[p-1]
    unsigned int nmodbpl;
    unsigned int shift;
    unsigned int bshift;
    mp_limb_t mask;
    unsigned int *nbestslegalindexes;
    unsigned int nbestlegalindexes;
    unsigned int mayBeGood;
    
    /*Initialisation*/
    limballoc = pmax;
    work0 = calloc(limballoc, sizeof(mp_limb_t));
    work1 = calloc(limballoc, sizeof(mp_limb_t));
    sfpartition = calloc(pmax, sizeof(mp_limb_t *));  //Tableau contenant la partition
    sfpartitioninvert = calloc(pmax, sizeof(mp_limb_t *));    //Tableau contenant les ensembles "inverses" de la partition
    nbestslegalindexes = calloc(pmax, sizeof(int)); // Tableau dont les entiers indiquent les huches susceptibles de contenir nbest[p-1]
    for (i=0; i<pmax; i++) {
        sfpartition[i] = calloc(limballoc, sizeof(mp_limb_t));
        sfpartitioninvert[i] = calloc(limballoc, sizeof(mp_limb_t));
        nbests[i] = i+1;
        nbestslegalindexes[i] = 1;
    }
    
    /*Création de la première partition {{1}}*/
    n = 1;  // Taille de l'intervalle
    i = 0;  // Huche où placer l'entier suivant
    p = 1;  // Nombre de huches non vides
    *sfpartition[0] = (mp_limb_t)1;
    sfpartitioninvert[0][limballoc - 1] = (mp_limb_t)1<<(mp_bits_per_limb-1);
    limbsize = 1;
    wlimbsize = 1;    // Taille en limbe de la seconde moitié
    wblimbsize = 1;
    nmodbpl = 1;
    shift = mp_bits_per_limb - nmodbpl;
    bshift = mp_bits_per_limb - 1;
    nbestlegalindexes = *nbestslegalindexes;
    
    /*Itération jusqu'à énumérer toutes les partions sans-somme à au plus pmax huches*/
    while (n>0) {
        // Placer n+1 dans une des huches en conservant la sans-sommité
        while (i < p) {
            // Tester si l'ensemble obtenu en ajoutant n à la huche i est sans-somme
            mpz_add_ui(iternum, iternum, 1);
            set = sfpartitioninvert[i];
            mpn_rshift(work1, &set[limballoc - wlimbsize], wlimbsize, shift);
            set = sfpartition[i];
            mpn_and_n(work0, set, work1, wlimbsize);
            isSumFree = mpn_zero_p(work0, wlimbsize);
            if (isSumFree) {
                // Tester si la partition peut contenir nbest[p-1]+1
                if (/*nbestslegalindexes[p-1]*/nbestlegalindexes & (1<<i)) {
                    set = sfpartitioninvert[i];
                    mpn_rshift(work1, &set[limballoc - wblimbsize], wblimbsize, bshift);
                    set = sfpartition[i];
                    mpn_and_n(work0, set, work1, wblimbsize);
                    isSumFree = !mpn_zero_p(work0, wblimbsize);
                }
                mayBeGood = nbestlegalindexes & ~((unsigned int)isSumFree<<i);
                if (mayBeGood) {
                    // Ajouter n+1 à la huche i
                    mask = (mp_limb_t)1<<nmodbpl;
                    sfpartition[i][limbsize -1] |= mask;
                    mask = (mp_limb_t)1<<(shift - 1);
                    sfpartitioninvert[i][limballoc - limbsize] |= mask;
                    n++;
                    nmodbpl = n%mp_bits_per_limb;
                    shift = mp_bits_per_limb - nmodbpl;
                    wlimbsize = ((n+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                    nbestlegalindexes = mayBeGood;
                    if (n > nbests[p-1]) {
                        nbests[p-1] = n;
                        nbestlegalindexes = (INT_MAX>>(31-p))<<i;
                        //nbestslegalindexes[p-1] = INT_MAX>>(31-p);
                        wblimbsize = wlimbsize;
                        bshift = shift;
                        if (n >= mp_bits_per_limb*limballoc) {
                            limballoc *= 2;
                            work0 = realloc(work0, limballoc);
                            work1 = realloc(work1, limballoc);
                            for (i=0; i<pmax; i++) {
                                sfpartition[i] = realloc(sfpartition[i], limballoc);
                                sfpartitioninvert[i] = realloc(sfpartitioninvert[i], limballoc);
                            }
                        }
                    }
                    if (!nmodbpl) {
                        for (i=0; i<pmax; i++) {
                            sfpartition[i][limbsize] = (mp_limb_t)0;
                            sfpartitioninvert[i][limballoc - limbsize -1] = (mp_limb_t)0;
                        }
                        limbsize++;
                    }
                    i = 0;
                } else {
                    i++;
                }
            } else {
                i++;
            }
        }
        if (i == p && p < pmax) {
            // Cas particulier où il faut remplir une nouvelle huche
            mpz_add_ui(iternum, iternum, 1);
            mask = (mp_limb_t)1<<nmodbpl;
            sfpartition[i][limbsize -1] |= mask; // Ajoute n+1 à la huche i
            mask = (mp_limb_t)1<<(shift - 1);
            sfpartitioninvert[i][limballoc - limbsize] |= mask;
            n++;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            wlimbsize = ((n+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
            nbestslegalindexes[p-1] = nbestlegalindexes;
            if (n > nbests[p]) {
                nbests[p] = n;
                nbestlegalindexes = 1<<p;
                //nbestslegalindexes[p] = INT_MAX>>(31-p);
                wblimbsize = wlimbsize;
                bshift = shift;
                if (n >= 64*limballoc -1) {
                    limballoc *= 2;
                    work0 = realloc(work0, limballoc);
                    work1 = realloc(work1, limballoc);
                    for (i=0; i<pmax; i++) {
                        sfpartition[i] = realloc(sfpartition[i], limballoc);
                        sfpartitioninvert[i] = realloc(sfpartitioninvert[i], limballoc);
                    }
                }
            } else {
                wblimbsize = ((nbests[p]+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                bshift = mp_bits_per_limb - (nbests[p])%mp_bits_per_limb;
                nbestlegalindexes = 0;
                for (i=0; i<=p; i++) {
                    set = sfpartitioninvert[i];
                    mpn_rshift(work1, &set[limballoc - wblimbsize], wblimbsize, bshift);
                    set = sfpartition[i];
                    mpn_and_n(work0, set, work1, wblimbsize);
                    nbestlegalindexes |= mpn_zero_p(work0, wblimbsize)<<i;
                }
            }
            if (!nmodbpl) {
                for (i=0; i<pmax; i++) {
                    sfpartition[i][limbsize] = (mp_limb_t)0;
                    sfpartitioninvert[i][limballoc - limbsize -1] = (mp_limb_t)0;
                }
                limbsize++;
            }
            i = 0;
            p++;
        } else {
            // Dépiler
            if (!nmodbpl) {
                limbsize--;
            }
            n--;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            wlimbsize = ((n+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
            mask = (mp_limb_t)1<<nmodbpl;
            i = 0;
            while (i < p && !(sfpartition[i][limbsize-1] & mask)) {
                i++;
            }
            sfpartition[i][limbsize-1] ^= mask;
            mask = (mp_limb_t)1<<(shift - 1);
            sfpartitioninvert[i][limballoc - limbsize] ^= mask;
//            if (!(nbestslegalindexes[p-1] & (1<<i))) {
//                set = sfpartitioninvert[i];
//                mpn_rshift(work1, &set[limballoc - wblimbsize], wblimbsize, bshift);
//                set = sfpartition[i];
//                mpn_and_n(work0, set, work1, wblimbsize);
//                isSumFree = mpn_zero_p(work0, wblimbsize);
//                nbestslegalindexes[p-1] |= isSumFree<<i;
//            }
            if (i == p-1) {
                if (mpn_zero_p(sfpartition[i], limbsize)) {
                    p--;
                    nbestlegalindexes = nbestslegalindexes[p-1];
                    wblimbsize = ((nbests[p-1]+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                    bshift = mp_bits_per_limb - (nbests[p-1])%mp_bits_per_limb;
                } else if (!(nbestlegalindexes & (1<<i))) {
                    set = sfpartitioninvert[i];
                    mpn_rshift(work1, &set[limballoc - wblimbsize], wblimbsize, bshift);
                    set = sfpartition[i];
                    mpn_and_n(work0, set, work1, wblimbsize);
                    isSumFree = mpn_zero_p(work0, wblimbsize);
                    nbestlegalindexes |= isSumFree<<i;
                }
            } else if (!(nbestlegalindexes & (1<<i))) {
                set = sfpartitioninvert[i];
                mpn_rshift(work1, &set[limballoc - wblimbsize], wblimbsize, bshift);
                set = sfpartition[i];
                mpn_and_n(work0, set, work1, wblimbsize);
                isSumFree = mpn_zero_p(work0, wblimbsize);
                nbestlegalindexes |= isSumFree<<i;
            }
            i++;
        }
    }
    
    /*Nettoyage*/
    for (i=0; i<pmax; i++) {
        free(sfpartition[i]);
        free(sfpartitioninvert[i]);
    }
    free(sfpartition);
    free(sfpartitioninvert);
    free(work0);
    free(work1);
    free(nbestslegalindexes);
    return nbests[pmax-1];
}

unsigned long schurNumberIterBound2(unsigned long pmax, unsigned long *nbests, mpz_t iternum) {
    unsigned long n;
    unsigned long i;
    unsigned long p;
    char isSumFree;
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t *set;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    mp_size_t limballoc;    // Nombre de limbes alloué à chaque ensemble de sfpartition et à work0 et work1
    mp_size_t limbsize;     // Nombre de limbes utilisés par les ensembles de sfpartition
    mp_size_t wlimbsize;    // Nombre de limbes utilisés par work0 et work1
    mp_size_t wblimbsize;   // Idem mais pour tester nbests[p-1]
    unsigned int nmodbpl;
    unsigned int shift;
    unsigned int bshift;
    mp_limb_t mask;
    unsigned int *nbestslegalindexes;
    unsigned int nbestlegalindexes;
    unsigned int mayBeGood;
    
    /*Initialisation*/
    limballoc = pmax;
    work0 = calloc(limballoc, sizeof(mp_limb_t));
    work1 = calloc(limballoc, sizeof(mp_limb_t));
    sfpartition = calloc(pmax, sizeof(mp_limb_t *));  //Tableau contenant la partition
    sfpartitioninvert = calloc(pmax, sizeof(mp_limb_t *));    //Tableau contenant les ensembles "inverses" de la partition
    nbestslegalindexes = calloc(pmax, sizeof(int)); // Tableau dont les entiers indiquent les huches susceptibles de contenir nbest[p-1]
    for (i=0; i<pmax; i++) {
        sfpartition[i] = calloc(limballoc, sizeof(mp_limb_t));
        sfpartitioninvert[i] = calloc(limballoc, sizeof(mp_limb_t));
        nbests[i] = i+1;
        nbestslegalindexes[i] = 1;
    }
    
    /*Création de la première partition {{1}}*/
    n = 1;  // Taille de l'intervalle
    i = 0;  // Huche où placer l'entier suivant
    p = 1;  // Nombre de huches non vides
    *sfpartition[0] = (mp_limb_t)1;
    sfpartitioninvert[0][limballoc - 1] = (mp_limb_t)1<<(mp_bits_per_limb-1);
    limbsize = 1;
    wlimbsize = 1;    // Taille en limbe de la seconde moitié
    wblimbsize = 1;
    nmodbpl = 1;
    shift = mp_bits_per_limb - nmodbpl;
    bshift = mp_bits_per_limb - 1;
    nbestlegalindexes = *nbestslegalindexes;
    
    /*Itération jusqu'à énumérer toutes les partions sans-somme à au plus pmax huches*/
    while (n>0) {
        // Placer n+1 dans une des huches en conservant la sans-sommité
        while (i < p) {
            // Tester si la partition peut contenir nbest[p-1]+1
            mpz_add_ui(iternum, iternum, 1);
            isSumFree = 0;
            if (nbestlegalindexes & (1<<i)) {
                set = sfpartitioninvert[i];
                mpn_rshift(work1, &set[limballoc - wblimbsize], wblimbsize, bshift);
                set = sfpartition[i];
                mpn_and_n(work0, set, work1, wblimbsize);
                isSumFree = !mpn_zero_p(work0, wblimbsize);
            }
            mayBeGood = nbestlegalindexes & ~((unsigned int)isSumFree<<i);
            if (mayBeGood) {
                // Tester si l'ensemble obtenu en ajoutant n à la huche i est sans-somme
                set = sfpartitioninvert[i];
                mpn_rshift(work1, &set[limballoc - wlimbsize], wlimbsize, shift);
                set = sfpartition[i];
                mpn_and_n(work0, set, work1, wlimbsize);
                isSumFree = mpn_zero_p(work0, wlimbsize);
                if (isSumFree) {
                    // Ajouter n+1 à la huche i
                    mask = (mp_limb_t)1<<nmodbpl;
                    sfpartition[i][limbsize -1] |= mask;
                    mask = (mp_limb_t)1<<(shift - 1);
                    sfpartitioninvert[i][limballoc - limbsize] |= mask;
                    n++;
                    nmodbpl = n%mp_bits_per_limb;
                    shift = mp_bits_per_limb - nmodbpl;
                    wlimbsize = ((n+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                    nbestlegalindexes = mayBeGood;
                    if (n > nbests[p-1]) {
                        nbests[p-1] = n;
                        nbestlegalindexes = (INT_MAX>>(31-p))<<i;
                        wblimbsize = wlimbsize;
                        bshift = shift;
                        if (n >= mp_bits_per_limb*limballoc) {
                            limballoc *= 2;
                            work0 = realloc(work0, limballoc);
                            work1 = realloc(work1, limballoc);
                            for (i=0; i<pmax; i++) {
                                sfpartition[i] = realloc(sfpartition[i], limballoc);
                                sfpartitioninvert[i] = realloc(sfpartitioninvert[i], limballoc);
                            }
                        }
                    }
                    if (!nmodbpl) {
                        for (i=0; i<pmax; i++) {
                            sfpartition[i][limbsize] = (mp_limb_t)0;
                            sfpartitioninvert[i][limballoc - limbsize -1] = (mp_limb_t)0;
                        }
                        limbsize++;
                    }
                    i = 0;
                } else {
                    i++;
                }
            } else {
                i++;
            }
        }
        if (i == p && p < pmax) {
            // Cas particulier où il faut remplir une nouvelle huche
            mpz_add_ui(iternum, iternum, 1);
            mask = (mp_limb_t)1<<nmodbpl;
            sfpartition[i][limbsize -1] |= mask; // Ajoute n+1 à la huche i
            mask = (mp_limb_t)1<<(shift - 1);
            sfpartitioninvert[i][limballoc - limbsize] |= mask;
            n++;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            wlimbsize = ((n+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
            nbestslegalindexes[p-1] = nbestlegalindexes;
            if (n > nbests[p]) {
                nbests[p] = n;
                nbestlegalindexes = 1<<p;
                wblimbsize = wlimbsize;
                bshift = shift;
                if (n >= 64*limballoc -1) {
                    limballoc *= 2;
                    work0 = realloc(work0, limballoc);
                    work1 = realloc(work1, limballoc);
                    for (i=0; i<pmax; i++) {
                        sfpartition[i] = realloc(sfpartition[i], limballoc);
                        sfpartitioninvert[i] = realloc(sfpartitioninvert[i], limballoc);
                    }
                }
            } else {
                wblimbsize = ((nbests[p]+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                bshift = mp_bits_per_limb - (nbests[p])%mp_bits_per_limb;
                nbestlegalindexes = 0;
                for (i=0; i<=p; i++) {
                    set = sfpartitioninvert[i];
                    mpn_rshift(work1, &set[limballoc - wblimbsize], wblimbsize, bshift);
                    set = sfpartition[i];
                    mpn_and_n(work0, set, work1, wblimbsize);
                    nbestlegalindexes |= mpn_zero_p(work0, wblimbsize)<<i;
                }
            }
            if (!nmodbpl) {
                for (i=0; i<pmax; i++) {
                    sfpartition[i][limbsize] = (mp_limb_t)0;
                    sfpartitioninvert[i][limballoc - limbsize -1] = (mp_limb_t)0;
                }
                limbsize++;
            }
            i = 0;
            p++;
        } else {
            // Dépiler
            if (!nmodbpl) {
                limbsize--;
            }
            n--;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            wlimbsize = ((n+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
            mask = (mp_limb_t)1<<nmodbpl;
            i = 0;
            while (i < p && !(sfpartition[i][limbsize-1] & mask)) {
                i++;
            }
            sfpartition[i][limbsize-1] ^= mask;
            mask = (mp_limb_t)1<<(shift - 1);
            sfpartitioninvert[i][limballoc - limbsize] ^= mask;
            if (i == p-1) {
                if (mpn_zero_p(sfpartition[i], limbsize)) {
                    p--;
                    nbestlegalindexes = nbestslegalindexes[p-1];
                    wblimbsize = ((nbests[p-1]+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                    bshift = mp_bits_per_limb - (nbests[p-1])%mp_bits_per_limb;
                } else if (!(nbestlegalindexes & (1<<i))) {
                    set = sfpartitioninvert[i];
                    mpn_rshift(work1, &set[limballoc - wblimbsize], wblimbsize, bshift);
                    set = sfpartition[i];
                    mpn_and_n(work0, set, work1, wblimbsize);
                    isSumFree = mpn_zero_p(work0, wblimbsize);
                    nbestlegalindexes |= isSumFree<<i;
                }
            } else if (!(nbestlegalindexes & (1<<i))) {
                set = sfpartitioninvert[i];
                mpn_rshift(work1, &set[limballoc - wblimbsize], wblimbsize, bshift);
                set = sfpartition[i];
                mpn_and_n(work0, set, work1, wblimbsize);
                isSumFree = mpn_zero_p(work0, wblimbsize);
                nbestlegalindexes |= isSumFree<<i;
            }
            i++;
        }
    }
    
    /*Nettoyage*/
    for (i=0; i<pmax; i++) {
        free(sfpartition[i]);
        free(sfpartitioninvert[i]);
    }
    free(sfpartition);
    free(sfpartitioninvert);
    free(work0);
    free(work1);
    free(nbestslegalindexes);
    return nbests[pmax-1];
}
