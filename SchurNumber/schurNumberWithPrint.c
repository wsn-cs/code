//
//  schurNumberWithPrint.c
//  SchurNumber
//
//  Created by rubis on 17/04/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#ifndef SCHUR_NUMBER

#define SCHUR_NUMBER

#include <stdlib.h>
#include <gmp.h>

#if GMP_NUMB_BITS == 64
#define GMP_2EXP_NUMB_BITS 6
#elif GMP_NUMB_BITS == 128
#define GMP_2EXP_NUMB_BITS 7
#elif GMP_NUMB_BITS == 32
#define GMP_2EXP_NUMB_BITS 5
#endif

#endif

#include "rscan1.c"

void printPartition(unsigned long p, unsigned long n, mp_limb_t **partition);

unsigned long schurNumberIterWithPrintWithUnstack(unsigned long pmax, unsigned long *nbests, mpz_t iternum) {
    /*
     Cette fonction calcule successivement les nombres de Schur S(p) pour p<= pmax.
     Elle remplit le tableau nbests.
     
     La partition est représentée comme un tableau de grands entiers sfpartition.
     Chaque grand entier représente un ensemble grâce à ses bits: si le bit k vaut 1,
     alors k appartient à l'ensemble; sinon il n'y appartient pas.
     
     Pour tester si il est possible d'ajouter n+1 à l'ensemble P sans-somme,
     on effectue un & entre les [(n+1)/2] premiers bits avec les [(n+1)/2] derniers.
     Si le résultat est 0, il est possible d'ajouter n+1.
     
     Lorsqu'un entier n n'a pu être inséré dans aucune huches, on revient au sup sur les huches
     du plus petit n'≥[(n+1)/2] tel que n'  et n-n' appartiennent simultanément à cette huche.
     */
    unsigned long n;
    unsigned long nblocking;
    unsigned long nblockingmax;
    unsigned long nmed;
    unsigned long i, j;
    unsigned long iblocking;
    unsigned long p;
    char isSumFree;
    char isAppendable;  //Vrai si il est possible de former une partition sans-somme contenant n
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t *set;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    mp_size_t limballoc;    // Nombre de limbes alloué à chaque ensemble de sfpartition et à work0 et work1
    mp_size_t limbsize;     // Nombre de limbes utilisés par les ensembles de sfpartition
    mp_size_t wlimbsize;    // Nombre de limbes utilisés par work0
    unsigned int nmodbpl;
    unsigned int shift;
    mp_limb_t mask;
    mp_limb_t mask2;
    
    mp_limb_t **sfpartitionsbest;  // Tableau contenant les partitions les plus grandes
    size_t sfpartitionbestcount;
    size_t sfpartitionbestsize;
    
    /*Initialisation*/
    limballoc = pmax;
    work0 = calloc(limballoc, sizeof(mp_limb_t));
    work1 = calloc(limballoc, sizeof(mp_limb_t));
    sfpartition = calloc(pmax, sizeof(mp_limb_t *));  //Tableau contenant la partition
    sfpartitioninvert = calloc(pmax, sizeof(mp_limb_t *));    //Tableau contenant les ensembles "inverses" de la partition
    for (i=0; i<pmax; i++) {
        sfpartition[i] = calloc(limballoc, sizeof(mp_limb_t));
        sfpartitioninvert[i] = calloc(limballoc, sizeof(mp_limb_t));
        nbests[i] = i+1;
    }
    sfpartitionbestsize = 64*pmax;
    sfpartitionbestcount = 0;
    sfpartitionsbest = calloc(sfpartitionbestsize, pmax * sizeof(mp_limb_t *));
    for (j=0; j<sfpartitionbestsize; j++) {
        for (i=0; i<pmax; i++) {
            sfpartitionsbest[pmax*j + i] = calloc(limballoc, sizeof(mp_limb_t));
        }
    }
    
    /*Création de la première partition {{1}}*/
    n = 1;  // Taille de l'intervalle
    i = 0;  // Huche où placer l'entier suivant
    p = 1;  // Nombre de huches non vides
    *sfpartition[0] = (mp_limb_t)1;
    sfpartitioninvert[0][limballoc - 1] = (mp_limb_t)1<<(mp_bits_per_limb-1);
    limbsize = 1;
    wlimbsize = 1;    // Taille en limbe de la seconde moitié
    nmodbpl = 1;
    shift = mp_bits_per_limb - nmodbpl;
    nmed = 0;
    nblockingmax = 0;
    isAppendable = 0;
    
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
                // Ajouter n+1 à la huche i
                mask = (mp_limb_t)1<<nmodbpl;
                sfpartition[i][limbsize -1] |= mask;
                mask = (mp_limb_t)1<<(shift - 1);
                sfpartitioninvert[i][limballoc - limbsize] |= mask;
                n++;
                nmodbpl = n%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                nmed = (n+1)>>1;
                wlimbsize = (nmed>>GMP_2EXP_NUMB_BITS) + 1;
                nmed--;
                nblockingmax = 0;
                if (n > nbests[p-1]) {
                    nbests[p-1] = n;
                    
                    if (p == pmax) {
                        sfpartitionbestcount = 0;
                    }
                    
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
                isAppendable = 0;
            } else {
                nblocking = n - mpn_rscan1(work0, nmed);
                if (nblocking > nblockingmax) {
                    nblockingmax = nblocking;
                    iblocking = i;
                }
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
            nmed = (n+1)>>1;
            wlimbsize = (nmed>>GMP_2EXP_NUMB_BITS) + 1;
            nmed--;
            nblockingmax = 0;
            if (n > nbests[p]) {
                nbests[p] = n;
                
                if (p+1 == pmax) {
                    sfpartitionbestcount = 0;
                }
                
                if (n >= 64*limballoc -1) {
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
            isAppendable = 0;
            p++;
        } else if(isAppendable) {
            // Sauvegarder éventuellement la partition
            if ((p == pmax) && (n == nbests[pmax - 1])) {
                for (i=0; i<pmax; i++) {
                    mpn_copyi(sfpartitionsbest[pmax*sfpartitionbestcount + i], sfpartition[i], limbsize);
                }
                sfpartitionbestcount++;
                if (sfpartitionbestcount == sfpartitionbestsize) {
                    sfpartitionbestsize *= 2;
                    sfpartitionsbest = realloc(sfpartitionsbest, sfpartitionbestsize*pmax);
                }
            }
            // Dépiler d'un entier
            if (!nmodbpl) {
                limbsize--;
            }
            n--;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            nmed = (n+1)>>1;
            wlimbsize = (nmed>>GMP_2EXP_NUMB_BITS) + 1;
            nmed--;
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
                }
            }
            i++;
            isAppendable = 1;
        } else {
            // Sauvegarder éventuellement la partition
            if ((p == pmax) && (n == nbests[pmax - 1])) {
                for (i=0; i<pmax; i++) {
                    mpn_copyi(sfpartitionsbest[pmax*sfpartitionbestcount + i], sfpartition[i], limbsize);
                }
                sfpartitionbestcount++;
                if (sfpartitionbestcount == sfpartitionbestsize) {
                    sfpartitionbestsize *= 2;
                    sfpartitionsbest = realloc(sfpartitionsbest, sfpartitionbestsize*pmax);
                }
            }
            // Dépiler jusqu'à revenir à nblockingmax
            n = nblockingmax;
            limbsize = (n>>GMP_2EXP_NUMB_BITS) + 1;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            nmed = (n+1)>>1;
            wlimbsize = (nmed>>GMP_2EXP_NUMB_BITS) + 1;
            nmed--;
            mask = (mp_limb_t)GMP_NUMB_MAX >> shift;
            mask2 = (mp_limb_t)GMP_NUMB_MAX << (shift);
            for (i=0; i<pmax; i++) {
                sfpartition[i][limbsize-1] &= mask;
                sfpartitioninvert[i][limballoc - limbsize] &= mask2;
            }
            while (mpn_zero_p(sfpartition[p-1], limbsize)) {
                p--;
            }
            i = iblocking + 1;
            isAppendable = 1;
        }
    }
    
    for (j=0; j<sfpartitionbestsize; j++) {
        /*Afficher les partitions optimales puis les libérer*/
        printPartition(pmax, nbests[pmax-1], &sfpartitionsbest[j*pmax]);
        for (i=0; i<pmax; i++) {
            free(sfpartitionsbest[pmax*j + i]);
        }
    }
    free(sfpartitionsbest);
    
    /*Nettoyage*/
    for (i=0; i<pmax; i++) {
        free(sfpartition[i]);
        free(sfpartitioninvert[i]);
    }
    free(sfpartition);
    free(sfpartitioninvert);
    free(work0);
    free(work1);
    return nbests[pmax-1];
}

unsigned long schurNumberIterWithPrint(unsigned long pmax, unsigned long *nbests, mpz_t iternum) {
    /*
     Cette fonction calcule successivement les nombres de Schur S(p) pour p<= pmax.
     Elle remplit le tableau nbests.
     
     La partition est représentée comme un tableau de grands entiers sfpartition.
     Chaque grand entier représente un ensemble grâce à ses bits: si le bit k vaut 1,
     alors k appartient à l'ensemble; sinon il n'y appartient pas.
     
     Pour tester si il est possible d'ajouter n+1 à l'ensemble P sans-somme,
     on effectue un & entre les [(n+1)/2] premiers bits avec les [(n+1)/2] derniers.
     Si le résultat est 0, il est possible d'ajouter n+1.
     */
    unsigned long n;
    unsigned long i;
    unsigned long p;
    char isSumFree;
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t *set;
    mp_limb_t *lastlimbptr;
    mp_limb_t *lastlimbinvertptr;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    mp_size_t limballoc;    // Nombre de limbes alloué à chaque ensemble de sfpartition et à work0 et work1
    mp_size_t limbsize;     // Nombre de limbes utilisés par les ensembles de sfpartition
    mp_size_t wlimbsize;    // Nombre de limbes utilisés par work0
    unsigned int nmodbpl;
    unsigned int shift;
    mp_limb_t mask;
    
    /*Initialisation*/
    limballoc = pmax;
    work0 = calloc(limballoc, sizeof(mp_limb_t));
    work1 = calloc(limballoc, sizeof(mp_limb_t));
    sfpartition = calloc(pmax, sizeof(mp_limb_t *));  //Tableau contenant la partition
    sfpartitioninvert = calloc(pmax, sizeof(mp_limb_t *));    //Tableau contenant les ensembles "inverses" de la partition
    for (i=0; i<pmax; i++) {
        sfpartition[i] = calloc(limballoc, sizeof(mp_limb_t));
        sfpartitioninvert[i] = calloc(limballoc, sizeof(mp_limb_t));
        nbests[i] = i+1;
    }
    
    /*Création de la première partition {{1}}*/
    n = 1;  // Taille de l'intervalle
    i = 0;  // Huche où placer l'entier suivant
    p = 1;  // Nombre de huches non vides
    (*sfpartition)[0] = (mp_limb_t)1;
    (*sfpartitioninvert)[limballoc - 1] = (mp_limb_t)1<<(mp_bits_per_limb-1);
    limbsize = 1;
    wlimbsize = 1;    // Taille en limbe de la seconde moitié
    nmodbpl = 1;
    shift = mp_bits_per_limb - nmodbpl;
    
    /*Itération jusqu'à énumérer toutes les partions sans-somme à au plus pmax huches*/
    while (n>0) {
        // Placer n+1 dans une des huches en conservant la sans-sommité
        while (i < p) {
            // Tester si l'ensemble obtenu en ajoutant n à la huche i est sans-somme
            mpz_add_ui(iternum, iternum, 1);
            set = sfpartition[i];
            lastlimbptr = &(set[limbsize -1]);
            lastlimbinvertptr = &(sfpartitioninvert[i][limballoc - wlimbsize]);
            mpn_rshift(work1, lastlimbinvertptr, wlimbsize, shift);
            mpn_and_n(work0, set, work1, wlimbsize);
            isSumFree = mpn_zero_p(work0, wlimbsize);
            if (isSumFree) {
                // Ajouter n+1 à la huche i
                mask = (mp_limb_t)1<<nmodbpl;
                *lastlimbptr |= mask;
                mask = (mp_limb_t)1<<(shift - 1);
                *lastlimbinvertptr |= mask;
                n++;
                nmodbpl = n%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                wlimbsize = ((n+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                if (n > nbests[p-1]) {
                    nbests[p-1] = n;
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
                //printPartition(p, n, sfpartition);
                //printPartition(p, mp_bits_per_limb*limballoc, sfpartitioninvert);
            } else {
                i++;
            }
        }
        if (i == p && p < pmax) {
            // Cas particulier où il faut remplir une nouvelle huche
            mpz_add_ui(iternum, iternum, 1);
            set = sfpartition[i];
            lastlimbptr = &(set[limbsize -1]);
            lastlimbinvertptr = &(sfpartitioninvert[i][limballoc - wlimbsize]);
            mask = (mp_limb_t)1<<nmodbpl;
            *lastlimbptr |= mask; // Ajoute n+1 à la huche i
            mask = (mp_limb_t)1<<(shift - 1);
            *lastlimbinvertptr |= mask;
            n++;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            wlimbsize = ((n+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
            if (n > nbests[p]) {
                nbests[p] = n;
                if (n >= GMP_NUMB_BITS*limballoc) {
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
            p++;
            //printPartition(p, n, sfpartition);
            //printPartition(p, mp_bits_per_limb*limballoc, sfpartitioninvert);
        } else {
            if (n == 40) {
                printPartition(p, n, sfpartition);
            }
            // Dépiler
            if (!nmodbpl) {
                limbsize--;
            }
            n--;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            wlimbsize = ((n)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
            mask = (mp_limb_t)1<<nmodbpl;
            i = 0;
            while (i < p && !(sfpartition[i][limbsize-1] & mask)) {
                i++;
            }
            set = sfpartition[i];
            set[limbsize-1] ^= mask;
            mask = (mp_limb_t)1<<(shift - 1);
            sfpartitioninvert[i][limballoc - limbsize] ^= mask;
            if (i == p-1) {
                if (mpn_zero_p(set, limbsize)) {
                    p--;
                }
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
    return nbests[pmax-1];
}
