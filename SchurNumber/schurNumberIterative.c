//
//  schurNumberIterative.c
//  SchurNumber
//
//  Created by Gabriel Merlin on 25/02/2019.
//

#ifndef SCHUR_NUMBER

#define SCHUR_NUMBER

#include <stdlib.h>
#include <gmp.h>
#include <limits.h>

#if GMP_NUMB_BITS == 64
#define GMP_2EXP_NUMB_BITS 6
#elif GMP_NUMB_BITS == 128
#define GMP_2EXP_NUMB_BITS 7
#elif GMP_NUMB_BITS == 32
#define GMP_2EXP_NUMB_BITS 5
#endif

#endif

//void printPartition(unsigned long p, mpz_t *partition) {
//    /*Affiche une partition.*/
//    unsigned long i;
//    unsigned long j;
//    printf("Partition:\n");
//    for (i=0; i<p; i++) {
//        printf("\t");
//        for (j=0; j < 64*(partition[i]->_mp_size); j++) {
//            if (mpz_tstbit(partition[i], j)) {
//                printf(" %lu", j+1);
//            }
//        }
//        printf("\n");
//    }
//}

unsigned long schurNumberIterative(unsigned long pmax, unsigned long *nbests, mpz_t iternum) {
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
    unsigned long limbnum;
    unsigned long rem;
    unsigned long limbmax;
    char isSumFree;
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t *set;
    mpz_t *sfpartition;
    mpz_t *sfpartitioninvert;
    
    /*Initialisation*/
    limbmax = pmax; //Taille de l'entier work0 et work1
    work0 = calloc(limbmax, sizeof(mp_limb_t));
    work1 = calloc(limbmax, sizeof(mp_limb_t));
    sfpartition = calloc(pmax, sizeof(mpz_t));  //Tableau contenant la partition
    sfpartitioninvert = calloc(pmax, sizeof(mpz_t));    //Tableau contenant les ensembles "inverses" de la partition
    for (i=0; i<pmax; i++) {
        mpz_init2(sfpartition[i], 64*pmax);
        mpz_init2(sfpartitioninvert[i], 64*pmax);
        nbests[i] = i+1;
    }
    
    /*Création de la première partition {{1}}*/
    n = 1;  // Taille de l'intervalle
    i = 0;  // Position où placer l'entier suivant
    p = 1;  // Nombre de huches non vides
    mpz_set_ui(sfpartition[0], 1);
    mpz_setbit(sfpartitioninvert[0], 64*(sfpartitioninvert[0]->_mp_alloc) -1);
    
    /*Itération jusqu'à énumérer toutes les partions sum-free à au plus pmax huches*/
    while (n>0) {
        // Placer n+1 dans une des huches en conservant la sans-sommité
        limbnum = ((n+1)>>7)+1;
        rem = 64 - (n%64);
        if (i < p) {
            // Tester si l'ensemble obtenu en ajoutant n à la huche i est sans-somme
            mpz_add_ui(iternum, iternum, 1);
            set = sfpartitioninvert[i]->_mp_d;
            mpn_copyd(work1, &set[sfpartitioninvert[i]->_mp_alloc - limbnum], limbnum);
            mpn_rshift(work0, work1, limbnum, rem%0x10000000);
            mpn_rshift(work1, work0, limbnum, ((unsigned long)rem>>16));
            set = sfpartition[i]->_mp_d;
            mpn_and_n(work0, set, work1, limbnum);
            isSumFree = mpn_zero_p(work0, limbnum);
            if (isSumFree) {
                mpz_setbit(sfpartition[i], n); // Ajoute n+1 à la huche i
                mpz_setbit(sfpartitioninvert[i], 64*(sfpartitioninvert[i]->_mp_alloc) - n - 1);
                n++;
                limbnum = ((n+1)>>7)+1;
                rem = 64 - (n%64);
                if (n > nbests[p-1]) {
                    nbests[p-1] = n;
                    if (n >= 64*limbmax -1) {
                        work0 = realloc(work0, 2*limbmax);
                        work1 = realloc(work1, 2*limbmax);
                    }
                }
                i = 0;
            } else {
                i++;
            }
        } else if (i == p && p < pmax) {
            // Cas particulier où il faut remplir une nouvelle huche
            mpz_add_ui(iternum, iternum, 1);
            mpz_setbit(sfpartition[i], n); // Ajoute n+1 à la huche i
            mpz_setbit(sfpartitioninvert[i], 64*(sfpartitioninvert[i]->_mp_alloc) - n - 1);
            n++;
            if (n > nbests[p]) {
                nbests[p] = n;
                if (n >= 64*limbmax -1) {
                    limbmax *= 2;
                    work0 = realloc(work0, limbmax);
                    work1 = realloc(work1, limbmax);
                }
            }
            i = 0;
            p++;
        } else {
            // Dépiler
            n --;
            i = 0;
            while (i < p && !mpz_tstbit(sfpartition[i], n)) {
                i++;
            }
            mpz_clrbit(sfpartition[i], n);
            mpz_clrbit(sfpartitioninvert[i], 64*(sfpartitioninvert[i]->_mp_alloc) - n -1);
            if (i == p-1) {
                if (!mpz_size(sfpartition[i])) {
                    p--;
                }
            }
            i++;
        }
    }
    
    /*Nettoyage*/
    for (i=0; i<pmax; i++) {
        mpz_clear(sfpartition[i]);
        mpz_clear(sfpartitioninvert[i]);
    }
    free(sfpartition);
    free(sfpartitioninvert);
    free(work0);
    free(work1);
    return nbests[pmax-1];
}

unsigned long schurNumberIterative1(unsigned long pmax, unsigned long *nbests, mpz_t iternum) {
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
    *sfpartition[0] = (mp_limb_t)1;
    sfpartitioninvert[0][limballoc - 1] = (mp_limb_t)1<<(mp_bits_per_limb-1);
    limbsize = 1;
    wlimbsize = 1;    // Taille en limbe de la seconde moitié
    nmodbpl = 1;
    shift = mp_bits_per_limb - nmodbpl;
    
    /*Itération jusqu'à énumérer toutes les partions sans-somme à au plus pmax huches*/
    while (n>0) {
        // Placer n+1 dans une des huches en conservant la sans-sommité
        if (i < p) {
            // Tester si l'ensemble obtenu en ajoutant n à la huche i est sans-somme
            mpz_add_ui(iternum, iternum, 1);
            set = sfpartitioninvert[i];
            mpn_copyd(work0, &set[limballoc - wlimbsize], wlimbsize);
            mpn_rshift(work1, work0, wlimbsize, shift);
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
                wlimbsize = ((n+1)>>(3*sizeof(mp_limb_t) + 1)) + 1;
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
            } else {
                i++;
            }
        } else if (i == p && p < pmax) {
            // Cas particulier où il faut remplir une nouvelle huche
            mpz_add_ui(iternum, iternum, 1);
            mask = (mp_limb_t)1<<nmodbpl;
            sfpartition[i][limbsize -1] |= mask; // Ajoute n+1 à la huche i
            mask = (mp_limb_t)1<<(shift - 1);
            sfpartitioninvert[i][limballoc - limbsize] |= mask;
            n++;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            wlimbsize = ((n+1)>>(3*sizeof(mp_limb_t) + 1)) + 1;
            if (n > nbests[p]) {
                nbests[p] = n;
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
            p++;
        } else {
            // Dépiler
            if (!nmodbpl) {
                limbsize--;
            }
            n--;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            wlimbsize = ((n+1)>>(3*sizeof(mp_limb_t) + 1)) + 1;
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

unsigned long schurNumberIterative2(unsigned long pmax, unsigned long *nbests, mpz_t iternum) {
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
    *sfpartition[0] = (mp_limb_t)1;
    sfpartitioninvert[0][limballoc - 1] = (mp_limb_t)1<<(mp_bits_per_limb-1);
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
            set = sfpartitioninvert[i];
            //mpn_copyd(work0, &set[limballoc - wlimbsize], wlimbsize);
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
            if (n > nbests[p]) {
                nbests[p] = n;
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

unsigned long schurNumberIterative21(unsigned long pmax, unsigned long *nbests, mpz_t iternum) {
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
    mp_limb_t **setptr;
    mp_limb_t **setinvertptr;
    mp_limb_t *setlastlimbptr;
    mp_limb_t *setlastlimbinvertptr;
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
    setptr = sfpartition;
    setinvertptr = sfpartitioninvert;
    (*setptr)[0] = (mp_limb_t)1;
    (*setinvertptr)[limballoc - 1] = (mp_limb_t)1<<(mp_bits_per_limb-1);
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
            setlastlimbptr = &((*setptr)[limbsize -1]);
            setlastlimbinvertptr = &((*setinvertptr)[limballoc - wlimbsize]);
            mpn_rshift(work1, setlastlimbinvertptr, wlimbsize, shift);
            mpn_and_n(work0, *setptr, work1, wlimbsize);
            isSumFree = mpn_zero_p(work0, wlimbsize);
            if (isSumFree) {
                // Ajouter n+1 à la huche i
                mask = (mp_limb_t)1<<nmodbpl;
                *setlastlimbptr |= mask;
                mask = (mp_limb_t)1<<(shift - 1);
                *setlastlimbinvertptr |= mask;
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
                setptr = sfpartition;
                setinvertptr = sfpartitioninvert;
            } else {
                i++;
                setptr++;
                setinvertptr++;
            }
        }
        if (i == p && p < pmax) {
            // Cas particulier où il faut remplir une nouvelle huche
            mpz_add_ui(iternum, iternum, 1);
            setlastlimbptr = &((*setptr)[limbsize -1]);
            setlastlimbinvertptr = &((*setinvertptr)[limballoc - wlimbsize]);
            mask = (mp_limb_t)1<<nmodbpl;
            *setlastlimbptr |= mask; // Ajoute n+1 à la huche i
            mask = (mp_limb_t)1<<(shift - 1);
            *setlastlimbinvertptr |= mask;
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
            setptr = sfpartition;
            setinvertptr = sfpartitioninvert;
            setlastlimbptr = &((*setptr)[limbsize -1]);
            setlastlimbinvertptr = &((*setinvertptr)[limballoc - wlimbsize]);
            p++;
            //printPartition(p, n, sfpartition);
            //printPartition(p, mp_bits_per_limb*limballoc, sfpartitioninvert);
        } else {
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
            setptr = sfpartition;
            while (i < p && !(*setptr[limbsize-1] & mask)) {
                i++;
                setptr++;
            }
            (*setptr)[limbsize-1] ^= mask;
            mask = (mp_limb_t)1<<(shift - 1);
            setinvertptr = &sfpartitioninvert[i];
            (*setinvertptr)[limballoc - limbsize] ^= mask;
            if (i == p-1) {
                if (mpn_zero_p(*setptr, limbsize)) {
                    p--;
                }
            }
            i++;
            setptr++;
            setinvertptr++;
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

unsigned long schurNumberIterative22(unsigned long pmax, unsigned long *nbests, mpz_t iternum) {
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

//int main(int argc, const char * argv[]) {
//    int i;
//    unsigned long schurNumbers[4];
//    unsigned long partitionNumbers[4];
//    partitionNumbers[0] = 0;
//    partitionNumbers[1] = 0;
//    partitionNumbers[2] = 0;
//    partitionNumbers[3] = 0;
//    schurNumberIterative(4, schurNumbers, partitionNumbers);
//    for (i=0; i<4; i++) {
//        printf("Schur Number of %i : %lu en nombre %lu \n", i+1, schurNumbers[i], partitionNumbers[i]);
//    }
//    return 0;
//}
