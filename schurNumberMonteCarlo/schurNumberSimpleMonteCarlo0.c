//
//  schurNumberSimpleMonteCarlo.c
//  schurNumberMonteCarlo
//
//  Created by rubis on 21/03/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include <stdlib.h>
#include <gmp.h>
#include <stdio.h>

#if GMP_NUMB_BITS == 64
#define GMP_2EXP_NUMB_BITS 6
#elif GMP_NUMB_BITS == 128
#define GMP_2EXP_NUMB_BITS 7
#elif GMP_NUMB_BITS == 32
#define GMP_2EXP_NUMB_BITS 5
#endif

#define REALLOC(ptr, reptr, size, resize) \
do { \
reptr = realloc(ptr, resize);\
while (reptr == NULL && resize > size) { \
resize -= resize >> 1;\
 reptr = realloc(ptr, resize);\
} \
ptr = reptr;\
} while(0)

void printPartition(unsigned long p, unsigned long n, mp_limb_t **partition) {
    /*Affiche une partition.*/
    unsigned long limbn;
    mp_size_t limbsize;
    unsigned long i;
    mp_bitcnt_t j;
    mp_limb_t *set;
    mp_limb_t limb;
    limbsize = (n>>6) + 1;
    printf("Partition:\n");
    for (i=0; i<p; i++) {
        printf("\t");
        set = partition[i];
        for (limbn=0; limbn<limbsize; limbn++) {
            limb = set[limbn];
            for (j=0; j<mp_bits_per_limb; j++) {
                if (limb & ((unsigned long)1<<j)) {
                    printf(" %lu", limbn*mp_bits_per_limb + j + 1);
                }
            }
        }
        printf("\n");
    }
}

unsigned long schurNumberSimpleNestedMonteCarlo1(unsigned int pmax) {
    unsigned long nglobal, ntest;
    unsigned long nbest, nbestlocal;
    unsigned long i, itest, index, indexbestlocal;
    unsigned long pglobal, prand, ptest;
    unsigned long *ivalid;
    unsigned long inum;
    char isSumFree, isNotFinished;
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t *set;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    mp_limb_t **sfpartitionbest;
    mp_size_t limballoc;    // Nombre de limbes alloué à chaque ensemble de sfpartition et à work0 et work1
    mp_size_t limbsizeglobal, limbsizetest;     // Nombre de limbes utilisés par les ensembles de sfpartition
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
    sfpartitionbest = calloc(pmax, sizeof(mp_limb_t *));  //Tableau contenant la meilleure partition trouvée
    for (i=0; i<pmax; i++) {
        sfpartition[i] = calloc(limballoc, sizeof(mp_limb_t));
        sfpartitioninvert[i] = calloc(limballoc, sizeof(mp_limb_t));
        sfpartitionbest[i] = calloc(limballoc, sizeof(mp_limb_t));
    }
    ivalid = calloc(1 + (pmax>>GMP_2EXP_NUMB_BITS), sizeof(unsigned long));  //Tableau contenant les indices des huches prolongeables
    
    /*Création de la première partition {{1}}*/
    nglobal = 1;  // Taille de l'intervalle
    nbest = 1; // Meilleur n trouvé
    i = 0;  // Huche où placer l'entier suivant
    pglobal = 1;  // Nombre de huches non vides
    prand = 2; // Nombre de huches à sélectionner
    *sfpartition[0] = (mp_limb_t)1;
    sfpartitioninvert[0][limballoc - 1] = (mp_limb_t)1<<(mp_bits_per_limb-1);
    limbsizeglobal = 1;
    wlimbsize = 1;    // Taille en limbe de la seconde moitié
    nmodbpl = 1;
    shift = mp_bits_per_limb - nmodbpl;
    isNotFinished = 1;
    
    while (1) {
         // Trouver les huches pouvant accepter n+1
        inum = 0;
        for (i=0; i<pglobal; i++) {
            // Tester si l'ensemble obtenu en ajoutant n+1 à la huche i est sans-somme
            set = sfpartitioninvert[i];
            mpn_copyd(work0, &set[limballoc - wlimbsize], wlimbsize);
            mpn_rshift(work1, work0, wlimbsize, shift);
            set = sfpartition[i];
            mpn_and_n(work0, set, work1, wlimbsize);
            isSumFree = mpn_zero_p(work0, wlimbsize);
            if (isSumFree) {
                // Ajouter l'indice au tableau des indices valides
                ivalid[inum] = i;
                inum++;
            }
        }
        if (pglobal < pmax) {
            // Quand il est possible de remplir une nouvelle huche
            ivalid[inum] = pglobal;
            inum++;
        }
        if (!inum) {
            // Aucune huche n'est prolongeable
            break;
        }
        
        
        // Initialisation de la descente aléatoire
        nbestlocal = nglobal;
        indexbestlocal = 0;
        ntest = nglobal;
        nmodbpl = ntest%mp_bits_per_limb;
        shift = mp_bits_per_limb - nmodbpl;
        
        if (inum != 1) {
            for (index=0; index<inum; index++) {
                // Descente aléatoire de niveau 0
                itest = ivalid[index];
                // Placer n+1 dans la huche ivalid[index]
                mask = (mp_limb_t)1<<nmodbpl;
                sfpartition[itest][limbsizeglobal - 1] |= mask;
                mask = (mp_limb_t)1<<(shift - 1);
                sfpartitioninvert[itest][limballoc - limbsizeglobal] |= mask;
                ntest++;
                ptest = pglobal;
                if (itest == pglobal) {
                    ptest++;
                }
                
                // Définition de prand
                if (ptest == pmax) {
                    prand = ptest;
                } else {
                    prand = ptest+1;
                }
                limbsizetest = limbsizeglobal;
                isSumFree = 1;
                while (isSumFree) {
                    // Sélectionner une huche aléatoirement
                    itest = arc4random_uniform(prand);
                    nmodbpl = ntest%mp_bits_per_limb;
                    shift = mp_bits_per_limb - nmodbpl;
                    wlimbsize = ((ntest+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                    if (itest < ptest) {
                        // Quand la huche n'est pas vide, tester si l'ensemble obtenu en ajoutant ntest+1 à la huche itest est sans-somme
                        set = sfpartitioninvert[itest];
                        mpn_copyd(work0, &set[limballoc - wlimbsize], wlimbsize);
                        mpn_rshift(work1, work0, wlimbsize, shift);
                        set = sfpartition[itest];
                        mpn_and_n(work0, set, work1, wlimbsize);
                        isSumFree = mpn_zero_p(work0, wlimbsize);
                    } else {
                        // Placer ntest+1 dans une huche vide
                        isSumFree = 1;
                        ptest++;
                        if (ptest < pmax) {
                            prand++;
                        }
                    }
                    if (isSumFree) {
                        // Ajouter ntest+1 à la huche itest
                        mask = (mp_limb_t)1<<nmodbpl;
                        sfpartition[itest][limbsizetest -1] |= mask;
                        mask = (mp_limb_t)1<<(shift - 1);
                        sfpartitioninvert[itest][limballoc - limbsizetest] |= mask;
                        ntest++;
                        if (ntest > nbest) {
                            //nbest = ntest;
                            if (ntest >= mp_bits_per_limb*limballoc) {
                                limballoc *= 2;
                                work0 = realloc(work0, limballoc);
                                work1 = realloc(work1, limballoc);
                                for (itest=0; itest<pmax; itest++) {
                                    sfpartition[itest] = realloc(sfpartition[itest], limballoc);
                                    sfpartitioninvert[itest] = realloc(sfpartitioninvert[itest], limballoc);
                                    mpn_lshift(sfpartitioninvert[itest], sfpartitioninvert[itest], limballoc, limballoc>>1);
                                    sfpartitionbest[itest] = realloc(sfpartitionbest[itest], limballoc);
                                }
                            }
                        }
                        if (!nmodbpl) {
                            for (itest=0; itest<pmax; itest++) {
                                sfpartition[itest][limbsizetest] = (mp_limb_t)0;
                                sfpartitioninvert[itest][limballoc - limbsizetest -1] = (mp_limb_t)0;
                            }
                            limbsizetest++;
                        }
                    }
                }
                ntest--;
                if (ntest > nbestlocal) {
                    indexbestlocal = index;
                    nbestlocal = ntest;
                }
                if (ntest > nbest) {
                    // Copier sfpartition dans sfpartitionbest
                    nbest = ntest;
                    for (itest=0; itest<ptest; itest++) {
                        mpn_copyi(sfpartitionbest[itest], sfpartition[itest], limbsizetest);
                    }
                }
                
                // Vider les huches nouvellement remplies
                for (itest=pglobal; itest<ptest; itest++) {
                    mpn_zero_p(sfpartition[itest], limbsizetest);
                    mpn_zero_p(sfpartitioninvert[itest], limbsizetest);
                }
                // Nettoyer les autres et réinitialiser ntest
                ntest = nglobal;
                nmodbpl = ntest%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                mask = (mp_limb_t)GMP_NUMB_MAX>>(shift);
                for (itest=0; itest<pglobal; itest++) {
                    sfpartition[itest][limbsizeglobal -1] &= mask;
                    //mpn_zero_p(sfpartition[i] + limbsize, limbsizetest - limbsize);
                }
                mask = (mp_limb_t)GMP_NUMB_MAX<<(nmodbpl);
                for (itest=0; itest<pglobal; itest++) {
                    sfpartitioninvert[itest][limballoc - limbsizeglobal] &= mask;
                    //mpn_zero_p(sfpartitioninvert[i] + limballoc - limbsizetest - 1, limbsizetest - limbsize);
                }
            }
        }
        
        /*if (nbestlocal == nglobal) {
            // Aucune partition plus grande n'a été trouvée
            break;
        }*/
        
        // Placer n+1 dans la huche indexbestlocal
        nmodbpl = nglobal%mp_bits_per_limb;
        shift = mp_bits_per_limb - nmodbpl;
        //mask = (mp_limb_t)GMP_NUMB_MAX>>(shift);
        //sfpartition[index][limbsizeglobal -1] &= mask;
        mask = (mp_limb_t)1<<nmodbpl;
        i = ivalid[indexbestlocal];
        sfpartition[i][limbsizeglobal - 1] |= mask;
        //mask = (mp_limb_t)GMP_NUMB_MAX<<(nmodbpl);
        //sfpartitioninvert[index][limballoc - limbsizeglobal] &= mask;
        mask = (mp_limb_t)1<<(shift - 1);
        sfpartitioninvert[i][limballoc - limbsizeglobal] |= mask;
        nglobal++;
        nmodbpl = nglobal%mp_bits_per_limb;
        shift = mp_bits_per_limb - nmodbpl;
        wlimbsize = ((nglobal+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
        if (i == pglobal) {
            pglobal++;
        }
        if (!nmodbpl) {
            limbsizeglobal++;
        }
    }
    
    /*Nettoyage*/
    for (i=0; i<pmax; i++) {
        free(sfpartition[i]);
        free(sfpartitioninvert[i]);
        free(sfpartitionbest[i]);
    }
    free(sfpartition);
    free(sfpartitioninvert);
    free(sfpartitionbest);
    free(work0);
    free(work1);
    free(ivalid);
    
    return nbest;
}

unsigned long schurNumberSimpleNestedMonteCarlo(unsigned int pmax, unsigned int nestinglevel, mpz_t iternum) {
    unsigned long nglobal, nnest, ntest;
    unsigned long nbest, nbestlocal;
    unsigned long i, itest, index, indexbestlocal;
    unsigned long pglobal, pnest, prand, ptest;
    //unsigned long *ivalid;
    //unsigned long inum;
    unsigned long l;
    char isSumFree, isNotFinished;
    mp_limb_t *work0, *work1, *reallocptr;
    mp_limb_t *set;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    mp_limb_t **sfpartitionbest;
    mp_size_t limballoc, limbrealloc;    // Nombre de limbes alloué à chaque ensemble de sfpartition et à work0 et work1
    mp_size_t limbsizeglobal, limbsizenest, limbsizetest;     // Nombre de limbes utilisés par les ensembles de sfpartition
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
    sfpartitionbest = calloc(pmax, sizeof(mp_limb_t *));  //Tableau contenant la meilleure partition trouvée
    for (i=0; i<pmax; i++) {
        sfpartition[i] = calloc(limballoc, sizeof(mp_limb_t));
        sfpartitioninvert[i] = calloc(limballoc, sizeof(mp_limb_t));
        sfpartitionbest[i] = calloc(limballoc, sizeof(mp_limb_t));
    }
    //ivalid = calloc(pmax * nestinglevel, sizeof(unsigned long));  //Tableau contenant les indices des huches prolongeables
    
    /*Création de la première partition {{1}}*/
    nglobal = 1;  // Taille de l'intervalle
    nbest = 1; // Meilleur n trouvé
    i = 0;  // Huche où placer l'entier suivant
    pglobal = 1;  // Nombre de huches non vides
    prand = 1; // Nombre de huches à sélectionner
    l = nestinglevel; //Niveau d'imbrication actuelle
    *sfpartition[0] = (mp_limb_t)1;
    sfpartitioninvert[0][limballoc - 1] = (mp_limb_t)1<<(mp_bits_per_limb-1);
    limbsizeglobal = 1;
    wlimbsize = 1;    // Taille en limbe de la seconde moitié
    nmodbpl = 1;
    shift = mp_bits_per_limb - nmodbpl;
    isNotFinished = 1;
    
    while (1) {
        /*
        // Trouver les huches pouvant accepter n+1
        inum = 0;
        for (i=0; i<p; i++) {
            // Tester si l'ensemble obtenu en ajoutant n à la huche i est sans-somme
            mpz_add_ui(iternum, iternum, 1);
            set = sfpartitioninvert[i];
            mpn_copyd(work0, &set[limballoc - wlimbsize], wlimbsize);
            mpn_rshift(work1, work0, wlimbsize, shift);
            set = sfpartition[i];
            mpn_and_n(work0, set, work1, wlimbsize);
            isSumFree = mpn_zero_p(work0, wlimbsize);
            if (isSumFree) {
                // Ajouter l'indice au tableau des indices valides
                ivalid[inum] = i;
                inum++;
            }
        }
        if (p < pmax) {
            // Quand il est possible de remplir une nouvelle huche
            ivalid[inum] = p;
            inum++;
        }
        if (!inum) {
            // Aucune huche n'est prolongeable
            break;
        }*/
        
        nbestlocal = nglobal;
        indexbestlocal = 0;
        
        /*Trouver les suites de placement de longueur nestinglevel valides*/
        l = 1;
        nnest = nglobal;
        pnest = pglobal;
        limbsizenest = limbsizeglobal;
        i = 0;
        index = 0;
        
        while (l > 0 && l <= nestinglevel) {
            if (i < pnest) {
                // Tester si l'ensemble obtenu en ajoutant nnest = n+l à la huche i est sans-somme
                mpz_add_ui(iternum, iternum, 1);
                set = sfpartitioninvert[i];
                mpn_copyd(work0, &set[limballoc - wlimbsize], wlimbsize);
                mpn_rshift(work1, work0, wlimbsize, shift);
                set = sfpartition[i];
                mpn_and_n(work0, set, work1, wlimbsize);
                isSumFree = mpn_zero_p(work0, wlimbsize);
                if (isSumFree) {
                    // Ajouter nnest à la huche i
                    mask = (mp_limb_t)1<<nmodbpl;
                    sfpartition[i][limbsizenest -1] |= mask;
                    mask = (mp_limb_t)1<<(shift - 1);
                    sfpartitioninvert[i][limballoc - limbsizenest] |= mask;
                    if (l == 1) {
                        index = i;
                    }
                    l++;
                    nnest++;
                    nmodbpl = (nnest)%mp_bits_per_limb;
                    shift = mp_bits_per_limb - nmodbpl;
                    wlimbsize = ((nnest+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                    if (nnest > nbest) {
                        if (nnest >= mp_bits_per_limb*limballoc) {
                            limbrealloc = 2*limballoc;
                            /*reallocptr = realloc(work0, limbrealloc);
                            while (reallocptr == NULL && limbrealloc > limballoc) {
                                limbrealloc =- limbrealloc>>1;
                                reallocptr = realloc(work0, limbrealloc);
                            }
                            work0 = reallocptr;*/
                            REALLOC(work0, reallocptr, limballoc, limbrealloc);
                            REALLOC(work1, reallocptr, limballoc, limbrealloc);
                            for (i=0; i<pmax; i++) {
                                REALLOC(sfpartition[i], reallocptr, limballoc, limbrealloc);
                                REALLOC(sfpartitionbest[i], reallocptr, limballoc, limbrealloc);
                                REALLOC(sfpartitioninvert[i], reallocptr, limballoc, limbrealloc);
                                mpn_lshift(sfpartitioninvert[i], sfpartitioninvert[i], limbrealloc, limbrealloc - limballoc);
                            }
                            limballoc = limbrealloc;
                        }
                    }
                    if (!nmodbpl) {
                        for (i=0; i<pmax; i++) {
                            sfpartition[i][limbsizenest] = (mp_limb_t)0;
                            sfpartitioninvert[i][limballoc - limbsizenest -1] = (mp_limb_t)0;
                        }
                        limbsizenest++;
                    }
                    i = 0;
                } else {
                    i++;
                }
            } else if (i == pnest && pnest < pmax) {
                // Cas particulier où il faut remplir une nouvelle huche
                mpz_add_ui(iternum, iternum, 1);
                mask = (mp_limb_t)1<<nmodbpl;
                sfpartition[i][limbsizenest -1] |= mask; // Ajoute nnest à la huche i
                mask = (mp_limb_t)1<<(shift - 1);
                sfpartitioninvert[i][limballoc - limbsizenest] |= mask;
                if (l == 1) {
                    index = i;
                }
                l++;
                nnest++;
                nmodbpl = nnest%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                wlimbsize = ((nnest+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                if (nnest > nbest) {
                    if (nnest >= 64*limballoc -1) {
                        limbrealloc = 2*limballoc;
                        REALLOC(work0, reallocptr, limballoc, limbrealloc);
                        REALLOC(work1, reallocptr, limballoc, limbrealloc);
                        for (i=0; i<pmax; i++) {
                            REALLOC(sfpartition[i], reallocptr, limballoc, limbrealloc);
                            REALLOC(sfpartitionbest[i], reallocptr, limballoc, limbrealloc);
                            REALLOC(sfpartitioninvert[i], reallocptr, limballoc, limbrealloc);
                            mpn_lshift(sfpartitioninvert[i], sfpartitioninvert[i], limbrealloc, limbrealloc - limballoc);
                        }
                        limballoc = limbrealloc;
                    }
                }
                if (!nmodbpl) {
                    for (i=0; i<pmax; i++) {
                        sfpartition[i][limbsizenest] = (mp_limb_t)0;
                        sfpartitioninvert[i][limballoc - limbsizenest -1] = (mp_limb_t)0;
                    }
                    limbsizenest++;
                }
                i = 0;
                pnest++;
            } else {
                // Dépiler
                if (!nmodbpl) {
                    limbsizenest--;
                }
                l--;
                if (l) {
                    nnest--;
                    nmodbpl = nnest%mp_bits_per_limb;
                    shift = mp_bits_per_limb - nmodbpl;
                    wlimbsize = ((nnest+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                    mask = (mp_limb_t)1<<nmodbpl;
                    i = 0;
                    while (i < pnest && !(sfpartition[i][limbsizenest-1] & mask)) {
                        i++;
                    }
                    sfpartition[i][limbsizenest-1] ^= mask;
                    mask = (mp_limb_t)1<<(shift - 1);
                    sfpartitioninvert[i][limballoc - limbsizenest] ^= mask;
                    if (i == pnest-1) {
                        if (mpn_zero_p(sfpartition[i], limbsizenest)) {
                            pnest--;
                        }
                    }
                    i++;
                }
            }
            
            if (l == nestinglevel) {
                /*Réaliser la descente aléatoire de niveau 0*/
                indexbestlocal = 0;
                ntest = nnest;
                nmodbpl = ntest%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                ptest = pnest;
                
                // Définition de prand
                if (ptest == pmax) {
                    prand = ptest;
                } else {
                    prand = ptest + 1;
                }
                limbsizetest = limbsizenest;
                isSumFree = 1;
                while (isSumFree) {
                    // Sélectionner une huche aléatoirement
                    itest = arc4random_uniform(prand);
                    nmodbpl = ntest%mp_bits_per_limb;
                    shift = mp_bits_per_limb - nmodbpl;
                    wlimbsize = ((ntest+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                    if (itest < ptest) {
                        // Quand la huche n'est pas vide, tester si l'ensemble obtenu en ajoutant ntest+1 à la huche i est sans-somme
                        mpz_add_ui(iternum, iternum, 1);
                        set = sfpartitioninvert[itest];
                        mpn_copyd(work0, &set[limballoc - wlimbsize], wlimbsize);
                        mpn_rshift(work1, work0, wlimbsize, shift);
                        set = sfpartition[itest];
                        mpn_and_n(work0, set, work1, wlimbsize);
                        isSumFree = mpn_zero_p(work0, wlimbsize);
                    } else {
                        // Placer ntest+1 dans une huche vide
                        isSumFree = 1;
                        ptest++;
                        if (ptest < pmax) {
                            prand++;
                        }
                    }
                    if (isSumFree) {
                        // Ajouter ntest+1 à la huche i
                        mask = (mp_limb_t)1<<nmodbpl;
                        sfpartition[itest][limbsizetest -1] |= mask;
                        mask = (mp_limb_t)1<<(shift - 1);
                        sfpartitioninvert[itest][limballoc - limbsizetest] |= mask;
                        ntest++;
                        if (ntest > nbest) {
                            nbest = ntest;
                            if (ntest >= mp_bits_per_limb*limballoc) {
                                limbrealloc = 2*limballoc;
                                REALLOC(work0, reallocptr, limballoc, limbrealloc);
                                REALLOC(work1, reallocptr, limballoc, limbrealloc);
                                for (itest=0; itest<pmax; itest++) {
                                    REALLOC(sfpartition[itest], reallocptr, limballoc, limbrealloc);
                                    REALLOC(sfpartitionbest[itest], reallocptr, limballoc, limbrealloc);
                                    REALLOC(sfpartitioninvert[itest], reallocptr, limballoc, limbrealloc);
                                    mpn_lshift(sfpartitioninvert[itest], sfpartitioninvert[i], limbrealloc, limbrealloc - limballoc);
                                }
                                limballoc = limbrealloc;
                            }
                        }
                        if (!nmodbpl) {
                            for (itest=0; itest<pmax; itest++) {
                                sfpartition[itest][limbsizetest] = (mp_limb_t)0;
                                sfpartitioninvert[itest][limballoc - limbsizetest -1] = (mp_limb_t)0;
                            }
                            limbsizetest++;
                        }
                    }
                }
                if (ntest > nbestlocal) {
                    indexbestlocal = index;
                    nbestlocal = ntest;
                }
                if (ntest > nbest) {
                    // Copier sfpartition dans sfpartitionbest
                    for (itest=0; itest<ptest; itest++) {
                        mpn_copyi(sfpartitionbest[itest], sfpartition[itest], limbsizetest);
                    }
                    ntest = nbest;
                }
                
                // Vider les huches nouvellement remplies
                for (itest=pnest; itest<ptest; itest++) {
                    mpn_zero_p(sfpartition[itest], limbsizetest);
                    mpn_zero_p(sfpartitioninvert[itest], limbsizetest);
                }
                // Nettoyer les autres et réinitialiser ntest
                ntest = nnest;
                nmodbpl = ntest%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                mask = (mp_limb_t)GMP_NUMB_MAX>>(shift);
                //printPartition(pnest, ntest, sfpartition);
                for (itest=0; itest<pnest; itest++) {
                    sfpartition[itest][limbsizenest -1] &= mask;
                    //mpn_zero_p(sfpartition[i] + limbsize, limbsizetest - limbsize);
                }
                //printPartition(pnest, ntest, sfpartition);
                //mask = (mp_limb_t)GMP_NUMB_MAX<<(nmodbpl);
                mask = (mp_limb_t)GMP_NUMB_MAX<<shift;
                for (itest=0; itest<pnest; itest++) {
                    sfpartitioninvert[itest][limballoc - limbsizenest] &= mask;
                    //mpn_zero_p(sfpartitioninvert[i] + limballoc - limbsizetest - 1, limbsizetest - limbsize);
                }
                
                /*for (index=0; index<inum; index++) {
                    // Descente aléatoire de niveau 0
                    itest = ivalid[index];
                    // Placer n+1 dans la huche ivalid[index]
                    mask = (mp_limb_t)1<<nmodbpl;
                    sfpartition[itest][limbsize - 1] |= mask;
                    mask = (mp_limb_t)1<<(shift - 1);
                    sfpartitioninvert[i][limballoc - limbsize] |= mask;
                    ntest++;
                    ptest = p;
                }*/
                
                // Dépiler
                if (!(nnest%mp_bits_per_limb)) {
                    limbsizenest--;
                }
                l--;
                nnest--;
                nmodbpl = nnest%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                wlimbsize = ((nnest+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                mask = (mp_limb_t)1<<nmodbpl;
                i = 0;
                while (i < pnest && !(sfpartition[i][limbsizenest-1] & mask)) {
                    i++;
                }
                sfpartition[i][limbsizenest-1] ^= mask;
                mask = (mp_limb_t)1<<(shift - 1);
                sfpartitioninvert[i][limballoc - limbsizenest] ^= mask;
                if (i == pnest-1) {
                    if (mpn_zero_p(sfpartition[i], limbsizenest)) {
                        pnest--;
                    }
                }
                //printPartition(pnest, nnest, sfpartition);
                i++;
            }
            
        }
        
        if (nbestlocal == nglobal) {
            // Aucune partition plus grande n'a été trouvée
            break;
        }
        
        // Placer n+1 dans la huche indexbestlocal
        nmodbpl = nglobal%mp_bits_per_limb;
        shift = mp_bits_per_limb - nmodbpl;
        //mask = (mp_limb_t)GMP_NUMB_MAX>>(shift);
        //sfpartition[index][limbsizeglobal -1] &= mask;
        mask = (mp_limb_t)1<<nmodbpl;
        sfpartition[indexbestlocal][limbsizeglobal - 1] |= mask;
        //mask = (mp_limb_t)GMP_NUMB_MAX<<(nmodbpl);
        //sfpartitioninvert[index][limballoc - limbsizeglobal] &= mask;
        mask = (mp_limb_t)1<<(shift - 1);
        sfpartitioninvert[indexbestlocal][limballoc - limbsizeglobal] |= mask;
        nglobal++;
        if (index == pglobal) {
            pglobal++;
        }
        if (!nmodbpl) {
            limbsizeglobal++;
        }
        nmodbpl = nglobal%mp_bits_per_limb;
        shift = mp_bits_per_limb - nmodbpl;
        wlimbsize = ((nglobal+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
    }
    
    /*Nettoyage*/
    for (i=0; i<pmax; i++) {
        free(sfpartition[i]);
        free(sfpartitioninvert[i]);
        free(sfpartitionbest[i]);
    }
    free(sfpartition);
    free(sfpartitioninvert);
    free(sfpartitionbest);
    free(work0);
    free(work1);
    //free(ivalid);
    
    return nbest;
}
