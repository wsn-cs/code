//
//  schur_number_unstack_memory.c
//  SchurNumber
//
//  Created by rubis on 07/03/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include <stdlib.h>
#include <gmp.h>

struct sftest_struc {
    mp_size_t limbsize;
    unsigned int shift;
};

typedef struct sftest_struc sftest_t;

void testinsert(sftest_t *sftests, unsigned long nblocking) {
    /*Insère le nouvel n bloquant dans le tableau sftests en respectant un ordre décroissant*/
    sftest_t *p;
    mp_size_t limbsize[2];
    unsigned int shift[2];
    char i;
    p = sftests;
    *limbsize = ((nblocking + 1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
    *shift = mp_bits_per_limb - nblocking%mp_bits_per_limb;
    while (p->limbsize > *limbsize) {
        p++;
    }
    limbsize[1] = p->limbsize;
    shift[1] = p->shift;
    p->limbsize = *limbsize;
    p->shift = *shift;
    i = 1;
    while (limbsize[i]) {
        p++;
        i ^= 1;
        limbsize[i] = p->limbsize;
        shift[i] = p->shift;
        p->limbsize = limbsize[i^1];
        p->shift = shift[i^1];
    }
}

char testperform(sftest_t *sftests, mp_limb_t *set, mp_limb_t *setinvert, mp_size_t limballoc, mp_limb_t *work0, mp_limb_t *work1) {
    /*Parcourt le tableau sftests et vérifie si il est possible d'ajouter n à set en conservant la sans-sommité*/
    sftest_t *p;
    mp_size_t wlimbsize;
    unsigned int shift;
    char isSumFree;
    p = sftests;
    isSumFree = 1;
    while (isSumFree && (wlimbsize = p->limbsize)) {
        shift = p->shift;
        mpn_rshift(work1, &setinvert[limballoc - wlimbsize], wlimbsize, shift);
        mpn_and_n(work0, set, work1, wlimbsize);
        isSumFree = mpn_zero_p(work0, wlimbsize);
        p++;
    }
    return isSumFree;
}

void testdelete(sftest_t *sftests, mp_limb_t *set, mp_limb_t *setinvert, mp_size_t limballoc, mp_limb_t *work0, mp_limb_t *work1) {
    /*Supprime du tableau les tests vérifiés*/
    sftest_t *p;
    sftest_t *deletep;
    mp_size_t wlimbsize;
    unsigned int shift;
    p = sftests;
    deletep = p;
    while ((wlimbsize = p->limbsize)) {
        shift = p->shift;
        mpn_rshift(work1, &setinvert[limballoc - wlimbsize], wlimbsize, shift);
        mpn_and_n(work0, set, work1, wlimbsize);
        if (!mpn_zero_p(work0, wlimbsize)) {
            if (p != deletep) {
                deletep->limbsize = wlimbsize;
                deletep->shift = shift;
            }
            deletep++;
        }
        p++;
    }
    while (deletep != p) {
        deletep->limbsize = 0;
        deletep->shift = 0;
        deletep++;
    }
}

void testclear(sftest_t *sftests) {
    sftest_t *p;
    p = sftests;
    while (p->limbsize) {
        p->limbsize = 0;
        p->shift = 0;
        p++;
    }
}

unsigned long schurNumberIterUnstackMemory(unsigned long pmax, unsigned long *nbests, mpz_t iternum) {
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
    unsigned long nblocking;
    unsigned long nblockingmax;
    unsigned long nbestlocal;
    unsigned long nmed;
    unsigned long i;
    unsigned long j;
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
    //mp_size_t wblimbsize;
    unsigned int nmodbpl;
    unsigned int shift;
    //unsigned int bshift;
    mp_limb_t mask;
    mp_limb_t mask2;
    sftest_t **memory;
    
    /*Initialisation*/
    limballoc = pmax;
    work0 = calloc(limballoc, sizeof(mp_limb_t));
    work1 = calloc(limballoc, sizeof(mp_limb_t));
    sfpartition = calloc(pmax, sizeof(mp_limb_t *));  //Tableau contenant la partition
    sfpartitioninvert = calloc(pmax, sizeof(mp_limb_t *));    //Tableau contenant les ensembles "inverses" de la partition
    memory = calloc(pmax, sizeof(sftest_t *));
    for (i=0; i<pmax; i++) {
        sfpartition[i] = calloc(limballoc, sizeof(mp_limb_t));
        sfpartitioninvert[i] = calloc(limballoc, sizeof(mp_limb_t));
        memory[i] = calloc(mp_bits_per_limb*pmax, sizeof(sftest_t));
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
                // Vérifier si la huche i peut toujours contenir nbest
                isSumFree = testperform(memory[i], sfpartition[i], sfpartitioninvert[i], limballoc, work0, work1);
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
                    //printPartition1(p, GMP_NUMB_BITS * limbsize, sfpartition);
                    //printPartition1(p, GMP_NUMB_BITS * limballoc, sfpartitioninvert);
                } else {
                    i++;
                }
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
            for (i=0; i<p; i++) {
                testclear(memory[i]);
            }
            i = 0;
            isAppendable = 0;
            p++;
            //printPartition1(p, GMP_NUMB_BITS * limbsize, sfpartition);
            //printPartition1(p, GMP_NUMB_BITS * limballoc, sfpartitioninvert);
        } else if(isAppendable) {
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
            // Vérifier si la huche i peut désormais contenir nbest
            for (j=0; j<p; j++) {
                if (j != i) {
                    testdelete(memory[j], sfpartition[i], sfpartitioninvert[i], limballoc, work0, work1);
                }
            }
            if (i == p-1) {
                if (mpn_zero_p(sfpartition[i], limbsize)) {
                    p--;
                }
            }
            i++;
            isAppendable = 1;
        } else {
            // Dépiler jusqu'à revenir à nblockingmax
            nbestlocal = n;
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
                if (i < p) {
                    for (j=0; j<p; j++) {
                        if (j != i) {
                            testdelete(memory[j], sfpartition[i], sfpartitioninvert[i], limballoc, work0, work1);
                        }
                    }
                }
            }
            testinsert(memory[iblocking], nbestlocal);
            while (mpn_zero_p(sfpartition[p-1], limbsize)) {
                p--;
            }
            i = iblocking + 1;
            isAppendable = 1;
            //printPartition1(p, GMP_NUMB_BITS * limbsize, sfpartition);
            //printPartition1(p, GMP_NUMB_BITS * limballoc, sfpartitioninvert);
        }
    }
    
    /*Nettoyage*/
    for (i=0; i<pmax; i++) {
        free(sfpartition[i]);
        free(sfpartitioninvert[i]);
        free(memory[i]);
    }
    free(sfpartition);
    free(sfpartitioninvert);
    free(memory);
    free(work0);
    free(work1);
    return nbests[pmax-1];
}
