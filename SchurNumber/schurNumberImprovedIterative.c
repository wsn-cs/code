//
//  schurNumberImprovedIterative.c
//  SchurNumber
//
//  Created by Gabriel Merlin on 28/02/2019.
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

unsigned long schurNumberImprovedIterative0(unsigned long pmax, unsigned long *nbests) {
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
    unsigned long nmod64;
    unsigned long i;
    unsigned long j;
    unsigned long p;
    unsigned long limb0;
    unsigned long limbsize;
    unsigned long wlimbsize;
    unsigned int shift1;
    unsigned int shift2;
    unsigned long limbmax;
    unsigned long nblocking;
    unsigned long iblocking;
    unsigned long n2unstack;
    char isSumFree;
    char joinable;
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t *work2;
    mp_limb_t *set;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    
    /*Initialisation*/
    limbmax = pmax; //Taille de l'entier work0 et work1
    work0 = calloc(limbmax, sizeof(mp_limb_t));
    work1 = calloc(limbmax, sizeof(mp_limb_t));
    work2 = calloc(limbmax, sizeof(mp_limb_t));
    sfpartition = calloc(pmax, sizeof(mpz_t));  //Tableau contenant la partition
    sfpartitioninvert = calloc(pmax, sizeof(mpz_t));    //Tableau contenant les ensembles "inverses" de la partition
    for (i=0; i<pmax; i++) {
        sfpartition[i] = calloc(pmax, sizeof(mp_limb_t));
        sfpartitioninvert[i] = calloc(pmax, sizeof(mp_limb_t));
        nbests[i] = i+1;
    }
    
    /*Création de la première partition {{1}}*/
    n = 1;  // Taille de l'intervalle
    i = 0;  // Huche où placer l'entier suivant
    p = 1;  // Nombre de huches non vides
    shift1 = 0;
    limb0 = 0;
    iblocking = 0;
    limbsize = 1; //Nombre de limbe d'un ensemble de la partition
    joinable = 0;
    *sfpartition[0] = 1;
    //printPartition(p, n, sfpartition);
    sfpartitioninvert[0][limbmax - 1] = (unsigned long)1<<(mp_bits_per_limb-1);
    
    /*Itération jusqu'à énumérer toutes les partions sum-free à au plus pmax huches*/
    while (n>0) {
        // Placer n+1 dans une des huches en conservant la sans-sommité
        //limbsize = (n>>mp_bits_per_limb) +1;
        limb0 = ((n+1)>>7); //Indice du limbe contenant n/2
        wlimbsize = limbsize - limb0;   //Taille du demi-entier
        shift1 = ((n>>1)%mp_bits_per_limb);
        //shift2 = mp_bits_per_limb - (((n+1)>>1)%mp_bits_per_limb);
        shift2 = mp_bits_per_limb - n + shift1;
        n2unstack = 0;
        iblocking = p;
        while (i < p) {
            // Tester si l'ensemble obtenu en ajoutant n à la huche i est sans-somme
            set = sfpartition[i];
            mpn_copyd(work0, &set[limb0], wlimbsize);
            mpn_rshift(work2, work0, wlimbsize, shift1);
            //mpn_rshift(work2, work0, wlimbsize, ((unsigned long)shift1>>16));
            set = sfpartitioninvert[i];
            mpn_copyd(work0, &set[limbmax - wlimbsize], wlimbsize);
            mpn_rshift(work1, work0, wlimbsize, shift2);
            //mpn_rshift(work1, work0, wlimbsize, ((unsigned long)shift2>>16));
            mpn_and_n(work0, work1, work2, wlimbsize);
            isSumFree = mpn_zero_p(work0, wlimbsize);
            if (isSumFree) {
                // Ajoute n+1 à la huche i
                nmod64 = n%64;
                if (!nmod64) {
                    //Commencer un nouveau limbe
                    limbsize++;
                    for (j=0; j<pmax; j++) {
                        sfpartition[i][limbsize -1] = 0;
                        sfpartitioninvert[i][limbmax - limbsize] = 0;
                    }
                    sfpartition[i][limbsize-1] = (unsigned long)1;
                    sfpartitioninvert[i][limbmax - limbsize] = (unsigned long)1<<(mp_bits_per_limb-1);
                } else {
                    sfpartition[i][limbsize-1] |= (unsigned long)1<<nmod64;
                    sfpartitioninvert[i][limbmax - limbsize] |= (unsigned long)1<<(64 - nmod64 -1);
                }
                n++;
                limb0 = ((n+1)>>7);
                wlimbsize = limbsize - limb0;
                shift1 = ((n>>1)%mp_bits_per_limb);
                //shift2 = mp_bits_per_limb - (((n+1)>>1)%mp_bits_per_limb);
                shift2 = mp_bits_per_limb - n + shift1;
                n2unstack = 0;
                iblocking = p;
                joinable = 0;
                //shift = ((n>>1)%mp_bits_per_limb);
                if (n > nbests[p-1]) {
                    nbests[p-1] = n;
                    if (n >= mp_bits_per_limb*limbmax -1) {
                        limbmax *= 2;
                        work0 = realloc(work0, limbmax);
                        work1 = realloc(work1, limbmax);
                        work2 = realloc(work2, limbmax);
                    }
                }
                i = 0;
                //printPartition(p, limbsize, sfpartition);
            } else {
                //nblocking = ((n+1)>>1) + mpn_scan1(work0, 0) -1;
                nblocking = mp_bits_per_limb*limb0 + shift1 + mpn_scan1(work0, 0);
                //printf("%lu %lu %lu \n", *work0, nblocking, n2unstack);
                if (nblocking > n2unstack) {
                    n2unstack = nblocking;
                    iblocking = i;
                }
                i++;
            }
        }
        if (i == p && p < pmax) {
            // Cas particulier où il faut remplir une nouvelle huche
            nmod64 = n%64;
            if (!nmod64) {
                // Commencer un nouveau limbe
                limbsize++;
                for (j=0; j<pmax; j++) {
                    sfpartition[i][limbsize -1] = (unsigned long)0;
                    sfpartitioninvert[i][limbmax - limbsize] = (unsigned long)0;
                }
                sfpartition[i][limbsize-1] = (unsigned long)1;
                sfpartitioninvert[i][limbmax - limbsize] = (unsigned long)1<<(mp_bits_per_limb-1);
            } else {
                sfpartition[i][limbsize - 1] |= (unsigned long)1<<nmod64;
                sfpartitioninvert[i][limbmax - limbsize] |= (unsigned long)1<<(64 - nmod64 -1);
            }
            n++;
            if (n > nbests[p]) {
                nbests[p] = n;
                if (n >= mp_bits_per_limb*limbmax -1) {
                    limbmax *= 2;
                    work0 = realloc(work0, limbmax);
                    work1 = realloc(work1, limbmax);
                    work2 = realloc(work2, limbmax);
                }
            }
            i = 0;
            joinable = 0;
            p++;
            //printPartition(p, limbsize, sfpartition);
        } else if (joinable) {
            // Dépiler
            n2unstack = n-1;
            nmod64 = (n2unstack)%64;
            limb0 = n2unstack>>6;
            i = 0;
            while (i < p && !(sfpartition[i][limb0] & ((unsigned long)1<<nmod64))) {
                i++;
            }
            iblocking = i;
            for (j=0; j<pmax; j++) {
                sfpartition[j][limb0] &= (unsigned long)0xffffffffffffffff>>(64 - nmod64);
                sfpartitioninvert[j][limbmax - limb0 - 1] &= (unsigned long)0xffffffffffffffff<<(64 - nmod64);
            }
            limbsize = limb0 + 1;
            while (mpn_zero_p(sfpartition[p-1], limbsize)) {
                p--;
            }
            n = n2unstack;
            i = iblocking + 1;
            joinable = 1;
            //printf("A dépiler: %lu\n", n+1);
            //printf("Nouvelle ");
            //printPartition(p, n, sfpartition);
            //printf("Inverse ");
            //printPartition(p, pmax*64, sfpartitioninvert);
        } else {
            // Dépiler
            nmod64 = (n2unstack)%64;
            limb0 = n2unstack>>6;
            for (j=0; j<pmax; j++) {
                sfpartition[j][limb0] &= (unsigned long)0xffffffffffffffff>>(64 - nmod64);
                sfpartitioninvert[j][limbmax - limb0 - 1] &= (unsigned long)0xffffffffffffffff<<(64 - nmod64);
            }
            limbsize = limb0 + 1;
            while (mpn_zero_p(sfpartition[p-1], limbsize)) {
                p--;
            }
            n = n2unstack;
            i = iblocking + 1;
            joinable = 1;
            //printf("A dépiler: %lu\n", n+1);
            //printf("Nouvelle ");
            //printPartition(p, n, sfpartition);
            //printf("Inverse ");
            //printPartition(p, pmax*64, sfpartitioninvert);
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
    free(work2);
    return nbests[pmax-1];
}



unsigned long schurNumberIterWithUnstack(unsigned long pmax, unsigned long *nbests, mpz_t iternum) {
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
    unsigned long i;
    unsigned long iblocking;
    unsigned long p;
    char isSumFree;
    char isAppendable;  //Vrai si il est possible de former une partition sans-somme contenant n
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    mp_size_t limballoc;    // Nombre de limbes alloué à chaque ensemble de sfpartition et à work0 et work1
    mp_size_t limbsize;     // Nombre de limbes utilisés par les ensembles de sfpartition
    mp_size_t wlimbsize;    // Nombre de limbes utilisés par work0
    unsigned int nmodbpl;
    unsigned int shift;
    mp_limb_t mask;
    mp_limb_t mask2;
    
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
    nmed = 0;
    nblockingmax = 0;
    isAppendable = 0;
    
    /*Itération jusqu'à énumérer toutes les partions sans-somme à au plus pmax huches*/
    while (n>0) {
        // Placer n+1 dans une des huches en conservant la sans-sommité
        while (i < p) {
            // Tester si l'ensemble obtenu en ajoutant n à la huche i est sans-somme
            mpz_add_ui(iternum, iternum, 1);
            //set = sfpartitioninvert[i];
            mpn_rshift(work1, sfpartitioninvert[i] + limballoc - limbsize, limbsize, shift);
            //set = sfpartition[i];
            mpn_and_n(work0, sfpartition[i], work1, wlimbsize);
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
//    schurNumberImprovedIterative(4, schurNumbers);
//    for (i=0; i<4; i++) {
//        printf("Schur Number of %i : %lu\n", i+1, schurNumbers[i]);
//    }
//    return 0;
//}
