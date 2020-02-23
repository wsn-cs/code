//
//  schurNumberWeakImprovedIterative.c
//  SchurNumber
//
//  Created by wsn-cs on 08/05/2019.
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

unsigned long schurNumberWeakIterWithUnstack(unsigned long pmax, unsigned long *nbests, mpz_t iternum) {
    /*
     Cette fonction calcule successivement les nombres de Schur faibles WS(p) pour p<= pmax.
     Elle remplit le tableau nbests, et compte le nombre d'itérations dans iternum.
     
     La partition est représentée comme un tableau de grands entiers sfpartition.
     Chaque grand entier représente un ensemble grâce à ses bits: si le bit k vaut 1,
     alors k appartient à l'ensemble; sinon il n'y appartient pas.
     
     Pour tester si il est possible d'ajouter n+1 à l'ensemble P sans-somme,
     on effectue un & entre les [n/2] premiers bits et les [n/2] derniers.
     Si le résultat est 0, il est possible d'ajouter n+1.
     
     Lorsqu'un entier n n'a pu être inséré dans aucune huches, on revient au sup sur les huches
     du plus petit n'≥[n/2]+1 tel que n' et n-n' appartiennent simultanément à cette huche.
     */
    unsigned long n;
    unsigned long nblocking;
    unsigned long nblockingmax;
    unsigned long nmed;
    unsigned long i;
    unsigned long iblocking;
    unsigned long p;
    char isSumFree;
    char isAppendable;  //Vrai si il est possible de former une partition sans-somme contenant n+1
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
    limbsize = 1;       // Nombre de limbes utilisés par chaque ensemble de la partition
    wlimbsize = 1;      // Taille en limbe de la seconde moitié
    nmodbpl = 1;        // Reste modulo mp_bits_per_limb de n
    shift = mp_bits_per_limb - nmodbpl;
    nmed = 0;           // Quotient de n par 2
    nblockingmax = 0;
    isAppendable = 0;
    
    /*Itération jusqu'à énumérer toutes les partions sans-somme à au plus pmax huches*/
    while (n>0) {
        // Placer n+1 dans une des huches en conservant le caractère faiblement sans-somme
        while (i < p) {
            // Tester si l'ensemble obtenu en ajoutant n+1 à la huche i est faiblement sans-somme
            mpz_add_ui(iternum, iternum, 1);
            mpn_rshift(work1, sfpartitioninvert[i] + limballoc - limbsize, limbsize, shift);
            mpn_and_n(work0, sfpartition[i], work1, wlimbsize);
            mask = (mp_limb_t)1 << (nmed % mp_bits_per_limb);
            work0[wlimbsize - 1] &= ~mask;    // Ce masque permet d'éliminer l'enventuelle somme double n/2 + n/2
            isSumFree = mpn_zero_p(work0, wlimbsize);
            if (isSumFree) {
                // Ajouter n+1 à la huche i
                mask = (mp_limb_t)1<<nmodbpl;
                sfpartition[i][limbsize -1] |= mask;
                mask = (mp_limb_t)1<<(shift - 1);
                sfpartitioninvert[i][limballoc - limbsize] |= mask;
                // Incrémenter n
                n++;
                nmodbpl = n%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                nmed = n>>1;
                wlimbsize = (nmed>>GMP_2EXP_NUMB_BITS) + 1;
                //nmed--;
                nblockingmax = 0;
                
                if (n > nbests[p-1]) {
                    // Mettre à jour le meilleur n obtenu.
                    nbests[p-1] = n;
                    
                    if (n >= mp_bits_per_limb*limballoc) {
                        // Augmentation de la taille
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
                    /* mp_bits_per_limb divisant n, il faut commencer un nouveau limbe.*/
                    for (i=0; i<pmax; i++) {
                        sfpartition[i][limbsize] = (mp_limb_t)0;
                        sfpartitioninvert[i][limballoc - limbsize -1] = (mp_limb_t)0;
                    }
                    limbsize++;
                }
                
                i = 0;
                isAppendable = 0;
                
            } else {
                /* Déterminer le plus petit entier nblocking ≥ n/2 bloquant pour la huche i,
                 i.e. tel que (n - nblocking + 1) appartienne à la huche i.*/
                nblocking = n - mpn_rscan1(work0, nmed+1);
                if (nblocking > nblockingmax) {
                    /* Mettre à jour nblockingmax = max_{1≤i≤p} nblocking_i.*/
                    nblockingmax = nblocking;
                    iblocking = i;
                }
                i++;
            }
        }
        if (i == p && p < pmax) {
            // Cas particulier où il faut remplir une nouvelle huche
            mpz_add_ui(iternum, iternum, 1);
            // Ajouter n+1 à la nouvelle huche p.
            mask = (mp_limb_t)1<<nmodbpl;
            sfpartition[i][limbsize -1] |= mask; // Ajoute n+1 à la huche i
            mask = (mp_limb_t)1<<(shift - 1);
            sfpartitioninvert[i][limballoc - limbsize] |= mask;
            // Incrémenter n.
            n++;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            nmed = n>>1;
            wlimbsize = (nmed>>GMP_2EXP_NUMB_BITS) + 1;
            //nmed--;
            nblockingmax = 0;
            
            if (n > nbests[p]) {
                // Mettre à jour le meilleur n obtenu.
                nbests[p] = n;
                
                if (n >= 64*limballoc -1) {
                    // Augmenter la taille
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
                /* mp_bits_per_limb divisant n, il faut commencer un nouveau limbe.*/
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
            // Décrémenter n
            n--;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            nmed = n>>1;
            wlimbsize = (nmed>>GMP_2EXP_NUMB_BITS) + 1;
            //nmed--;
            mask = (mp_limb_t)1<<nmodbpl;
            i = 0;
            // Trouver la huche i contenant n
            while (i < p && !(sfpartition[i][limbsize-1] & mask)) {
                i++;
            }
            // Retirer n de cette huche
            sfpartition[i][limbsize-1] ^= mask;
            mask = (mp_limb_t)1<<(shift - 1);
            sfpartitioninvert[i][limballoc - limbsize] ^= mask;
            if (i == p-1) {
                // Cas où n appartient à la dernière huche non vide
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
            nmed = n>>1;
            wlimbsize = (nmed>>GMP_2EXP_NUMB_BITS) + 1;
            //nmed--;
            mask = (mp_limb_t)GMP_NUMB_MAX >> shift;
            mask2 = (mp_limb_t)GMP_NUMB_MAX << (shift);
            for (i=0; i<pmax; i++) {
                sfpartition[i][limbsize-1] &= mask;
                sfpartitioninvert[i][limballoc - limbsize] &= mask2;
            }
            
            while (mpn_zero_p(sfpartition[p-1], limbsize)) {
                // Tant que la dernière huche potentiellement non vide est vide
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
