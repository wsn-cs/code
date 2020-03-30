//
//  schurNumberIO.c
//  schurNumberRecursive
//
//  Created by rubis on 28/01/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberIO.h"

#define ADD_POINT(set, x) set[(x) / mp_bits_per_limb] |= ((unsigned long)1 << ((x) % mp_bits_per_limb))

#define GET_POINT(set, x) (set[(x) / mp_bits_per_limb] & ((unsigned long)1 << ((x) % mp_bits_per_limb)))

unsigned long schurNumberGetSetMaximum(char *str) {
    char *ptr0 = str;
    char *ptr1 = str;
    
    unsigned long nmax = 0;
    while (*ptr1 != '\0') {
        unsigned long n = strtoul(ptr0, &ptr1, 10);
        if (ptr0 != ptr1) {
            // Chiffres valides trouvés dans ptr0
            if (n > nmax) {
                nmax = n;
            }
        }
        ptr0 = ptr1 + 1;
    }
    return nmax;
}

void schurNumberGetSet(char *str, mp_limb_t *set, mp_limb_t *setinvert, mp_size_t limballoc) {
    /*Remplit l'ensemble d'entiers set grâce à la chaîne de caractères str et renvoie le plus grand élément trouvé. La variable set doit pointée vers un tableau de taille suffisante pour contenir l'ensemble, et de même pour setinvert sauf si il s'agit du pointeur nul.*/
    char *ptr0 = str;
    char *ptr1 = str;
    
    unsigned long nalloc = mp_bits_per_limb * limballoc;
    
    while (*ptr1 != '\0') {
        unsigned long k = strtoul(ptr0, &ptr1, 10);
        if (ptr0 != ptr1) {
            // Chiffres valides trouvés dans ptr0
            ADD_POINT(set, k);
            if (setinvert) {
                ADD_POINT(setinvert, nalloc - k);
            }
        }
        ptr0 = ptr1 + 1;
    }
}

unsigned long schurNumberGetPartition(unsigned long p, char **str, mp_limb_t **partition, mp_limb_t **partitioninvert, mp_size_t limballoc) {
    /*Remplit la variable partition à partir du tableau de p chaînes de caractères *str. La variable partition doit pointée vers un tableau à p entrées, et de même pour partitioninvert sauf si il s'agit du pointeur nul.
     Spécifier une valeur de limballoc assure que chaque ensemble aura une taille d'au moins limballoc limbes.*/
    unsigned long n = 0;
    
    for (unsigned long i = 0; i < p; i++) {
        unsigned long m = schurNumberGetSetMaximum(str[i]);
        if (m > n) {
            n = m;
        }
    }
    
    mp_size_t limbsize = (n>>6) + 1;
    if (limbsize < limballoc) {
        limbsize = limballoc;
    }
    for (unsigned long i = 0; i < p; i++) {
        partition[i] = calloc(sizeof(mp_limb_t), limbsize);
        if (partitioninvert) {
            partitioninvert[i] = calloc(sizeof(mp_limb_t), limbsize);
            schurNumberGetSet(str[i], partition[i], partitioninvert[i], limbsize);
        } else {
            schurNumberGetSet(str[i], partition[i], NULL, limbsize);
        }
    }
    
    return n;
}

void schurNumberPrintSet(int fd, unsigned long n, mp_limb_t *set) {
    /*Affiche le contenu de l'ensemble d'entiers set inclus dans l'intervalle [1, n].*/
    
    for (unsigned long i = 1; i <= n; i++) {
        if (GET_POINT(set, i)) {
            dprintf(fd, " %lu", i);
        }
    }
}

void schur_number_dprint_partition(int fd, unsigned long p, unsigned long n, mp_limb_t **partition) {
    /*Affiche une partition à p ensembles de [1, n] dans le fichier spécifié par fd.*/
    dprintf(fd, "Partition:\n");
    
    for (unsigned long i=0; i<p; i++) {
        dprintf(fd, "\t");
        schurNumberPrintSet(fd, n, partition[i]);
        dprintf(fd, "\n");
    }
}
