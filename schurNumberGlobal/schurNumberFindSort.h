//
//  schurNumberFindSort.h
//  SchurNumber
//
//  Created by rubis on 24/06/2020.
//  Copyright © 2020 rubis. All rights reserved.
//
//  Cet en-tête définit des fonctions permettant de trouver rapidement grâce à des méthodes dichotomiques un élément dans un tableau trié.

#ifndef schurNumberFindSort_h
#define schurNumberFindSort_h

typedef size_t (^indexing_block_t)(size_t i);

__attribute__((flatten)) static inline unsigned long find_inf_index(mp_limb_t *set, mp_limb_t *set_array, mp_size_t limbsize, indexing_block_t ind_b, size_t imin, size_t isup) {
    /* Renvoie le plus grand indice imin ≤ i1 ≤ isup tel que set_array[limbsize * sorted_indexes[i1]] < set, sachant que set_array est un tableau de grands entiers de taille limbsize.
     Cette fonction utilise une méthode dichotomique, et suppose que set_array[limbsize * sorted_indexes[isup]] ≥ set. */
    size_t i1 = imin;
    size_t i2 = isup;
    
    while (i2 - i1 > 1) {
        unsigned long imedian = (i2 + i1) >> 1;
        
        if (mpn_cmp(&set_array[ind_b(imedian)], set, limbsize) < 0) {
            // set_array[imedian] < set
            i1 = imedian;
        } else {
            // set_array[imedian] >= set
            i2 = imedian;
            // Que set_array[i2] >= set est un invariant de boucle.
        }
    }
    return i1;
}

__attribute__((flatten)) static inline unsigned long find_sup_index(mp_limb_t *set, mp_limb_t *set_array, mp_size_t limbsize, indexing_block_t ind_b, size_t imin, size_t isup) {
    /* Renvoie le plus petit indice imin ≤ i2 ≤ isup tel que set_array[limbsize * sorted_indexes[i2]] > set, sachant que set_array est un tableau de grands entiers de taille limbsize.
     Cette fonction utilise une méthode dichotomique, et suppose que set_array[limbsize * sorted_indexes[imin]] ≤ set. */
    size_t i1 = imin;
    size_t i2 = isup;
    
    while (i2 - i1 > 1) {
        unsigned long imedian = (i2 + i1) >> 1;
        
        if (mpn_cmp(&set_array[ind_b(imedian)], set, limbsize) > 0) {
            // set_array[imedian] > set
            i2 = imedian;
        } else {
            // set_array[imedian] <= set
            i1 = imedian;
            // Que set_array[i1] <= set est un invariant de boucle.
        }
    }
    return i2;
}

__attribute__((flatten)) static inline size_t find_extremal_indices(mp_limb_t *set, mp_limb_t *set_array, mp_size_t limbsize, indexing_block_t ind_b, size_t *i1_p, size_t *i2_p) {
    /* Renvoie les deux indices extremaux i1 ≤ i2 tels que set_array[ind_b(i1)] < set < set_array[ind_b(i2)] , sachant que set_array est un tableau de grands entiers de taille limbsize.
     Cette fonction utilise une méthode dichotomique. */
    
    if (mpn_cmp(&set_array[ind_b(*i1_p)], set, limbsize) > 0) {
        // Tout élément de set_array est > set
        *i2_p = *i1_p;
        return 0;
    }
    
    if (mpn_cmp(&set_array[ind_b(*i2_p)], set, limbsize) < 0) {
        // Tout élément de set_array est < set
        *i1_p = *i2_p;
        return 0;
    }
    
    while (*i2_p - *i1_p > 1) {
        unsigned long imedian = (*i2_p + *i1_p) >> 1;
        
        int cmp = mpn_cmp(&set_array[ind_b(imedian)], set, limbsize);
        
        if (!cmp) {
            // set_array[imedian] == set
            *i1_p = find_inf_index(set, set_array, limbsize, ind_b, *i1_p, imedian);
            *i2_p = find_inf_index(set, set_array, limbsize, ind_b, imedian, *i2_p);
            break;
        } else if (cmp > 0) {
            // set_array[imedian] > set
            *i2_p = imedian;
        } else {
            // set_array[imedian] < set
            *i1_p = imedian;
        }
    }
    
    return *i2_p - *i1_p;
}

#endif /* schurNumberFindSort_h */
