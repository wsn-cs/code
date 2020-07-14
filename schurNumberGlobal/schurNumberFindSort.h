//
//  schurNumberFindSort.h
//  SchurNumber
//
//  Created by rubis on 24/06/2020.
//  Copyright © 2020 rubis. All rights reserved.
//
//  Cet en-tête définit des fonctions permettant de trouver rapidement grâce à des méthodes dichotomiques
//  un élément dans un tableau trié.
//  Il faut au préalable définir des macros INDEXING_PARAMS, INDEXING_FUNC et CMP_FUNC pour décrire
//  respectivement l'indexage de la liste et la comparaison à effectuer.
//  La fonction CMP_FUNC doit renvoyer 0 en cas d'égalité.

#if !defined(FIND_INF_INDEX_FUNC) || !defined(FIND_SUP_INDEX_FUNC) || !defined(FIND_EXTREMAL_INDICES_FUNC)
#error "Manque les noms des fonctions de recherche."
#endif


#ifndef FUNC_PARAMS
#error "Aucun paramètre spécial défini pour les fonctions de recherche"
#endif
#ifndef FUNC_ARGS
#error "Aucun argument défini pour les fonctions de recherche"
#endif

#ifndef INDEXING_FUNC
#define INDEXING_FUNC(i) (i)
#endif

#ifndef CMP_FUNC
#error "Aucune fonction de comparaison définie."
#endif

#ifndef schurNumberFindSort_h
#define schurNumberFindSort_h

//typedef size_t (^indexing_block_t)(size_t i);

//typedef long (^cmp_block_t)(mp_limb_t *s1, mp_limb_t *s2);

__attribute__((flatten)) static inline unsigned long FIND_INF_INDEX_FUNC(mp_limb_t *set, mp_limb_t *set_array, FUNC_PARAMS, size_t imin, size_t isup) {
    /* Renvoie le plus petit indice imin ≤ i1 ≤ isup tel que cmp_b(set_array[ind_b(i1)], set) == 0, en utilisant une méthode dichotomique.
     Cela suppose que cmp_b(set_array[ind_b(isup)], set) == 0 et que la fonction croît avec l'indice. */
    long cmp_res;
    CMP_FUNC(cmp_res, &set_array[INDEXING_FUNC(imin)], set);
    if (!cmp_res) {
        return imin;
    }
    
    size_t i1 = imin;
    size_t i2 = isup;
    
    while (i2 - i1 > 1) {
        unsigned long imedian = (i2 + i1) >> 1;
        
        CMP_FUNC(cmp_res, &set_array[INDEXING_FUNC(imedian)], set);
        if (cmp_res) {
            // set_array[imedian] < set
            i1 = imedian;
        } else {
            // set_array[imedian] == set
            i2 = imedian;
            // Que cmp_b(set_array[i2], set) == 0 est un invariant de boucle.
        }
    }
    return i2;
}

//__attribute__((flatten)) static inline unsigned long find_sup_index(mp_limb_t *set, mp_limb_t *set_array, mp_size_t limbsize, indexing_block_t ind_b, size_t imin, size_t isup) {
//    /* Renvoie le plus petit indice imin ≤ i2 ≤ isup tel que set_array[limbsize * sorted_indexes[i2]] > set, sachant que set_array est un tableau de grands entiers de taille limbsize.
//     Cette fonction utilise une méthode dichotomique, et suppose que set_array[limbsize * sorted_indexes[imin]] ≤ set. */
//    size_t i1 = imin;
//    size_t i2 = isup;
//
//    while (i2 - i1 > 1) {
//        unsigned long imedian = (i2 + i1) >> 1;
//
//        if (mpn_cmp(&set_array[ind_b(imedian)], set, limbsize) > 0) {
//            // set_array[imedian] > set
//            i2 = imedian;
//        } else {
//            // set_array[imedian] <= set
//            i1 = imedian;
//            // Que set_array[i1] <= set est un invariant de boucle.
//        }
//    }
//    return i2;
//}

__attribute__((flatten)) static inline unsigned long FIND_SUP_INDEX_FUNC(mp_limb_t *set, mp_limb_t *set_array, FUNC_PARAMS, size_t imin, size_t isup) {
    /* Renvoie le plus petit indice imin ≤ i1 ≤ isup tel que cmp_b(set_array[ind_b(i1)], set) == 0, en utilisant une méthode dichotomique.
     Cela suppose que cmp_b(set_array[ind_b(isup)], set) == 0 et que la fonction croît avec l'indice. */
    long cmp_res;
    CMP_FUNC(cmp_res, &set_array[INDEXING_FUNC(isup)], set);
    if (!cmp_res) {
        return isup;
    }
    
    size_t i1 = imin;
    size_t i2 = isup;
    
    while (i2 - i1 > 1) {
        unsigned long imedian = (i2 + i1) >> 1;
        
        CMP_FUNC(cmp_res, &set_array[INDEXING_FUNC(imedian)], set);
        if (cmp_res) {
            // set_array[imedian] > set
            i2 = imedian;
        } else {
            // set_array[imedian] == set
            i1 = imedian;
            // Que cmp_b(set_array[i1], set) == 0 est un invariant de boucle.
        }
    }
    return i2;
}

//__attribute__((flatten)) static inline size_t find_extremal_indices(mp_limb_t *set, mp_limb_t *set_array, mp_size_t limbsize, indexing_block_t ind_b, size_t *i1_p, size_t *i2_p) {
//    /* Renvoie les deux indices extremaux i1 ≤ i2 tels que set_array[ind_b(i1)] < set < set_array[ind_b(i2)] , sachant que set_array est un tableau de grands entiers de taille limbsize.
//     Cette fonction utilise une méthode dichotomique. */
//
//    if (mpn_cmp(&set_array[ind_b(*i1_p)], set, limbsize) > 0) {
//        // Tout élément de set_array est > set
//        *i2_p = *i1_p;
//        return 0;
//    }
//
//    if (mpn_cmp(&set_array[ind_b(*i2_p)], set, limbsize) < 0) {
//        // Tout élément de set_array est < set
//        *i1_p = *i2_p;
//        return 0;
//    }
//
//    while (*i2_p - *i1_p > 1) {
//        unsigned long imedian = (*i2_p + *i1_p) >> 1;
//
//        int cmp = mpn_cmp(&set_array[ind_b(imedian)], set, limbsize);
//
//        if (!cmp) {
//            // set_array[imedian] == set
//            *i1_p = find_inf_index(set, set_array, limbsize, ind_b, *i1_p, imedian);
//            *i2_p = find_inf_index(set, set_array, limbsize, ind_b, imedian, *i2_p);
//            break;
//        } else if (cmp > 0) {
//            // set_array[imedian] > set
//            *i2_p = imedian;
//        } else {
//            // set_array[imedian] < set
//            *i1_p = imedian;
//        }
//    }
//
//    return *i2_p - *i1_p;
//}

__attribute__((flatten)) static inline size_t FIND_EXTREMAL_INDICES_FUNC(mp_limb_t *set, mp_limb_t *set_array, FUNC_PARAMS, size_t *i1_p, size_t *i2_p) {
    /* Place les deux indices extremaux i1 ≤ i2 tels que cmp_b(set_array[ind_b(i1)], set) ==  cmp_b(set_array[ind_b(i2)], set) == 0 , en utilisant une méthode dichotomique, et renvoie le nombre d'éléments tels que cmp_b(set_array[ind_b(i)], set) == 0.
     Cela suppose que la fonction cmp_b croît avec les indices, et peut donc prendre des valeurs positives et négatives. */
    long cmp_res;
    
    CMP_FUNC(cmp_res, &set_array[INDEXING_FUNC(*i1_p)], set);
    if (cmp_res > 0) {
        // Tout élément de set_array est > set
        *i2_p = *i1_p;
        return 0;
    } else if (!cmp_res) {
        // set_array[ind_b(*i1_p)] == set
        *i2_p = FIND_SUP_INDEX_FUNC(set, set_array, FUNC_ARGS, *i1_p, *i2_p);
        return *i2_p - *i1_p + 1;
    }
    
    CMP_FUNC(cmp_res, &set_array[INDEXING_FUNC(*i2_p)], set);
    if (cmp_res < 0) {
        // Tout élément de set_array est < set
        //printf("%li %zu\n", cmp_b(&set_array[ind_b(*i2_p)], set), *i2_p);
        *i1_p = ++(*i2_p);
        return 0;
    } else if (!cmp_res) {
        // set_array[ind_b(*i2_p)] == set
        *i1_p = FIND_INF_INDEX_FUNC(set, set_array, FUNC_ARGS, *i1_p, *i2_p);
        return *i2_p - *i1_p + 1;
    }
    
    size_t card = 0;
    
    while (*i2_p - *i1_p > 1) {
        unsigned long imedian = (*i2_p + *i1_p) >> 1;
        
        CMP_FUNC(cmp_res, &set_array[INDEXING_FUNC(imedian)], set);
        
        if (!cmp_res) {
            // set_array[imedian] == set
            *i1_p = FIND_INF_INDEX_FUNC(set, set_array, FUNC_ARGS, *i1_p, imedian);
            *i2_p = FIND_SUP_INDEX_FUNC(set, set_array, FUNC_ARGS, imedian, *i2_p);
            card = *i2_p - *i1_p + 1;
            break;
        } else if (cmp_res > 0) {
            // set_array[imedian] > set
            *i2_p = imedian;
        } else {
            // set_array[imedian] < set
            *i1_p = imedian;
        }
    }
    
    // Si set_array contient au moins un élément == set, alors [i1, i2] correspond exactement à l'intervalle cherché.
    // Sinon, i2 - i1 < 2 et set_array[i1] < set < set_array[i2]
    
    return card;
}

#endif /* schurNumberFindSort_h */
