//
//  schurNumberSaving.c
//  SchurNumber
//
//  Created by rubis on 14/03/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberSaving.h"

static inline unsigned long fast_powl(unsigned long base, unsigned long exp) {
    /* Renvoie base^exp en utilisant l'exponentiation rapide. */
    unsigned long cur_exp = exp;
    unsigned long res = 1;
    unsigned long mul = base;
    
    while (cur_exp) {
        if (cur_exp & (unsigned long)1) {
            res *= mul;
        }
        mul *= mul;
        cur_exp >>= 1;
    }
    
    return res;
}

void schurNumberSaveInit(schur_number_intermediate_save_t *save, unsigned long p, unsigned long part_pool_count, unsigned long n0) {
    /* Cette fonction crée un fichier temporaire où sauvegarder les partitions intermédiaires générées au cours de l'exécution du programme.
     Elle renvoie le descripteur de fichier associé. */
    char template[31] = "./schur_number_tmp_save_XXXXXX";
    int fd = mkstemp(template);
    
    if (fd == -1) {
        fprintf(stderr, "Impossible de créer un fichier temporaire pour sauvegarder les résultats\n");
    } else {
        fprintf(stderr, "Fichier temporaire pour sauvegarder les résultats créé : %s\n", template);
    }
    
    save->fd = fd;
    asprintf(&save->filename, "%s", template);
    
    save->p = p;
    
    unsigned long nbest_estimated = 1<<(2 * p - 1);     // borne sur S(p)
    
    // Allocation de best_partition
    save->nbest = 0;
    mp_size_t limballoc = INTEGER_2_LIMBSIZE(nbest_estimated);
    save->best_partition = calloc(sizeof(mp_limb_t *), p);
    for (unsigned long j = 0; j < p; j++) {
        save->best_partition[j] = calloc(sizeof(mp_limb_t), limballoc);
    }
    save->toprint = 0;
    
    // Estimation du nombre d'itérations
    save->nbest_estimated = nbest_estimated;
    unsigned long estimated_iternum = fast_powl(p, nbest_estimated - n0);
    save->estimated_iternum = estimated_iternum * part_pool_count;
    
    save->iremainding = 0;
    save->estimated_iternum_initial_partition = NULL;
    
    if (part_pool_count > 1) {
        save->iremainding = part_pool_count;
        save->estimated_iternum_initial_partition = calloc(sizeof(unsigned long), part_pool_count);
        
        for (unsigned long i = 0; i < part_pool_count; i++) {
            save->estimated_iternum_initial_partition[i] = estimated_iternum;
        }
    }
    
    pthread_mutex_init(&save->mutex_s, NULL);
}

void schurNumberSaveDealloc(schur_number_intermediate_save_t *save) {
    
    for (unsigned long j = 0; j < save->p; j++) {
        free(save->best_partition[j]);
    }
    free(save->best_partition);
    
    if (save->estimated_iternum_initial_partition) {
        free(save->estimated_iternum_initial_partition);
    }
    
    close(save->fd);
    remove(save->filename);
    free(save->filename);
    
    pthread_mutex_destroy(&save->mutex_s);
}

unsigned long schurNumberEstimatedRemainingIteration(unsigned long p, mp_limb_t **partition, unsigned long n, unsigned long n0, unsigned long nbound) {
    /* Cette fonction estime le nombre d'itérations qu'il restera à parcourir sur l'arbre.
     Cette estimation s'obtient en dénombrant le nombre de combinaison qu'il reste à effectuer pour énumérer les partitions de nbound sachant que nous sommes parvenus à partition. */
    
    unsigned long i = nbound;
    
    unsigned long num_configuration = 1; // Nombre de configurations possibles en partant de l'étape n pour atteindre nbound
    while (i > n) {
        num_configuration *= p;
        i--;
    }
    
    unsigned long iternum_estimated = 0;
    
    while (i > n0) {
        unsigned long j = 0;
        while (j < p && !GET_POINT(partition[j], i)) {
            j++;
        }
        iternum_estimated += (p - j) * num_configuration;
        num_configuration *= p;
        i--;
    }
    
    return iternum_estimated;
}

unsigned long schurNumberSaveUpgrade(schur_number_intermediate_save_t *save, unsigned long n, mp_limb_t **partition) {
    /* Met à jour la meilleure partition trouvée si il y a lieu.
     La fonction renvoie le nbest mis à jour, ce qui permet une éventuelle synchronisation entre threads. */
    
    pthread_mutex_lock(&save->mutex_s);
    
    unsigned long nbest = save->nbest;
    
    if (n > nbest) {
        unsigned long p = save->p;
        
        mp_size_t limbsize = INTEGER_2_LIMBSIZE(n);
        
        for (unsigned long i = 0; i < p; i++) {
            mpn_copyd(save->best_partition[i], partition[i], limbsize);
        }
        
        nbest = n;
        save->nbest = nbest;
        save->toprint = 1;
    }
    
    pthread_mutex_unlock(&save->mutex_s);
    
    return nbest;
}

void schurNumberSaveToFile(schur_number_intermediate_save_t *save) {
    /* Cette fonction sauvegarde l'état de la recherche dans le fichier fd. */
    
    pthread_mutex_lock(&save->mutex_s);
    
    int fd = save->fd;
    
    dprintf(fd, "\nNombre de partitions initiales restant : %lu\nNombre d'itérations restant estimées : %lu\n", save->iremainding, save->estimated_iternum);
    
    if (save->toprint) {
        unsigned long nbest = save->nbest;
        dprintf(fd, "Taille maximale trouvée : %lu\n", nbest);
        schur_number_dprint_partition(fd, save->p, nbest, save->best_partition);
        save->toprint = 0;
    }
    
    pthread_mutex_unlock(&save->mutex_s);
}
