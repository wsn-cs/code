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

int schurNumberSaveFileOpen(schur_number_save_file_t *save_file_p) {
    char template[31] = "./schur_number_tmp_save_XXXXXX";
    int fd = mkstemp(template);
    
    if (fd != -1) {
        asprintf(&save_file_p->filename, "%s", template);
        
        dprintf(fd, "Nombre de partitions initiales inexplorées :\t");
        save_file_p->irem_offset = lseek(fd, 0, SEEK_CUR);
        
        dprintf(fd, "%20lu\nEstimation du nombre d'itérations restants :\t", (unsigned long)0);
        save_file_p->iternum_offset = lseek(fd, 0, SEEK_CUR);
        
        dprintf(fd, "%20lu\nTaille maximale trouvée :\t", (unsigned long)0);
        save_file_p->nbest_offset = lseek(fd, 0, SEEK_CUR);
    }
    
    return fd;
}

int schurNumberSaveAlloc(schur_number_intermediate_save_t *save, unsigned long p, unsigned long n0) {
    /* Cette fonction crée un fichier temporaire où sauvegarder les partitions intermédiaires générées au cours de l'exécution du programme.
     Elle renvoie le descripteur de fichier associé. */
    int fd = schurNumberSaveFileOpen(&save->file);
    
    if (fd == -1) {
        fprintf(stderr, "Impossible de créer un fichier temporaire pour sauvegarder les résultats\n");
    } else {
        fprintf(stderr, "Fichier temporaire pour sauvegarder les résultats créé : %s\n", save->file.filename);
    }
    
    save->p = p;
    save->n0 = n0;
    
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
    save->estimated_iternum = estimated_iternum;
    
    save->iremainding = 0;
    save->branch_estimated_iternum = estimated_iternum;
    pthread_key_create(&save->key, NULL);
    
    save->tick = 7;
    
    pthread_mutex_init(&save->mutex_s, NULL);
    
    return fd;
}

void schurNumberSaveFileClose(schur_number_save_file_t *save_file_p) {
    
    close(save_file_p->fd);
    remove(save_file_p->filename);
    free(save_file_p->filename);
}

void schurNumberSaveDealloc(schur_number_intermediate_save_t *save) {
    
    for (unsigned long j = 0; j < save->p; j++) {
        free(save->best_partition[j]);
    }
    free(save->best_partition);
    
    pthread_key_delete(save->key);
    
    pthread_mutex_destroy(&save->mutex_s);
}

void schurNumberSavePartitionPoolRegister(schur_number_intermediate_save_t *save, size_t part_pool_count, unsigned long n0) {
    /*Cette fonction enregistre dans save la présence d'un partitionpool.*/
    
    unsigned long estimated_iternum = fast_powl(save->p, save->nbest_estimated - n0);
    
    save->branch_estimated_iternum = estimated_iternum;
    save->estimated_iternum = estimated_iternum * part_pool_count;
    save->iremainding = part_pool_count;
}

void schurNumberSaveNewExplorationRegister(schur_number_intermediate_save_t *save, unsigned long index) {
    /* Cette fonction est à appeler à chaque fois q'un thread entame l'exploration à partir d'une nouvelle partition.
     Elle enregistre l'indice associé au thread dans la clef save->key. */
    pthread_key_t key = save->key;
    
    unsigned long estimated_branch_iternum = pthread_getspecific(key);
    
    save->estimated_iternum -= estimated_branch_iternum;
    pthread_setspecific(key, save->branch_estimated_iternum);
    save->iremainding --;
}

unsigned long schurNumberEstimatedRemaindingIteration(unsigned long p, mp_limb_t **partition, unsigned long n, unsigned long n0, unsigned long nbound) {
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

unsigned long schurNumberSaveBestUpgrade(schur_number_intermediate_save_t *save, unsigned long n, mp_limb_t **partition) {
    /* Met à jour la meilleure partition trouvée si il y a lieu.
     La fonction renvoie le nbest mis à jour, ce qui permet une éventuelle synchronisation entre threads. */
    
    pthread_mutex_lock(&save->mutex_s);
    
    unsigned long nbest = save->nbest;
    
    if (n > nbest) {
        unsigned long p = save->p;
        
        if (partition) {
            // Sauvegarde de partition dans best_partition
            mp_size_t limbsize = INTEGER_2_LIMBSIZE(n);
            
            for (unsigned long i = 0; i < p; i++) {
                mpn_copyd(save->best_partition[i], partition[i], limbsize);
            }
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
    
    schur_number_save_file_t save_file = save->file;
    int fd = save_file.fd;
    
    // Ecriture de iremainding
    lseek(fd, save_file.irem_offset, SEEK_SET);
    dprintf(fd, "%20lu", save->iremainding);
    
    // Ecriture de estimated_iternum
    lseek(fd, save_file.iternum_offset, SEEK_SET);
    dprintf(fd, "%20lu", save->estimated_iternum);
    
    if (save->toprint) {
        unsigned long nbest = save->nbest;
        
        // Ecriture de nbest
        lseek(fd, save_file.nbest_offset, SEEK_SET);
        dprintf(fd, "%20lu\n", nbest);
        
        // Ecriture de best_partition
        schur_number_dprint_partition(fd, save->p, nbest, save->best_partition);
        save->toprint = 0;
    }
}

unsigned long schurNumberSaveProgressionUpdate(schur_number_intermediate_save_t *save, unsigned long n, mp_limb_t **partition) {
    /* Cette fonction met à jour les informations liées à la progression du programme.
     Elle renvoie le nbest de save, ce qui permet une éventuelle synchronisation entre threads. */
    
    // Nombre d'itérations précédemment estimé pour la branche courante
    unsigned long estimated_branch_iternum = pthread_getspecific(save->key);
    
    // Nombre d'itérations restant
    unsigned long estimated_rem_iternum = schurNumberEstimatedRemaindingIteration(save->p, partition, n, save->n0, save->nbest_estimated);
    pthread_setspecific(save->key, estimated_rem_iternum);
    
    // Nombre d'itérations déjà effectuées
    unsigned long executed_iternum = estimated_branch_iternum - estimated_rem_iternum;
    
    pthread_mutex_lock(&save->mutex_s);
    save->estimated_iternum -= executed_iternum;
    save->tick--;
    save->tick %= 8;
    if (!save->tick) {
        schurNumberSaveToFile(save);
    }
    pthread_mutex_unlock(&save->mutex_s);

    return save->nbest;
}
