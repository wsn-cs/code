//
//  schurNumberSaving.c
//  SchurNumber
//
//  Created by rubis on 14/03/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberSaving.h"
#include <fcntl.h>

int schurNumberSaveFileOpen(schur_number_save_file_t *save_file_p) {
    char template[31] = "./schur_number_tmp_save_XXXXXX";
    int fd = mkstemp(template);
    
    if (fd != -1) {
        save_file_p->filename = calloc(1, sizeof(template));
        snprintf(save_file_p->filename, sizeof(template), "%s", template);
        
        dprintf(fd, "Nombre de partitions initiales inexplorées : ");
        save_file_p->irem_offset = lseek(fd, 0, SEEK_CUR);
        
        dprintf(fd, "%16lu\nEstimation du nombre d'itérations restants : ", (unsigned long)0);
        save_file_p->iternum_offset = lseek(fd, 0, SEEK_CUR);
        
        dprintf(fd, "%16lu\nEtat de progression : ", (unsigned long)0);
        save_file_p->percentage_offset = lseek(fd, 0, SEEK_CUR);
        
        dprintf(fd, "%16lu%%\nTaille maximale trouvée : ", (unsigned long)0);
        save_file_p->nbest_offset = lseek(fd, 0, SEEK_CUR);
    }
    
    save_file_p->fd = fd;
    
    return fd;
}

int schur_number_save_alloc(schur_number_intermediate_save_t *save, unsigned long p, unsigned long n0) {
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

    mpz_init2(save->estimated_iternum, nbest_estimated - n0 + 1);
    mpz_ui_pow_ui(save->estimated_iternum, p, nbest_estimated - n0 + 1);
    mpz_init_set(save->total_estimated_iternum, save->estimated_iternum);
    
    save->iremainding = 0;
    mpz_init_set(save->branch_estimated_iternum, save->estimated_iternum);
    pthread_key_create(&save->key, mpz_clear);
    
    save->tick = 7;
    
    pthread_mutex_init(&save->mutex_s, NULL);
    
    return fd;
}

void schurNumberSaveFileClose(schur_number_save_file_t *save_file_p) {
    
    close(save_file_p->fd);
    remove(save_file_p->filename);
    free(save_file_p->filename);
}

void schur_number_save_dealloc(schur_number_intermediate_save_t *save) {
    mpz_clear(save->estimated_iternum);
    mpz_clear(save->total_estimated_iternum);
    mpz_clear(save->branch_estimated_iternum);
    
    for (unsigned long j = 0; j < save->p; j++) {
        free(save->best_partition[j]);
    }
    free(save->best_partition);
    
    mpz_t *estimated_branch_iternum_p = pthread_getspecific(save->key);
    mpz_clear(*estimated_branch_iternum_p);
    pthread_key_delete(save->key);

    schurNumberSaveFileClose(&save->file);
    
    pthread_mutex_destroy(&save->mutex_s);
}

void schur_number_save_thread_register(schur_number_intermediate_save_t *save) {
    /* Tout thread doit appeler cette fonction pour s'enregistrer avant d'effectuer la moindre sauvegarde.
     La variable thread_estimated_iternum doit pointer vers un entier mpz déjà initialisée. Celui-ci est enregistré dans la clef save->key.*/
    
    mpz_t *thread_estimated_iternum_p = calloc(1, sizeof(mpz_t));
    mpz_init(*thread_estimated_iternum_p);
    
    pthread_setspecific(save->key, thread_estimated_iternum_p);
}

void schur_number_save_partition_pool_register(schur_number_intermediate_save_t *save, size_t part_pool_count, unsigned long n0) {
    /*Cette fonction enregistre dans save la présence d'un partitionpool.*/
    
    mpz_ui_pow_ui(save->branch_estimated_iternum, save->p, save->nbest_estimated - n0 + 1);
    mpz_mul_ui(save->estimated_iternum, save->branch_estimated_iternum, part_pool_count);
    mpz_set(save->total_estimated_iternum, save->estimated_iternum);

    save->n0 = n0;
    save->iremainding = part_pool_count;
}

void schur_number_save_newexploration_register(schur_number_intermediate_save_t *save) {
    /* Cette fonction est à appeler à chaque fois q'un thread entame l'exploration à partir d'une nouvelle partition.
     Elle enregistre l'indice associé au thread dans la clef save->key. */
    pthread_key_t key = save->key;
    
    mpz_t *estimated_branch_iternum_p = pthread_getspecific(key);
    
    mpz_submul_ui(save->estimated_iternum, *estimated_branch_iternum_p, 1);
    mpz_set(*estimated_branch_iternum_p, save->branch_estimated_iternum);
    
    save->iremainding --;
}

void schurNumberEstimatedRemaindingIteration(mpz_t iternum_estimated, unsigned long p, mp_limb_t **partition, unsigned long n, unsigned long n0, unsigned long nbound) {
    /* Cette fonction estime le nombre d'itérations qu'il restera à parcourir sur l'arbre.
     Cette estimation s'obtient en dénombrant le nombre de combinaison qu'il reste à effectuer pour énumérer les partitions de nbound sachant que nous sommes parvenus à partition. */
    
    if (nbound <= n) {
        nbound = n+1;
    }
    
    mpz_t num_configuration;        // Nombre de configurations possibles en partant de l'étape n pour atteindre nbound
    mpz_init2(num_configuration, nbound - n - 1);
    mpz_ui_pow_ui(num_configuration, p, nbound - n - 1);
    
    mpz_set_ui(iternum_estimated, 0);
    
    for (unsigned long i = n; i > n0; i--) {
        unsigned long j = 0;
        while (j < p && !GET_POINT(partition[j], i)) {
            j++;
        }
        mpz_addmul_ui(iternum_estimated, num_configuration, j);
        mpz_mul_ui(num_configuration, num_configuration, p);
    }
}

unsigned long schur_number_save_best_upgrade(schur_number_intermediate_save_t *save, unsigned long n, mp_limb_t **partition) {
    /* Met à jour la meilleure partition trouvée si il y a lieu.
     La fonction renvoie le nbest mis à jour, ce qui permet une éventuelle synchronisation entre threads. */
    
    char should_update = 0;
    
    pthread_mutex_lock(&save->mutex_s);
    
    unsigned long nbest = save->nbest;
    
    if (n > nbest) {
        
        if (partition) {
            // Sauvegarde de partition dans best_partition
            mp_size_t limbsize = INTEGER_2_LIMBSIZE(n);
            unsigned long p = save->p;
            
            for (unsigned long i = 0; i < p; i++) {
                mpn_copyd(save->best_partition[i], partition[i], limbsize);
            }
            
            should_update = 1;
        }
        
        nbest = n;
        save->nbest = nbest;
        save->toprint = 1;
    }
    
    pthread_mutex_unlock(&save->mutex_s);
    
    if (should_update) {
        schur_number_save_progression_update(save, n, partition);
    }
    
    return nbest;
}

void schurNumberSaveToFile(schur_number_intermediate_save_t *save) {
    /* Cette fonction sauvegarde l'état de la recherche dans le fichier fd. */
    
    schur_number_save_file_t save_file = save->file;
    int fd = save_file.fd;
    
    // Ecriture de iremainding
    lseek(fd, save_file.irem_offset, SEEK_SET);
    dprintf(fd, "%16lu", save->iremainding);
    
    // Ecriture de estimated_iternum
    lseek(fd, save_file.iternum_offset, SEEK_SET);
    mpz_ptr estimated_iternum = save->estimated_iternum;
    double rem_iternum = mpz_get_d(estimated_iternum);
    
    if (mpz_cmp_ui(estimated_iternum, 65536) > 0) {
        dprintf(fd, "%16.2e", rem_iternum);
    } else {
        unsigned long digits = mpz_get_ui(estimated_iternum);
        dprintf(fd, "%16lu", digits);
    }
    
    // Ecriture de pourcentage de progression
    lseek(fd, save_file.percentage_offset, SEEK_SET);
    double total_iternum = mpz_get_d(save->total_estimated_iternum);
    dprintf(fd, "%16.2f", 100 * rem_iternum / total_iternum);
    
    if (save->toprint) {
        unsigned long nbest = save->nbest;
        
        // Ecriture de nbest
        lseek(fd, save_file.nbest_offset, SEEK_SET);
        dprintf(fd, "%24lu\n", nbest);
        
        // Ecriture de best_partition
        schur_number_dprint_partition(fd, save->p, nbest, save->best_partition);
        save->toprint = 0;
    }
    
    fsync(fd);
}

unsigned long schur_number_save_progression_update(schur_number_intermediate_save_t *save, unsigned long n, mp_limb_t **partition) {
    /* Cette fonction met à jour les informations liées à la progression du programme.
     Elle renvoie le nbest de save, ce qui permet une éventuelle synchronisation entre threads. */
    
    // Nombre d'itérations précédemment estimé pour la branche courante
    mpz_t *estimated_branch_iternum_p = pthread_getspecific(save->key);
    mpz_t executed_iternum;
    mpz_init_set(executed_iternum, *estimated_branch_iternum_p);
    
    // Nombre d'itérations restant
    schurNumberEstimatedRemaindingIteration(*estimated_branch_iternum_p, save->p, partition, n, save->n0, save->nbest_estimated);
    
    // Nombre d'itérations déjà effectuées
    mpz_sub(executed_iternum, executed_iternum, *estimated_branch_iternum_p);
    
    pthread_mutex_lock(&save->mutex_s);
    mpz_sub(save->estimated_iternum, save->estimated_iternum, executed_iternum);
    save->tick--;
    save->tick %= 8;
    if (!save->tick) {
        schurNumberSaveToFile(save);
    }
    pthread_mutex_unlock(&save->mutex_s);

    return save->nbest;
}
