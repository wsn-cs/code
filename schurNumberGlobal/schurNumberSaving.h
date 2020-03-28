//
//  schurNumberSaving.h
//  SchurNumber
//
//  Created by rubis on 14/03/2020.
//  Copyright © 2020 rubis. All rights reserved.
//
// Ce fichier déclare les fonctions permettant la sauvegarde des résultats intermédiaires dans un fichier temporaire

#ifndef schurNumberSaving_h
#define schurNumberSaving_h

#include <stdio.h>
#include <unistd.h>
#include <pthread.h>
#include "schurNumberIO.h"
#include "schurNumberPartitionStruc.h"

struct schur_number_intermediate_save_struc {
    unsigned long p;            // Nombre d'ensembles par partition
    
    unsigned long nbest;        // Plus grande taille de partition trouvée
    mp_limb_t **best_partition; // Partition réalisant cet optimum
    char toprint;               // Booléen spécifiant si cette partition doit être affichées
    
    unsigned long nbest_estimated;  // Plus grande taille de partition estimée
    unsigned long estimated_iternum;
    
    size_t iremainding;             // Nombre de partitions initiales restant à explorer
    unsigned long *estimated_iternum_initial_partition; // Tableau contenant le nombre d'itérations estimées pour chaque partition initiale
    
    int fd;                     // Descripteur du fichier où effectuer les sauvegardes intermédiaires
    char *filename;             // Nom du fichier temporaire
    pthread_mutex_t mutex_s;    // Mutex liée à l'écriture dans fd
};

typedef struct schur_number_intermediate_save_struc schur_number_intermediate_save_t;

void schurNumberSaveInit(schur_number_intermediate_save_t *save, unsigned long p, size_t part_pool_count, unsigned long n0);
void schurNumberSaveDealloc(schur_number_intermediate_save_t *save);

unsigned long schurNumberSaveUpgrade(schur_number_intermediate_save_t *save, unsigned long nbest, mp_limb_t **partition);

void schurNumberSaveToFile(schur_number_intermediate_save_t *save);

#endif /* schurNumberSaving_h */
