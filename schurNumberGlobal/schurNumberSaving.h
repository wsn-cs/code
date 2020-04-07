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

struct schur_number_save_file_struc {
    int fd;                     // Descripteur du fichier où effectuer les sauvegardes intermédiaires
    char *filename;             // Nom du fichier temporaire
    
    off_t irem_offset;          // Incrément repérant la position où écrire iremainding
    off_t iternum_offset;       // Incrément repérant la position où écrire estimated_iternum
    off_t nbest_offset;         // Incrément repérant la position où écrire nbest
};

typedef struct schur_number_save_file_struc schur_number_save_file_t;

struct schur_number_intermediate_save_struc {
    unsigned long p;            // Nombre d'ensembles par partition
    unsigned long n0;           // Taille de la partition initiale
    
    unsigned long nbest;        // Plus grande taille de partition trouvée
    mp_limb_t **best_partition; // Partition réalisant cet optimum
    char toprint;               // Booléen spécifiant si cette partition doit être affichée
    
    unsigned long nbest_estimated;  // Plus grande taille de partition estimée
    unsigned long estimated_iternum;// Nombre estimé d'itérations restant à effectuer
    
    size_t iremainding;         // Nombre de partitions initiales restant à explorer
    unsigned long branch_estimated_iternum; // Nombre estimé d'itérations pour chaque branche
    pthread_key_t key;          // Clé indiquant à chaque thread le nombre d'itération lui restant sur sa branche
    
    unsigned char tick;         // Compteur modulo 8 permettant de savoir si une sauvegarde doit être faite dans le fichier temporaire
    
    schur_number_save_file_t file;  // Fichier où effectuer les sauvegardes intermédiaires
    pthread_mutex_t mutex_s;        // Mutex permettant le multi-threading
};

typedef struct schur_number_intermediate_save_struc schur_number_intermediate_save_t;

int schurNumberSaveAlloc(schur_number_intermediate_save_t *save, unsigned long p, unsigned long n0);
void schurNumberSaveDealloc(schur_number_intermediate_save_t *save);

void schurNumberSavePartitionPoolRegister(schur_number_intermediate_save_t *save, size_t part_pool_count, unsigned long n0);
void schurNumberSaveNewExplorationRegister(schur_number_intermediate_save_t *save, unsigned long index);

unsigned long schurNumberSaveBestUpgrade(schur_number_intermediate_save_t *save, unsigned long n, mp_limb_t **partition);

unsigned long schurNumberSaveProgressionUpdate(schur_number_intermediate_save_t *save, unsigned long n, mp_limb_t **partition);

#endif /* schurNumberSaving_h */
