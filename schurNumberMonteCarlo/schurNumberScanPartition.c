//
//  schurNumberScanPartition.c
//  schurNumberMonteCarlo
//
//  Created by Gabriel Merlin on 01/05/2019.
//

#include <stdio.h>
#include "schurNumberNestedMonteCarloHeader.h"

unsigned int schurNumberScanPartitionFromFile(char *filename, partition_t *partitionstruc) {
    /*Crée une partition à partir d'un fichier texte.
     Des virgules séparent les entiers au sein d'un même ensemble
     et des points séparent les ensembles entre eux.
     La fonction renvoie le nombre d'ensembles p de la partition, ou 0 si un problème survient.*/
    FILE *fp;
    int c;
    char *intstr;
    unsigned long m, mmodbpl, mdivbpl;
    size_t len, maxlen;// Taille de la chaîne de caractères codant un entier
    unsigned int p;  // Nombre d'ensembles de la partition
    unsigned long n; // Nombre d'entiers
    mp_size_t limballoc;
    mp_limb_t **partition, **partitioninvert;
    
    
    /*Ouverture du flux*/
    fp = fopen(filename, "r");
    
    if (!fp) {
        /*Le fichier ne peut être ouvert.*/
        return 0;
    }
    
    /*Lecture préliminaire comptant les virgules et les points*/
    p = 0;
    n = 0;
    c = getc(fp);
    len = 0;
    maxlen = 0;
    while (c != EOF) {
        switch (c) {
            case ',':
                n++;
                if (len > maxlen) {
                    maxlen = len;
                }
                len = 0;
                break;
            
            case '.':
                n++;
                p++;
                if (len > maxlen) {
                    maxlen = len;
                }
                len = 0;
                break;
                
            default:
                len++;
                break;
        }
        c = fgetc(fp);
    }
    rewind(fp);
    
    /*Initialisation de la partition*/
//    limballoc = 1 + (n / mp_bits_per_limb);
//    partition = calloc(p, sizeof(mp_limb_t *));
//    partitioninvert = calloc(p, sizeof(mp_limb_t *));
//    for (i=0; i<p; i++) {
//        partition[i] = calloc(limballoc, sizeof(mp_limb_t));
//        partitioninvert[i] = calloc(limballoc, sizeof(mp_limb_t));
//    }
    partition_init(p, n, partitionstruc);
    limballoc = partitionstruc->limballoc;
    
    /*Mise en place des éléments*/
    intstr = calloc(maxlen+1, sizeof(char));
    c = fgetc(fp);
    len = 0;
    //p = 0;
    partition = partitionstruc->partition;
    partitioninvert = partitionstruc->partitioninvert;
    while (c != EOF) {
        switch (c) {
            case ',':
                intstr[len] = '\0';
                m = atol(intstr)-1;
                /*Placer m dans la bonne partition*/
                mdivbpl = m / mp_bits_per_limb;
                mmodbpl = m % mp_bits_per_limb;
                //partition[p][mmodbpl] |= (mp_limb_t)1<<mrem;
                //partitioninvert[p][limballoc - mmodbpl - 1] |= (mp_limb_t)1 << (mp_bits_per_limb - mrem - 1);
                (*partition)[mdivbpl] |= (mp_limb_t)1<<mmodbpl;
                (*partitioninvert)[limballoc - mdivbpl - 1] |= (mp_limb_t)1 << (mp_bits_per_limb - mmodbpl - 1);
                len = 0;
                break;
                
            case '.':
                intstr[len] = '\0';
                m = atol(intstr)-1;
                /*Placer m dans la bonne partition*/
                mdivbpl = m / mp_bits_per_limb;
                mmodbpl = m % mp_bits_per_limb;
                //partition[p][mmodbpl] |= (mp_limb_t)1<<mrem;
                //partitioninvert[p][limballoc - mmodbpl - 1] |= (mp_limb_t)1 << (mp_bits_per_limb - mrem - 1);
                //p++;
                (*partition)[mdivbpl] |= (mp_limb_t)1<<mmodbpl;
                (*partitioninvert)[limballoc - mdivbpl - 1] |= (mp_limb_t)1 << (mp_bits_per_limb - mmodbpl - 1);
                partition++;
                partitioninvert++;
                len = 0;
                break;
                
            default:
                intstr[len] = c;
                len++;
                break;
        }
        c = fgetc(fp);
    }
    free(intstr);
    
    /*Sauvegarde dans partition_t*/
    //partitionstruc->limballoc = limballoc;
    partitionstruc->limbsize = limballoc;
    //partitionstruc->partition = partition;
    //partitionstruc->partitioninvert = partitioninvert;
    //partitionstruc->pmax = p;
    partitionstruc->p = p;
    partitionstruc->n = n;
    
    /*Fermeture du flux*/
    fclose(fp);
    
    return p;
}

unsigned long schurNumberPrintPartitionToFile(char *filename, unsigned int p, unsigned long n, mp_limb_t **partition) {
    /*Sauvegarde une partition dans un fichier texte.
     Des virgules séparent les entiers au sein d'un même ensemble
     et des points séparent les ensembles entre eux.
     La fonction renvoie le nombre d'entiers n, ou 0 si un problème survient.*/
    FILE *fp;
    unsigned int j;
    unsigned long i, imodbpl, idivbpl;
    mp_size_t limballoc;
    mp_limb_t **set;
    
    
    /*Ouverture du flux*/
    fp = fopen(filename, "w");
    
    if (!fp) {
        /*Le fichier ne peut être ouvert.*/
        return 0;
    }
    
    /*Initialisation des variables*/
    limballoc = 1 + (n / mp_bits_per_limb);
    
    j = 0;
    set = partition;
    while (j < p) {
        
        if (**set & (mp_limb_t)1) {
            fprintf(fp, "1");
        }
        
        for (i=1; i<n; i++) {
            imodbpl = i % mp_bits_per_limb;
            idivbpl = i / mp_bits_per_limb;
            
            if ((*set)[idivbpl] & ((mp_limb_t)1 << imodbpl)) {
                fprintf(fp, ",%lu", i);
            }
        }
        
        fprintf(fp, ".");
        set++;
        j++;
    }
    
    /*Fermeture du flux*/
    fclose(fp);
    
    return n;
}
