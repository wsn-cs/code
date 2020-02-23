//
//  schurNumberAbbotHansonMonteCarlo.c
//  schurNumberMonteCarlo
//
//  Created by rubis on 01/05/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#include "schurNumberNestedMonteCarloHeader.h"

void AbbotHansonPartitionGeneration(partition_t *partitionstrucr, partition_t *partitionstrucs1, partition_t *partitionstrucs2) {
    /*Crée une nouvelle partition partitionr en fusionnant partitions1
     et partitions2 selon la méthode de AbbotHanson.
     Les partitions ne sont pas désallouées dans cette fonction,
     par contre, les partitions de partitionr sont allouées*/
    unsigned int i;
    unsigned long b;
    unsigned int pmax1, pmax2;
    unsigned int pr, p1, p2;
    unsigned long nr, n1, n2, N;
    unsigned long Nb, Nbmodbpl, Nbdivbpl;
    mp_size_t limballocr, limballoc1, limballoc2;
    mp_size_t limbsizer, limbsize1, limbsize2;
    mp_limb_t **partitionr, **partition1, **partition2;
    mp_limb_t **partitioninvertr, **partitioninvert1, **partitioninvert2;
    mp_limb_t *setr, *sets;
    mp_limb_t *work0, *work1;
    
    /*Initialistion des variables correspondant à partitionstrucs1*/
    n1 = partitionstrucs1->n;
    pmax1 = partitionstrucs1->pmax;
    p1 = partitionstrucs1->p;
    limballoc1 = partitionstrucs1->limballoc;
    limbsize1 = partitionstrucs1->limbsize;
    partition1 = partitionstrucs1->partition;     // Tableau contenant les ensembles de la partition
    partitioninvert1 = partitionstrucs1->partitioninvert; // Tableau contenant les inverses des ensembles de la partition
    
    /*Initialistion des variables correspondant à partitionstrucs2*/
    n2 = partitionstrucs2->n;
    pmax2 = partitionstrucs2->pmax;
    p2 = partitionstrucs2->p;
    limballoc2 = partitionstrucs2->limballoc;
    limbsize2 = partitionstrucs2->limbsize;
    partition2 = partitionstrucs2->partition;     // Tableau contenant les ensembles de la partition
    partitioninvert2 = partitionstrucs2->partitioninvert; // Tableau contenant les inverses des ensembles de la partition
    
    /*Initialisation des variables de partitionstrucr*/
    N = 2*n1 + 1;
    nr = n2 * N + n1;
    pr = p1 + p2;
    limballocr = limballoc1 + limballoc2;
    partitionstrucr->pmax = pmax1 + pmax2;
    partitionstrucr->p = pr;
    partitionstrucr->n = nr;
    partitionstrucr->limballoc = limballocr;
    partitionstrucr->limbsize = nr / GMP_NUMB_BITS;
    /*Allocation des partitions*/
    work0 = calloc(limballocr, sizeof(mp_limb_t));
    work1 = calloc(limballocr, sizeof(mp_limb_t));
    partitionr = calloc(pr, sizeof(mp_limb_t *));
    partitioninvertr = calloc(pr, sizeof(mp_limb_t *));
    for (i=0; i<pr; i++) {
        partitionr[i] = calloc(limballocr, sizeof(mp_limb_t));
        partitioninvertr[i] = calloc(limballocr, sizeof(mp_limb_t));
    }
    
    /*Construction des p1 premiers ensembles*/
    setr = *partitionr;
    sets = *partition1;
    for (i=0; i<p1; i++) {
        /* partitionr[i] = ∑_{0 ≤ b ≤ n2} partition1[i]<<(N*b) */
        mpn_copyd(work1, sets, limbsize1);
        for (Nb=0; Nb<=n2*N; Nb+=N) {
            Nbmodbpl = Nb / mp_bits_per_limb;
            Nbdivbpl = Nb % mp_bits_per_limb;
            mpn_rshift(work0, work1, limbsize1+1, Nbdivbpl);
            mpn_ior_n(setr + Nbmodbpl, setr, work0, limbsize1+1);
        }
        setr++;
        sets++;
    }
    
    /*Construction des p2 autres ensembles*/
    sets = *partition2;
    for (i=0; i<p2; i++) {
        /* partitionr[i] = ∑_{b \in partition2[i]} work1<<(Nb-n1) avec work1 = 1…1 (n1 fois) */
        Nbmodbpl = n1 / mp_bits_per_limb;
        Nbdivbpl = n1 % mp_bits_per_limb;
        mpn_zero(work0, Nbmodbpl);
        mpn_neg(work1, work0, Nbmodbpl);
        work1[Nbmodbpl] = GMP_NUMB_MAX>>(mp_bits_per_limb - Nbdivbpl);
        for (<#initialization#>; <#condition#>; <#increment#>) {
            <#statements#>
        }
    }
    
}
