//
//  schurNumberMethods.h
//  SchurNumber
//
//  Created by rubis on 10/02/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

/* Ce fichier déclare les méthodes pouvant être employées.
Elles sont construites selon le schéma suivant:
  unsigned long func(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit).
*/

#ifndef schurNumberMethods_h
#define schurNumberMethods_h

#include "../schurNumberGlobal/schurNumberPartitionStruc.h"
#include "../schurNumberGlobal/schurNumberIOAction.h"

mp_bitcnt_t mpn_rscan1 (mp_srcptr up, mp_bitcnt_t starting_bit);

unsigned long schur_number_weak_exhaustive(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit);

unsigned long schur_number_exhaustive(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit);

unsigned long schur_number_weak_branch_bound(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit);

unsigned long schur_number_branch_bound(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit);

unsigned long schur_number_stacked_branch_bound(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit);

unsigned long schur_number_weak_stacked_branch_bound(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit);

#endif /* schurNumberMethods_h */
