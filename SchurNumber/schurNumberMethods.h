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
#include "../schurNumberGlobal/schurNumberIO.h"

mp_bitcnt_t mpn_rscan1 (mp_srcptr up, mp_bitcnt_t starting_bit);

unsigned long schurNumberWeakExhaustive(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit);
//unsigned long schurNumberWeakExhaustive2(schur_number_partition_t *partitionstruc, schur_number_action_t *action);

unsigned long schurNumberExhaustive(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit);

unsigned long schurNumberWeakBranchBound(schur_number_partition_t *partitionstruc, schur_number_action_t *action, unsigned long nlimit);

#endif /* schurNumberMethods_h */
