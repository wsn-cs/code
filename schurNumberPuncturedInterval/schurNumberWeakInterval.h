//
//  schurNumberWeakInterval.h
//  schurNumberPuncturedInterval
//
//  Created by rubis on 16/03/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#ifndef schurNumberWeakInterval_h
#define schurNumberWeakInterval_h

#include "../schurNumberGlobal/schurNumberIO.h"
#include "../schurNumberGlobal/schurNumberPartitionStruc.h"

unsigned long schurNumberPuncturedInterval(schur_number_partition_t *partitionstruc, struct schurNumberIOAction *action, unsigned long nlimit, mp_limb_t **constraint_partition, mp_size_t constraint_size);

#endif /* schurNumberWeakInterval_h */
