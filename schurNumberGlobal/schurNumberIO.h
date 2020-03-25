//
//  schurNumberIO.h
//  schurNumberRecursive
//
//  Created by rubis on 28/01/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#ifndef schurNumberIO_h
#define schurNumberIO_h

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

unsigned long schurNumberGetSetMaximum(char *str);
void schurNumberGetSet(char *str, mp_limb_t *set, mp_limb_t *setinvert, mp_size_t limballoc);
unsigned long schurNumberGetPartition(unsigned long p, char **str, mp_limb_t **partition, mp_limb_t **partitioninvert, mp_size_t limballoc);

void schurNumberPrintSet(unsigned long n, mp_limb_t *set);
void schurNumberPrintPartition(unsigned long p, unsigned long n, mp_limb_t **partition);

#endif /* schurNumberIO_h */
