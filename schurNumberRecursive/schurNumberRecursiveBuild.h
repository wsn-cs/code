//
//  schurNumberRecursiveBuild.h
//  schurNumberRecursive
//
//  Created by rubis on 01/12/2019.
//  Copyright © 2019 rubis. All rights reserved.
//

#ifndef schurNumberRecursiveBuild_h
#define schurNumberRecursiveBuild_h

#include <stdlib.h>
#include <gmp.h>

void printPartition(unsigned long p, unsigned long n, mp_limb_t **partition);

char schurNumberRecursiveIteration(unsigned long *integer_stack, size_t stack_size, unsigned int p, mp_limb_t **partition, mp_limb_t **reverse_partition);

char schurNumberRecursiveBuild(unsigned long n, unsigned int p);

#endif /* schurNumberRecursiveBuild_h */
