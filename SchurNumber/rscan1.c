//
//  rscan1.c
//  SchurNumber
//
//  Created by Gabriel Merlin on 03/03/2019.
//  mpn_rscan1 -- Scan decreasingly from a given bit position for the last set bit.
//

#ifndef MPN_RSCAN1

#define MPN_RSCAN1

#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"

mp_bitcnt_t mpn_rscan1 (mp_srcptr up, mp_bitcnt_t starting_bit)
{
    mp_size_t starting_word;
    mp_limb_t alimb;
    int cnt;
    mp_srcptr p;
    
    /* Start at the word implied by STARTING_BIT. */
    starting_word = starting_bit / GMP_NUMB_BITS;
    p = up + starting_word;
    alimb = *p;
    
    /* Mask off any bits after STARTING_BIT in the first limb. */
    alimb &= - (mp_limb_t) 1 >> (GMP_NUMB_BITS - starting_bit % GMP_NUMB_BITS -1);
    
    while (alimb == 0) {
        p--;
        alimb = *p;
    }
    
    count_leading_zeros(cnt, alimb);
    return (p - up + 1) * GMP_NUMB_BITS - cnt;
}

#endif
