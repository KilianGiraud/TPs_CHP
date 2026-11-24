/* -------------------------------------------------------------------- */
/* mt19937ar.h       Mersenne Twister MT19937 (Matsumoto & Nishimura),  */
/*                   reentrant version of the reference code            */
/*                   mt19937ar.c (2002/1/26), taken from TP3            */
/*                                                                      */
/* Declarations only (ISIMA style guide): the code is in mt19937ar.c.   */
/* The generator state is a structure given as parameter, so that each  */
/* data set / thread can own its generator.                             */
/*                                                                      */
/* genrand_res53 is used to draw the doubles of the TP: 53 random bits, */
/* the full precision of a double. With genrand_real1 (32 random bits   */
/* only), the values would be multiples of 2^-32 and many sums would be */
/* exact by chance, hiding the rounding errors studied here.            */
/* -------------------------------------------------------------------- */
#ifndef MT19937AR_H
#define MT19937AR_H

#include <stdint.h>

#define MT_N        624
#define MT_M        397
#define MATRIX_A    0x9908b0dfU      /* constant vector a               */
#define UPPER_MASK  0x80000000U      /* most significant w-r bits       */
#define LOWER_MASK  0x7fffffffU      /* least significant r bits        */

typedef struct
{
   uint32_t mt[MT_N];                /* state vector (624 words)         */
   int      mti;                     /* next word to temper              */
} mtState_t;

void     init_genrand (mtState_t * g, uint32_t s);
void     init_by_array(mtState_t * g, const uint32_t key[], int keyLength);
uint32_t genrand_int32(mtState_t * g);
double   genrand_res53(mtState_t * g);

#endif /* MT19937AR_H */
