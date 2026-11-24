/* -------------------------------------------------------------------- */
/* mt19937ar.c       Mersenne Twister MT19937 (Matsumoto & Nishimura),  */
/*                   reentrant version of mt19937ar.c (2002/1/26)       */
/*                                                                      */
/* Same algorithm and outputs as the reference code (checked at the     */
/* start of sumBench: the first outputs of init_by_array({0x123,        */
/* 0x234, 0x345, 0x456}) are those of mt19937ar.out).                   */
/* -------------------------------------------------------------------- */
#include "mt19937ar.h"


/* -------------------------------------------------------------------- */
/* init_genrand      Initializes the state with a 32-bit seed           */
/*                                                                      */
/* Input : g  state, s  seed                                            */
/* Output: g initialized                                                */
/* -------------------------------------------------------------------- */
void init_genrand(mtState_t *g, uint32_t s)
{
   uint32_t * mt = g->mt;

   mt[0] = s;
   for (int i = 1; i < MT_N; i++)
      mt[i] = 1812433253U * (mt[i - 1] ^ (mt[i - 1] >> 30)) + (uint32_t) i;
   g->mti = MT_N;
}


/* -------------------------------------------------------------------- */
/* init_by_array     Initializes the state with an array of 32-bit keys */
/*                   (recommended initialization)                       */
/*                                                                      */
/* Input : g  state, key  array of keyLength words                      */
/* Output: g initialized                                                */
/* -------------------------------------------------------------------- */
void init_by_array(mtState_t *g, const uint32_t key[], int keyLength)
{
   uint32_t * mt = g->mt;
   int        i  = 1,
              j  = 0,
              k;

   init_genrand(g, 19650218U);
   k = (MT_N > keyLength ? MT_N : keyLength);
   for (; k; k--)
   {
      mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525U))
              + key[j] + (uint32_t) j;                     /* non linear */
      i++; j++;
      if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i = 1; }
      if (j >= keyLength) j = 0;
   }
   for (k = MT_N - 1; k; k--)
   {
      mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941U))
              - (uint32_t) i;                              /* non linear */
      i++;
      if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i = 1; }
   }
   mt[0] = 0x80000000U;     /* MSB is 1: assures a non-zero initial array */
   g->mti = MT_N;
}


/* -------------------------------------------------------------------- */
/* genrand_int32     Random 32-bit integer in [0, 0xffffffff]           */
/*                                                                      */
/* Input : g  initialized state                                         */
/* Output: next output of the generator                                 */
/* -------------------------------------------------------------------- */
uint32_t genrand_int32(mtState_t *g)
{
   static const uint32_t   mag01[2] = { 0x0U, MATRIX_A };
   uint32_t              * mt = g->mt;
   uint32_t                y;
   int                     kk;

   if (g->mti >= MT_N)        /* regenerate the whole block of 624 words */
   {
      for (kk = 0; kk < MT_N - MT_M; kk++)
      {
         y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
         mt[kk] = mt[kk + MT_M] ^ (y >> 1) ^ mag01[y & 0x1U];
      }
      for (; kk < MT_N - 1; kk++)
      {
         y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
         mt[kk] = mt[kk + (MT_M - MT_N)] ^ (y >> 1) ^ mag01[y & 0x1U];
      }
      y = (mt[MT_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1U];
      g->mti = 0;
   }

   y = mt[g->mti++];
   y ^= (y >> 11);                                         /* tempering */
   y ^= (y << 7)  & 0x9d2c5680U;
   y ^= (y << 15) & 0xefc60000U;
   y ^= (y >> 18);
   return y;
}


/* -------------------------------------------------------------------- */
/* genrand_res53     Uniform double in [0, 1) with 53-bit resolution    */
/*                   (function of the reference code)                   */
/*                                                                      */
/* Input : g  initialized state                                         */
/* Output: (a * 2^26 + b) / 2^53 with a, b of 27 and 26 random bits     */
/* -------------------------------------------------------------------- */
double genrand_res53(mtState_t *g)
{
   uint32_t a = genrand_int32(g) >> 5,
            b = genrand_int32(g) >> 6;

   return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
}
