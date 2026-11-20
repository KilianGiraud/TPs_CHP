/* ------------------------------------------------------------------ */
/* mt19937ar.h - Mersenne Twister MT19937 (Matsumoto & Nishimura)     */
/*                                                                    */
/* Reentrant version of the reference code mt19937ar.c (2002/1/26):   */
/* same algorithm, same outputs, but the generator state is stored in */
/* a structure given as parameter instead of global variables, so     */
/* that each OpenMP thread / replicate can own its own generator.     */
/*                                                                    */
/* Checked against reference outputs in benchRNG.c:                   */
/*  - init_by_array({0x123, 0x234, 0x345, 0x456}) gives the first     */
/*    values of mt19937ar.out (1067595299 955945823 ...)              */
/*  - init_genrand(5489): 10000th output = 4123659995 (value required */
/*    by the C++11 standard for std::mt19937)                         */
/*                                                                    */
/* A state MUST be initialized (init_genrand or init_by_array)        */
/* before the first call to genrand_int32 / genrand_real1.            */
/* ------------------------------------------------------------------ */
#ifndef MT19937AR_H
#define MT19937AR_H

#include <stdint.h>

#define MT_N        624
#define MT_M        397
#define MATRIX_A    0x9908b0dfU     /* constant vector a           */
#define UPPER_MASK  0x80000000U     /* most significant w-r bits   */
#define LOWER_MASK  0x7fffffffU     /* least significant r bits    */

typedef struct
{
    uint32_t mt[MT_N];              /* state vector (624 words)     */
    int      mti;                   /* next word to temper          */
} mt_state;

/* Initialize the state with a 32-bit seed                            */
static inline void init_genrand(mt_state *g, uint32_t s)
{
    uint32_t *mt = g->mt;

    mt[0] = s;
    for (int i = 1; i < MT_N; i++)
        mt[i] = 1812433253U * (mt[i - 1] ^ (mt[i - 1] >> 30)) + (uint32_t) i;
    g->mti = MT_N;
}

/* Initialize the state with an array of 32-bit keys (recommended)    */
static inline void init_by_array(mt_state *g, const uint32_t key[], int keyLength)
{
    uint32_t *mt = g->mt;
    int i = 1, j = 0, k;

    init_genrand(g, 19650218U);
    k = (MT_N > keyLength ? MT_N : keyLength);
    for (; k; k--)
    {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525U))
                + key[j] + (uint32_t) j;                    /* non linear */
        i++; j++;
        if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i = 1; }
        if (j >= keyLength) j = 0;
    }
    for (k = MT_N - 1; k; k--)
    {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941U))
                - (uint32_t) i;                             /* non linear */
        i++;
        if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i = 1; }
    }
    mt[0] = 0x80000000U;    /* MSB is 1: assures a non-zero initial array */
    g->mti = MT_N;
}

/* Random 32-bit integer in [0, 0xffffffff]                           */
static inline uint32_t genrand_int32(mt_state *g)
{
    static const uint32_t mag01[2] = { 0x0U, MATRIX_A };
    uint32_t *mt = g->mt;
    uint32_t  y;
    int       kk;

    if (g->mti >= MT_N)     /* regenerate the whole block of 624 words */
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
    /* tempering */
    y ^= (y >> 11);
    y ^= (y << 7)  & 0x9d2c5680U;
    y ^= (y << 15) & 0xefc60000U;
    y ^= (y >> 18);
    return y;
}

/* Uniform double in [0, 1]                                           */
static inline double genrand_real1(mt_state *g)
{
    return genrand_int32(g) * (1.0 / 4294967295.0);
}

#endif /* MT19937AR_H */
