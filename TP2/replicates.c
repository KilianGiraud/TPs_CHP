/* ------------------------------------------------------------------ */
/* replicates.c - n independent replicates of the Monte Carlo pi      */
/*                estimation, run sequentially.                        */
/*                                                                    */
/* Generator : Mersenne Twister MT19937 (Matsumoto & Nishimura),       */
/*             initialized ONCE with the original reference key        */
/*             {0x123, 0x234, 0x345, 0x456} (mt19937ar.c).             */
/* Independence between replicates is obtained by NOT re-seeding the   */
/* generator between two experiments: each replicate consumes the next */
/* consecutive part of the same long random sequence.                  */
/*                                                                    */
/* Compile : gcc -O2 -Wall -o replicates replicates.c -lm              */
/* Run     : ./replicates                 (10 replicates, 3 sizes)     */
/*           ./replicates 10 1000000      (nbReplicates nbPoints)      */
/* ------------------------------------------------------------------ */

/* Make math.h expose M_PI (POSIX / MSVC) so it comes from <math.h>    */
#define _DEFAULT_SOURCE
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* ================================================================== */
/*  Mersenne Twister MT19937                                           */
/* ================================================================== */
#define MT_N          624
#define MT_M          397
#define MATRIX_A      0x9908b0dfUL
#define UPPER_MASK    0x80000000UL
#define LOWER_MASK    0x7fffffffUL

static unsigned long mt[MT_N];
static int           mti = MT_N + 1;

static void init_genrand(unsigned long s)
{
    mt[0] = s & 0xffffffffUL;
    for (mti = 1; mti < MT_N; mti++)
    {
        mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
        mt[mti] &= 0xffffffffUL;
    }
}

static void init_by_array(unsigned long key[], int keyLength)
{
    int i = 1, j = 0, k;

    init_genrand(19650218UL);
    k = (MT_N > keyLength ? MT_N : keyLength);
    for (; k; k--)
    {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL))
                + key[j] + j;
        mt[i] &= 0xffffffffUL;
        i++; j++;
        if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i = 1; }
        if (j >= keyLength) j = 0;
    }
    for (k = MT_N - 1; k; k--)
    {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL)) - i;
        mt[i] &= 0xffffffffUL;
        i++;
        if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i = 1; }
    }
    mt[0] = 0x80000000UL;   /* MSB is 1: non-zero initial array */
}

static unsigned long genrand_int32(void)
{
    static const unsigned long mag01[2] = { 0x0UL, MATRIX_A };
    unsigned long y;
    int kk;

    if (mti >= MT_N)        /* regenerate the whole block of 624 words */
    {
        if (mti == MT_N + 1) init_genrand(5489UL);

        for (kk = 0; kk < MT_N - MT_M; kk++)
        {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; kk < MT_N - 1; kk++)
        {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (MT_M - MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[MT_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        mti = 0;
    }

    y = mt[mti++];
    /* tempering */
    y ^= (y >> 11);
    y ^= (y << 7)  & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return y;
}

/* Uniform double in [0, 1] */
static inline double genrand_real1(void)
{
    return genrand_int32() * (1.0 / 4294967295.0);
}

/* ================================================================== */
/*  Monte Carlo                                                        */
/* ================================================================== */

/* One experiment: estimate pi with nbPoints points (2 numbers each)   */
double simPi(long long nbPoints)
{
    long long inside = 0;
    double xr, yr;

    for (long long i = 0; i < nbPoints; i++)
    {
        xr = genrand_real1();
        yr = genrand_real1();
        if (xr * xr + yr * yr <= 1.0)
            inside++;
    }
    return 4.0 * (double) inside / (double) nbPoints;
}

/* Run nbReplicates independent experiments, store them, return mean   */
double runReplicates(int nbReplicates, long long nbPoints, double results[])
{
    double sum = 0.0;

    for (int r = 0; r < nbReplicates; r++)
    {
        /* NO re-initialization of the generator here: independence    */
        results[r] = simPi(nbPoints);
        sum += results[r];
    }
    return sum / nbReplicates;
}

static void experiment(int nbReplicates, long long nbPoints)
{
    double *results = malloc(nbReplicates * sizeof(double));
    if (results == NULL) { perror("malloc"); exit(EXIT_FAILURE); }

    clock_t start  = clock();
    double  meanPi = runReplicates(nbReplicates, nbPoints, results);
    double  secs   = (double)(clock() - start) / CLOCKS_PER_SEC;

    printf("\n=== %d replicates x %lld points ===\n", nbReplicates, nbPoints);
    for (int r = 0; r < nbReplicates; r++)
        printf("  replicate %2d : %.8f\n", r + 1, results[r]);

    double absErr = fabs(meanPi - M_PI);
    double relErr = absErr / M_PI;

    printf("  meanPi         = %.10f\n", meanPi);
    printf("  M_PI           = %.10f\n", M_PI);
    printf("  absolute error = %.3e\n", absErr);
    printf("  relative error = %.3e  (%.5f %%)\n", relErr, 100.0 * relErr);
    printf("  time           = %.2f s\n", secs);

    free(results);
}

int main(int argc, char *argv[])
{
    /* Original initialization from Makoto Matsumoto's reference code */
    unsigned long init[4] = { 0x123, 0x234, 0x345, 0x456 };
    init_by_array(init, 4);

    if (argc == 3)
    {
        experiment(atoi(argv[1]), atoll(argv[2]));
        return 0;
    }

    experiment(10, 1000LL);
    experiment(10, 1000000LL);
    experiment(10, 1000000000LL);
    return 0;
}