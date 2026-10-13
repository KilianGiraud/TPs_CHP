/* ------------------------------------------------------------------ */
/* simPI.c - One Monte Carlo replication of PI, seed given in input.   */
/*                                                                    */
/* SPMD model: the SAME program is launched N times in parallel by the */
/* shell, only the data (= the random stream) changes, through the     */
/* seed given on the command line.                                     */
/*                                                                    */
/* Usage   : ./simPI <seed> [nbPoints]     (default 1 000 000 000)     */
/*           ./simPI 19937 > PI1                                       */
/* Output  : only the estimated PI on stdout (easy to post-process)    */
/*                                                                    */
/* Compile : gcc -O2 -Wall -o simPI simPI.c                            */
/* ------------------------------------------------------------------ */
#include <stdio.h>
#include <stdlib.h>

/* ================================================================== */
/*  Mersenne Twister MT19937 (Matsumoto & Nishimura)                   */
/* ================================================================== */
#define MT_N        624
#define MT_M        397
#define MATRIX_A    0x9908b0dfUL
#define UPPER_MASK  0x80000000UL
#define LOWER_MASK  0x7fffffffUL

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

static unsigned long genrand_int32(void)
{
    static const unsigned long mag01[2] = { 0x0UL, MATRIX_A };
    unsigned long y;
    int kk;

    if (mti >= MT_N)
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

/* ------------------------------------------------------------------ */
/* simulPi    Calcul de PI par Monte Carlo                             */
/*                                                                    */
/* On tire au hasard des points dans le carre unite [0..1]x[0..1]      */
/* On compte tous les points qui sont sur le quart de disque unite     */
/* La proportion de ces points converge vers PI/4.                     */
/*                                                                    */
/* En entree : inNbPts le nombre de points pour l'estimation de PI     */
/* En sortie : une valeur approchee de PI en double precision          */
/* ------------------------------------------------------------------ */
double simulPi(long long inNbPts)
{
    long long i;
    double    xr     = 0.;
    double    yr     = 0.;
    double    pi     = 0.;
    long long inDisk = 0;

    for (i = 0; i < inNbPts; i++)
    {
        xr = genrand_real1();
        yr = genrand_real1();

        if ((xr * xr + yr * yr) <= 1)
        {
            inDisk++;
        }
    }

    /* Transtypage des entiers en double pour le calcul */
    pi = ((double) inDisk / (double) inNbPts) * 4.;

    return pi;
}

int main(int argc, char *argv[])
{
    unsigned long seed;
    long long     nbPts = 1000000000LL;

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s <seed> [nbPoints]\n", argv[0]);
        return EXIT_FAILURE;
    }

    seed = strtoul(argv[1], NULL, 10);
    if (argc >= 3)
        nbPts = atoll(argv[2]);

    init_genrand(seed);                 /* one stream per seed */

    printf("%.10f\n", simulPi(nbPts));  /* result on stdout only */
    return EXIT_SUCCESS;
}