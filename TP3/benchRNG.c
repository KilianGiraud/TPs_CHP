/* ------------------------------------------------------------------ */
/* benchRNG.c - Q1: check and time the two "fine" generators          */
/*                                                                    */
/*  1. Reference checks: the outputs of our MT19937 and MRG32k3a are  */
/*     compared with known reference values (original codes).         */
/*  2. Performance: nbCalls calls (default 10^9) to the functions     */
/*     drawing a pseudo-random number between 0 and 1:                */
/*       - rand() / RAND_MAX   (reference only, used in TP2 Q1)       */
/*       - genrand_real1()     Mersenne Twister MT19937, [0, 1]       */
/*       - MRG32k3a()          L'Ecuyer, ]0, 1[                       */
/*     The draws are summed so that the compiler cannot remove the    */
/*     loop; the mean must be close to 0.5.                           */
/*                                                                    */
/* Compile : gcc -O0 -Wall -o benchRNG_O0 benchRNG.c                  */
/*           gcc -O2 -Wall -o benchRNG_O2 benchRNG.c                  */
/*           gcc -O3 -march=native -Wall -o benchRNG_O3 benchRNG.c    */
/* Run     : ./benchRNG_O2          (10^9 calls per generator)        */
/*           ./benchRNG_O2 1e8      (nbCalls per generator)           */
/* ------------------------------------------------------------------ */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include "mt19937ar.h"
#include "mrg32k3a.h"

/* ================================================================== */
/*  Reference checks                                                  */
/* ================================================================== */

/* MT19937: first outputs of mt19937ar.out + C++11 std::mt19937 value */
static int checkMT(void)
{
    static const uint32_t ref[5] = { 1067595299U,  955945823U,  477289528U,
                                     4107218783U, 4228976476U };
    uint32_t key[4] = { 0x123, 0x234, 0x345, 0x456 };
    mt_state g;
    uint32_t y = 0;
    int      ok1 = 1, ok2;

    init_by_array(&g, key, 4);
    for (int i = 0; i < 5; i++)
        ok1 &= (genrand_int32(&g) == ref[i]);

    init_genrand(&g, 5489U);
    for (int i = 0; i < 10000; i++)
        y = genrand_int32(&g);
    ok2 = (y == 4123659995U);

    printf("  MT19937  init_by_array, 5 first outputs = mt19937ar.out : %s\n",
           ok1 ? "OK" : "FAILED");
    printf("  MT19937  init_genrand(5489), 10000th output = %10u  : %s\n",
           y, ok2 ? "OK" : "FAILED");
    return ok1 && ok2;
}

/* MRG32k3a: first output and seed of the 2nd stream (RngStreams)    */
static int checkMRG(void)
{
    /* first output with seed {12345 x 6}: p1 - p2 = 545508589 exactly */
    static const double refFirst = 545508589.0 * MRG_NORM;
    /* default seed advanced by 2^127 steps (2nd stream of RngStreams) */
    static const double refStream1[6] = { 3692455944.0, 1366884236.0,
                                          2968912127.0,  335948734.0,
                                          4161675175.0,  475798818.0 };
    mrg_state g;
    double    u;
    int       ok1, ok2;

    mrg_init_stream(&g, 0);
    u   = MRG32k3a(&g);
    ok1 = (u == refFirst);

    mrg_init_stream(&g, 1);
    ok2 = g.s10 == refStream1[0] && g.s11 == refStream1[1]
       && g.s12 == refStream1[2] && g.s20 == refStream1[3]
       && g.s21 == refStream1[4] && g.s22 == refStream1[5];

    printf("  MRG32k3a seed {12345 x 6}, first output = %.10f     : %s\n",
           u, ok1 ? "OK" : "FAILED");
    printf("  MRG32k3a jump 2^127, seed of stream 1 = RngStreams      : %s\n",
           ok2 ? "OK" : "FAILED");
    return ok1 && ok2;
}

/* ================================================================== */
/*  Performance                                                       */
/* ================================================================== */

static void report(const char *name, long long n, double secs, double sum)
{
    printf("  %-20s : %6.2f s | %6.2f ns/call | %7.1f M calls/s | mean = %.6f\n",
           name, secs, 1e9 * secs / (double) n, (double) n / secs / 1e6,
           sum / (double) n);
}

static void benchRand(long long n)
{
    double  sum = 0.0;
    clock_t start;

    srand(12345);
    start = clock();
    for (long long i = 0; i < n; i++)
        sum += (double) rand() / (double) RAND_MAX;
    report("rand() (TP2)", n, (double)(clock() - start) / CLOCKS_PER_SEC, sum);
}

static void benchMT(long long n)
{
    uint32_t key[4] = { 0x123, 0x234, 0x345, 0x456 };   /* Makoto's key */
    mt_state g;
    double   sum = 0.0;
    clock_t  start;

    init_by_array(&g, key, 4);
    start = clock();
    for (long long i = 0; i < n; i++)
        sum += genrand_real1(&g);
    report("MT19937 real1", n, (double)(clock() - start) / CLOCKS_PER_SEC, sum);
}

static void benchMRG(long long n)
{
    mrg_state g;
    double    sum = 0.0;
    clock_t   start;

    mrg_init_stream(&g, 0);
    start = clock();
    for (long long i = 0; i < n; i++)
        sum += MRG32k3a(&g);
    report("MRG32k3a", n, (double)(clock() - start) / CLOCKS_PER_SEC, sum);
}

int main(int argc, char *argv[])
{
    /* strtod accepts both 1000000000 and 1e9 */
    long long nbCalls = (argc > 1) ? (long long) strtod(argv[1], NULL)
                                   : 1000000000LL;
    if (nbCalls <= 0)
    {
        fprintf(stderr, "Usage: %s [nbCalls]\n", argv[0]);
        return EXIT_FAILURE;
    }

    printf("=== %s : reference checks ===\n", argv[0]);
    if (!(checkMT() & checkMRG()))
    {
        fprintf(stderr, "Reference check FAILED: generator not trusted.\n");
        return EXIT_FAILURE;
    }

    printf("\n=== %s : %lld calls per generator ===\n", argv[0], nbCalls);
    benchRand(nbCalls);
    benchMT(nbCalls);
    benchMRG(nbCalls);
    return EXIT_SUCCESS;
}
