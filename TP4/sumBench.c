/* -------------------------------------------------------------------- */
/* sumBench.c        Kahan summation and other compensated sums versus  */
/*                   the regular summation: results and computing time  */
/*                   on huge vectors of random floating point numbers   */
/*                   of different ranges ("To go a little further",     */
/*                   questions 2 and 3)                                 */
/*                                                                      */
/*  0. Reference checks: MT19937 outputs, exact sum (12 test cases of   */
/*     CPython's math.fsum), and two tiny sums where Kahan fails.       */
/*  1. Double precision: 4 data sets of N values (U[0,1], 10^U[-10,10], */
/*     U[-1,1], +-10^U[-10,10]): relative error of every algorithm      */
/*     against the exact sum, and time per value (best of NB_RUNS).     */
/*  2. Single precision: the same data sets rounded to float.           */
/*  3. Ill-conditioned sums: relative error versus the condition number */
/*     cond = sum|x_i| / |sum x_i|, from 4 to 1e58.                     */
/*                                                                      */
/* Each data set has its own MT19937 stream (key of TP3 + data set      */
/* number), so the results are reproducible.                            */
/*                                                                      */
/* Build   : make sumBench sumBench_fast   (see Makefile)               */
/* Run     : ./sumBench            N = 1e7 values per data set          */
/*           ./sumBench 1e8        N = 1e8 (needs ~1.3 GB of memory)    */
/*           ./sumBench_fast       summation.c built with -ffast-math   */
/* -------------------------------------------------------------------- */
#define _POSIX_C_SOURCE 200809L        /* clock_gettime                */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#include "mt19937ar.h"
#include "summation.h"

#define DEFAULT_N     10000000      /* values per data set (80 MB)      */
#define MIN_N         1000
#define NB_RUNS       5             /* timing: best of NB_RUNS runs     */
#define LOG_MIN      -10.0          /* 10^U[LOG_MIN, LOG_MAX]           */
#define LOG_MAX       10.0
#define N_COND        1000000       /* values of the ill-conditioned    */
#define COND_E_STEP   10            /* sums: |big values| up to 2^E,    */
#define COND_E_MAX    200           /* E = 0, 10, ..., 200              */
#define NB_KEYS       5

typedef enum
{
   UNIFORM01, LOGUNIFORM, UNIFORM11, LOGUNIFORM_PM, NB_DATASETS
} dataset_t;

typedef double (*sumDouble_t)(const double x[], size_t n);
typedef double (*sumFloat_t) (const float  x[], size_t n);

typedef struct
{
   const char  * name;
   const char  * shortName;             /* column of section 3          */
   sumDouble_t   f;
} algoDouble_t;

typedef struct
{
   const char  * name;
   sumFloat_t    f;
} algoFloat_t;

static const char * const DATASET_NAMES[NB_DATASETS] = {
   "U[0,1]          same sign, 1 order of magnitude",
   "10^U[-10,10]    same sign, 20 orders of magnitude",
   "U[-1,1]         both signs (cancellations)",
   "+-10^U[-10,10]  both signs, 20 orders of magnitude" };

static double floatNaive   (const float x[], size_t n);
static double floatPairwise(const float x[], size_t n);
static double floatKahan   (const float x[], size_t n);

static const algoDouble_t ALGOS[] = {
   { "naive",             "naive",    sumNaive      },
   { "naive, 8 acc.",     "naive8",   sumNaive8     },
   { "pairwise",          "pairwise", sumPairwise   },
   { "Kahan",             "Kahan",    sumKahan      },
   { "Kahan, 8 acc.",     "Kahan8",   sumKahan8     },
   { "Neumaier",          "Neumaier", sumNeumaier   },
   { "Sum2 (TwoSum)",     "Sum2",     sumTwoSum     },
   { "Klein (2nd order)", "Klein",    sumKlein      },
   { "long double",       "long dbl", sumLongDouble },
   { "exact (fsum)",      "exact",    sumExact      } };

static const algoFloat_t ALGOS_F[] = {
   { "naive (float acc.)",    floatNaive    },
   { "pairwise (float acc.)", floatPairwise },
   { "Kahan (float acc.)",    floatKahan    },
   { "naive (double acc.)",   sumNaiveFD    } };

#define NB_ALGOS    ((int) (sizeof ALGOS   / sizeof ALGOS[0]))
#define NB_ALGOS_F  ((int) (sizeof ALGOS_F / sizeof ALGOS_F[0]))


/* -------------------------------------------------------------------- */
/* floatNaive, floatPairwise, floatKahan   float sums returned as       */
/*                   double (same signature as sumNaiveFD)              */
/*                                                                      */
/* Input : x  array of n floats                                         */
/* Output: the float result converted to double (exact conversion)      */
/* -------------------------------------------------------------------- */
static double floatNaive(const float x[], size_t n)
{
   return sumNaiveF(x, n);
}

static double floatPairwise(const float x[], size_t n)
{
   return sumPairwiseF(x, n);
}

static double floatKahan(const float x[], size_t n)
{
   return sumKahanF(x, n);
}


/* -------------------------------------------------------------------- */
/* wallTime          Elapsed time in seconds (monotonic clock)          */
/*                                                                      */
/* Input : none                                                         */
/* Output: current time                                                 */
/* -------------------------------------------------------------------- */
static double wallTime(void)
{
   struct timespec ts;

   clock_gettime(CLOCK_MONOTONIC, &ts);
   return (double) ts.tv_sec + 1e-9 * (double) ts.tv_nsec;
}


/* -------------------------------------------------------------------- */
/* seedStream        MT19937 stream number k: key of TP3 + k            */
/*                                                                      */
/* Input : g  state, k  stream number                                   */
/* Output: g initialized                                                */
/* -------------------------------------------------------------------- */
static void seedStream(mtState_t *g, uint32_t k)
{
   uint32_t key[NB_KEYS] = { 0x123, 0x234, 0x345, 0x456, k };

   init_by_array(g, key, NB_KEYS);
}


/* -------------------------------------------------------------------- */
/* generate          Fills x with n random values of a data set         */
/*                                                                      */
/* Input : x  array of n doubles, d  data set                           */
/* Output: x filled (53 random bits per value)                          */
/* -------------------------------------------------------------------- */
static void generate(double x[], size_t n, dataset_t d)
{
   mtState_t g;
   double    u;

   seedStream(&g, (uint32_t) d);
   for (size_t i = 0; i < n; i++)
   {
      u = genrand_res53(&g);
      switch (d)
      {
         case UNIFORM01:
            x[i] = u;
            break;
         case LOGUNIFORM:
            x[i] = pow(10.0, LOG_MIN + (LOG_MAX - LOG_MIN) * u);
            break;
         case UNIFORM11:
            x[i] = 2.0 * u - 1.0;
            break;
         default:                              /* LOGUNIFORM_PM        */
            x[i] = pow(10.0, LOG_MIN + (LOG_MAX - LOG_MIN) * u);
            if (genrand_int32(&g) & 1U)
               x[i] = -x[i];
            break;
      }
   }
}


/* -------------------------------------------------------------------- */
/* generateIllConditioned  Sum with a large condition number: half of   */
/*                   the values are pairs (v, -v) with |v| = U[1,2) 2^k,*/
/*                   k uniform in 0..e; the other half are U[0,1). The  */
/*                   pairs cancel exactly, so the exact sum is the sum  */
/*                   of the small values, while sum|x_i| ~ 2^e n.       */
/*                   The array is shuffled (Fisher-Yates).              */
/*                                                                      */
/* Input : x  array of n doubles, e  largest exponent                   */
/* Output: x filled                                                     */
/* -------------------------------------------------------------------- */
static void generateIllConditioned(double x[], size_t n, int e)
{
   mtState_t g;
   size_t    half = n / 2, i, j;
   double    v, t;

   seedStream(&g, 100U + (uint32_t) e);
   for (i = 0; i + 1 < half; i += 2)
   {
      int k = (int) (genrand_int32(&g) % (uint32_t) (e + 1));

      v = ldexp(1.0 + genrand_res53(&g), k);
      if (genrand_int32(&g) & 1U)
         v = -v;
      x[i]     = v;
      x[i + 1] = -v;
   }
   for (; i < n; i++)
      x[i] = genrand_res53(&g);

   for (i = n - 1; i > 0; i--)                 /* Fisher-Yates shuffle  */
   {
      j    = (size_t) (genrand_res53(&g) * (double) (i + 1));
      t    = x[i];
      x[i] = x[j];
      x[j] = t;
   }
}


/* -------------------------------------------------------------------- */
/* sumAbs            Sum of |x_i| (long double, enough for cond)        */
/*                                                                      */
/* Input : x  array of n doubles                                        */
/* Output: sum |x_i|                                                    */
/* -------------------------------------------------------------------- */
static double sumAbs(const double x[], size_t n)
{
   long double s = 0.0L;

   for (size_t i = 0; i < n; i++)
      s += fabs(x[i]);
   return (double) s;
}


/* -------------------------------------------------------------------- */
/* relError          |computed - exact| / |exact|                       */
/*                                                                      */
/* Input : computed, exact                                              */
/* Output: relative error (absolute error if exact is 0)                */
/* -------------------------------------------------------------------- */
static double relError(double computed, double exact)
{
   double err = fabs(computed - exact);

   return (exact != 0.0) ? err / fabs(exact) : err;
}


/* -------------------------------------------------------------------- */
/* checkMT           First outputs of init_by_array({0x123, 0x234,      */
/*                   0x345, 0x456}) = mt19937ar.out                     */
/*                                                                      */
/* Input : none                                                         */
/* Output: 1 if the generator is the reference one                      */
/* -------------------------------------------------------------------- */
static int checkMT(void)
{
   static const uint32_t REF[5] = { 1067595299U,  955945823U, 477289528U,
                                    4107218783U, 4228976476U };
   uint32_t              key[4] = { 0x123, 0x234, 0x345, 0x456 };
   mtState_t             g;
   int                   ok = 1;

   init_by_array(&g, key, 4);
   for (int i = 0; i < 5; i++)
      ok &= (genrand_int32(&g) == REF[i]);
   printf("  MT19937 init_by_array, 5 first outputs = mt19937ar.out : %s\n",
          ok ? "OK" : "FAILED");
   return ok;
}


/* -------------------------------------------------------------------- */
/* checkExact        sumExact on the test cases of CPython's test_fsum  */
/*                   (expected values checked with Python math.fsum)    */
/*                                                                      */
/* Input : none                                                         */
/* Output: 1 if every result is the expected one, bit for bit           */
/* -------------------------------------------------------------------- */
static int checkExact(void)
{
   static double v[1100];
   const double  p53 = 9007199254740992.0;              /* 2^53        */
   const double  e54 = ldexp(1.0, -54), e53 = ldexp(1.0, -53);
   const double  e100 = ldexp(1.0, -100);
   int           ok = 1, nbCases = 0, n;

#define CASE(expected)                                              \
   ok &= (sumExact(v, (size_t) n) == (expected)); nbCases++

   n = 0;                                                 CASE(0.0);
   n = 1; v[0] = 0.0;                                     CASE(0.0);
   n = 7; v[0] = 1e100; v[1] = 1.0; v[2] = -1e100; v[3] = 1e-100;
          v[4] = 1e50;  v[5] = -1.0; v[6] = -1e50;        CASE(1e-100);
   n = 3; v[0] = p53; v[1] = -0.5; v[2] = -e54;           CASE(p53 - 1.0);
   n = 3; v[0] = p53; v[1] = 1.0; v[2] = e100;            CASE(p53 + 2.0);
   n = 3; v[0] = p53 + 10.0; v[1] = 1.0; v[2] = e100;     CASE(p53 + 12.0);
   n = 3; v[0] = p53 - 4.0; v[1] = 0.5; v[2] = e54;       CASE(p53 - 3.0);
   n = 1000;
   for (int i = 1; i <= n; i++) v[i - 1] = 1.0 / i;
                                         CASE(0x1.df11f45f4e61ap+2);
   for (int i = 1; i <= n; i++) v[i - 1] = ((i % 2) ? -1.0 : 1.0) / i;
                                         CASE(-0x1.62a2af1bd3624p-1);
   n = 3; v[0] = 1e16; v[1] = 1.0; v[2] = 1e-16;   CASE(10000000000000002.0);
   n = 4; v[0] = 1e16 - 2.0; v[1] = 1.0 - e53;
          v[2] = -(1e16 - 2.0); v[3] = -(1.0 - e53);      CASE(0.0);
   n = 0;                          /* many partials: 2^(n+50) .. 2^n    */
   for (int k = -1074; k < 972; k += 2)
      v[n++] = ldexp(1.0, k) - ldexp(1.0, k + 50) + ldexp(1.0, k + 52);
   v[n++] = -ldexp(1.0, 1022);           CASE(0x1.5555555555555p+970);
#undef CASE

   printf("  exact sum (fsum), %d test cases of CPython test_fsum      : %s\n",
          nbCases, ok ? "OK" : "FAILED");
   return ok;
}


/* -------------------------------------------------------------------- */
/* tinyExamples      Two small sums computed by every algorithm         */
/*                                                                      */
/* Input : none                                                         */
/* Output: prints the results ([1, 1e100, 1, -1e100] = 2 exactly, and   */
/*         ten times 0.1 = 1.0000000000000000555, i.e. 1.0 rounded)     */
/* -------------------------------------------------------------------- */
static void tinyExamples(void)
{
   const double neumaierCase[4] = { 1.0, 1e100, 1.0, -1e100 };
   double       tenths[10];

   for (int i = 0; i < 10; i++)
      tenths[i] = 0.1;

   printf("  %-20s %24s %24s\n", "algorithm", "1 + 1e100 + 1 - 1e100",
          "0.1 + ... + 0.1 (x10)");
   for (int a = 0; a < NB_ALGOS; a++)
      printf("  %-20s %24.17g %24.17g\n", ALGOS[a].name,
             ALGOS[a].f(neumaierCase, 4), ALGOS[a].f(tenths, 10));
}


/* -------------------------------------------------------------------- */
/* benchDouble       Section 1: every algorithm on the 4 data sets      */
/*                                                                      */
/* Input : x  work array of n doubles                                   */
/* Output: prints, per data set, the exact sum, cond, then for every    */
/*         algorithm the result, its relative error and its time        */
/* -------------------------------------------------------------------- */
static void benchDouble(double x[], size_t n)
{
   for (int d = 0; d < NB_DATASETS; d++)
   {
      generate(x, n, (dataset_t) d);
      double exact = sumExact(x, n);
      double cond  = sumAbs(x, n) / fabs(exact);

      printf("\n--- %s ---\n", DATASET_NAMES[d]);
      printf("exact sum = %.17g   cond = sum|x|/|sum x| = %.3g\n",
             exact, cond);
      printf("%-20s %26s %10s %10s\n", "algorithm", "computed sum",
             "rel. error", "ns/value");

      for (int a = 0; a < NB_ALGOS; a++)
      {
         double s = 0.0, best = 1e300, t;

         for (int r = 0; r < NB_RUNS; r++)
         {
            t = wallTime();
            s = ALGOS[a].f(x, n);
            t = wallTime() - t;
            if (t < best) best = t;
         }
         printf("%-20s %26.17g %10.1e %10.2f\n", ALGOS[a].name, s,
                relError(s, exact), 1e9 * best / (double) n);
      }
   }
}


/* -------------------------------------------------------------------- */
/* benchFloat        Section 2: the data sets rounded to float          */
/*                                                                      */
/* Input : x   work array of n doubles (receives the floats, exactly)   */
/*         xf  work array of n floats                                   */
/* Output: prints, per data set, the exact sum of the floats and, for   */
/*         every float algorithm, result, relative error and time       */
/* -------------------------------------------------------------------- */
static void benchFloat(double x[], float xf[], size_t n)
{
   for (int d = 0; d < NB_DATASETS; d++)
   {
      generate(x, n, (dataset_t) d);
      for (size_t i = 0; i < n; i++)
      {
         xf[i] = (float) x[i];
         x[i]  = (double) xf[i];               /* exact conversion     */
      }
      double exact = sumExact(x, n);           /* exact sum of floats  */

      printf("\n--- %s, float ---\n", DATASET_NAMES[d]);
      printf("exact sum = %.17g\n", exact);
      printf("%-22s %24s %10s %10s\n", "algorithm", "computed sum",
             "rel. error", "ns/value");

      for (int a = 0; a < NB_ALGOS_F; a++)
      {
         double s = 0.0, best = 1e300, t;

         for (int r = 0; r < NB_RUNS; r++)
         {
            t = wallTime();
            s = ALGOS_F[a].f(xf, n);
            t = wallTime() - t;
            if (t < best) best = t;
         }
         printf("%-22s %24.17g %10.1e %10.2f\n", ALGOS_F[a].name, s,
                relError(s, exact), 1e9 * best / (double) n);
      }
   }
}


/* -------------------------------------------------------------------- */
/* benchCondition    Section 3: relative error versus condition number  */
/*                                                                      */
/* Input : x  work array of at least N_COND doubles                     */
/* Output: one line per E: cond, then the relative error of every       */
/*         algorithm (1.0e+00 or more: no correct digit at all)         */
/* -------------------------------------------------------------------- */
static void benchCondition(double x[])
{
   printf("%4s %9s", "E", "cond");
   for (int a = 0; a < NB_ALGOS; a++)
      printf(" %9s", ALGOS[a].shortName);
   printf("\n");

   for (int e = 0; e <= COND_E_MAX; e += COND_E_STEP)
   {
      generateIllConditioned(x, N_COND, e);
      double exact = sumExact(x, N_COND);

      printf("%4d %9.1e", e, sumAbs(x, N_COND) / fabs(exact));
      for (int a = 0; a < NB_ALGOS; a++)
         printf(" %9.1e", relError(ALGOS[a].f(x, N_COND), exact));
      printf("\n");
   }
}


int main(int argc, char *argv[])
{
   size_t   n = DEFAULT_N;
   double * x;
   float  * xf;

   if (argc > 1)
      n = (size_t) strtod(argv[1], NULL);      /* 1e7 accepted          */
   if (n < MIN_N)
   {
      fprintf(stderr, "Usage: %s [N >= %d]   (default 1e7)\n", argv[0],
              MIN_N);
      return EXIT_FAILURE;
   }

   x  = malloc((n > N_COND ? n : N_COND) * sizeof *x);
   xf = malloc(n * sizeof *xf);
   if (x == NULL || xf == NULL)
   {
      perror("malloc");
      return EXIT_FAILURE;
   }

   printf("=== 0. Reference checks ===\n");
   printf("  summation.c compiled: %s\n", summationBuild());
   if (!(checkMT() & checkExact()))
   {
      fprintf(stderr, "Reference check FAILED.\n");
      return EXIT_FAILURE;
   }
   printf("\n  Two tiny sums (exact: 2 and 1):\n");
   tinyExamples();

   printf("\n=== 1. Double precision, N = %zu values per data set, "
          "time = best of %d runs ===\n", n, NB_RUNS);
   benchDouble(x, n);

   printf("\n=== 2. Single precision (float), N = %zu ===\n", n);
   benchFloat(x, xf, n);

   printf("\n=== 3. Relative error versus condition number, N = %d, "
          "|big values| up to 2^E ===\n", N_COND);
   benchCondition(x);

   free(x);
   free(xf);
   return EXIT_SUCCESS;
}
