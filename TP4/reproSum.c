/* -------------------------------------------------------------------- */
/* reproSum.c        Numerical reproducibility of a sum when the order  */
/*                   of the additions changes: sequential orders and    */
/*                   parallel reductions with OpenMP ("Out of order     */
/*                   execution" and "reduction phase" of the subject)   */
/*                                                                      */
/*  1. The same N values summed in 6 orders (forward, backward, sorted  */
/*     by increasing / decreasing |x|, 2 random permutations) with the  */
/*     naive, Kahan, Neumaier and exact sums.                           */
/*  2. OpenMP reduction with 1 .. 12 threads:                           */
/*       - naive   : #pragma omp parallel for reduction(+:s)            */
/*       - Kahan   : one Kahan sum per thread on a contiguous block,    */
/*                   the (s, c) pairs combined exactly in thread order  */
/*       - exact   : exact partials (fsumAdd) per thread, merged        */
/*                   exactly: same result for any number of threads     */
/*     with the result, its difference to the exact sum and the time.   */
/*  3. Same number of threads (12), NB_REPEAT runs: the reduction(+)    */
/*     combines the partial sums of the threads in the order they       */
/*     finish (and with schedule(dynamic) the blocks themselves change) */
/*     so the SAME binary on the SAME data gives different results.     */
/*                                                                      */
/* Build   : make reproSum                                              */
/* Run     : ./reproSum            N = 1e7 values U[-1,1]               */
/*           ./reproSum 1e8                                             */
/* -------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "mt19937ar.h"
#include "summation.h"

#define DEFAULT_N     10000000      /* values (80 MB)                   */
#define MIN_N         1000
#define MAX_THREADS   64
#define NB_RUNS       5             /* timing: best of NB_RUNS          */
#define NB_REPEAT     10            /* section 3: runs with 12 threads  */
#define DYN_CHUNK     10000         /* schedule(dynamic, DYN_CHUNK)     */
#define NB_ORDERS     6
#define NB_SEQ_ALGOS  4

typedef double (*sumDouble_t)(const double x[], size_t n);

static const int THREADS[] = { 1, 2, 3, 4, 6, 8, 12 };
#define NB_THREADS ((int) (sizeof THREADS / sizeof THREADS[0]))


/* -------------------------------------------------------------------- */
/* compareAbsUp, compareAbsDown   qsort comparators on |x|              */
/*                                                                      */
/* Input : a, b  pointers to doubles                                    */
/* Output: <0, 0, >0 for increasing (Up) or decreasing (Down) |x|       */
/* -------------------------------------------------------------------- */
static int compareAbsUp(const void *a, const void *b)
{
   double x = fabs(*(const double *) a), y = fabs(*(const double *) b);

   return (x > y) - (x < y);
}

static int compareAbsDown(const void *a, const void *b)
{
   return compareAbsUp(b, a);
}


/* -------------------------------------------------------------------- */
/* shuffle           Fisher-Yates random permutation                    */
/*                                                                      */
/* Input : x  array of n doubles, g  generator                          */
/* Output: x permuted                                                   */
/* -------------------------------------------------------------------- */
static void shuffle(double x[], size_t n, mtState_t *g)
{
   for (size_t i = n - 1; i > 0; i--)
   {
      size_t j = (size_t) (genrand_res53(g) * (double) (i + 1));
      double t = x[i];

      x[i] = x[j];
      x[j] = t;
   }
}


/* -------------------------------------------------------------------- */
/* countDistinct     Number of different values in a small array        */
/*                                                                      */
/* Input : v  array of n doubles                                        */
/* Output: number of distinct values (bit for bit)                      */
/* -------------------------------------------------------------------- */
static int countDistinct(const double v[], int n)
{
   int nb = 0;

   for (int i = 0; i < n; i++)
   {
      int seen = 0;

      for (int j = 0; j < i; j++)
         if (v[j] == v[i]) seen = 1;
      nb += !seen;
   }
   return nb;
}


/* -------------------------------------------------------------------- */
/* sequentialOrders  Section 1: one data set, 6 orders, 4 algorithms    */
/*                                                                      */
/* Input : x     the N values in their original order                   */
/*         work  work array of n doubles                                */
/*         exact the exact sum                                          */
/* Output: prints s - exact for every order and algorithm               */
/* -------------------------------------------------------------------- */
static void sequentialOrders(const double x[], double work[], size_t n,
                             double exact)
{
   static const char * const ORDERS[NB_ORDERS] = {
      "forward", "backward", "increasing |x|", "decreasing |x|",
      "random permutation 1", "random permutation 2" };
   static const char * const NAMES[NB_SEQ_ALGOS] = {
      "naive", "Kahan", "Neumaier", "exact" };
   const sumDouble_t         ALGOS[NB_SEQ_ALGOS] = {
      sumNaive, sumKahan, sumNeumaier, sumExact };
   double                    res[NB_SEQ_ALGOS][NB_ORDERS];
   mtState_t                 g;
   uint32_t                  key[2] = { 0x123, 0x456 };

   init_by_array(&g, key, 2);
   for (int o = 0; o < NB_ORDERS; o++)
   {
      memcpy(work, x, n * sizeof *work);
      switch (o)
      {
         case 1:
            for (size_t i = 0; i < n; i++) work[i] = x[n - 1 - i];
            break;
         case 2: qsort(work, n, sizeof *work, compareAbsUp);   break;
         case 3: qsort(work, n, sizeof *work, compareAbsDown); break;
         case 4:
         case 5: shuffle(work, n, &g);                         break;
         default: break;
      }
      for (int a = 0; a < NB_SEQ_ALGOS; a++)
         res[a][o] = ALGOS[a](work, n);
   }

   printf("%-22s", "order  (sum - exact)");
   for (int a = 0; a < NB_SEQ_ALGOS; a++)
      printf(" %12s", NAMES[a]);
   printf("\n");
   for (int o = 0; o < NB_ORDERS; o++)
   {
      printf("%-22s", ORDERS[o]);
      for (int a = 0; a < NB_SEQ_ALGOS; a++)
         printf(" %12.3e", res[a][o] - exact);
      printf("\n");
   }
   printf("%-22s", "distinct results");
   for (int a = 0; a < NB_SEQ_ALGOS; a++)
      printf(" %12d", countDistinct(res[a], NB_ORDERS));
   printf("\n");
}


/* -------------------------------------------------------------------- */
/* ompNaive          Naive parallel sum: OpenMP reduction(+)            */
/*                                                                      */
/* Input : x  array of n doubles, nt  number of threads                 */
/*         dynamic  0: schedule(static), 1: schedule(dynamic)           */
/* Output: the sum                                                      */
/* -------------------------------------------------------------------- */
static double ompNaive(const double x[], size_t n, int nt, int dynamic)
{
   double s = 0.0;

   if (dynamic)
   {
#pragma omp parallel for num_threads(nt) reduction(+:s) \
                         schedule(dynamic, DYN_CHUNK)
      for (size_t i = 0; i < n; i++)
         s += x[i];
   }
   else
   {
#pragma omp parallel for num_threads(nt) reduction(+:s) schedule(static)
      for (size_t i = 0; i < n; i++)
         s += x[i];
   }
   return s;
}


/* -------------------------------------------------------------------- */
/* ompKahan          Kahan sum per thread on a contiguous block, then   */
/*                   exact combination of the 2 nt values s_t and -c_t  */
/*                   in thread order (no dependence on the timing)      */
/*                                                                      */
/* Input : x  array of n doubles, nt  number of threads (<= MAX_THREADS)*/
/* Output: the sum                                                      */
/* -------------------------------------------------------------------- */
static double ompKahan(const double x[], size_t n, int nt)
{
   double part[2 * MAX_THREADS];

#pragma omp parallel num_threads(nt)
   {
      int    t  = omp_get_thread_num();
      int    T  = omp_get_num_threads();
      size_t lo = n * (size_t) t / (size_t) T;
      size_t hi = n * (size_t) (t + 1) / (size_t) T;
      double s  = 0.0, c = 0.0, y, u;

      for (size_t i = lo; i < hi; i++)
      {
         y = x[i] - c;
         u = s + y;
         c = (u - s) - y;
         s = u;
      }
      part[2 * t]     = s;
      part[2 * t + 1] = -c;
   }
   return sumExact(part, 2 * (size_t) nt);
}


/* -------------------------------------------------------------------- */
/* ompExact          Exact parallel sum: each thread builds the exact   */
/*                   partials of its block, the partials of all the     */
/*                   threads are then merged exactly                    */
/*                                                                      */
/* Input : x  array of n doubles, nt  number of threads (<= MAX_THREADS)*/
/* Output: the exact sum rounded once: independent of nt                */
/* -------------------------------------------------------------------- */
static double ompExact(const double x[], size_t n, int nt)
{
   static double partials[MAX_THREADS][FSUM_MAX_PARTIALS];
   static double merged[FSUM_MAX_PARTIALS];
   int           np[MAX_THREADS], nm = 0;

#pragma omp parallel num_threads(nt)
   {
      int    t  = omp_get_thread_num();
      int    T  = omp_get_num_threads();
      size_t lo = n * (size_t) t / (size_t) T;
      size_t hi = n * (size_t) (t + 1) / (size_t) T;

      np[t] = fsumAdd(partials[t], 0, x + lo, hi - lo);
   }
   for (int t = 0; t < nt; t++)
      nm = fsumAdd(merged, nm, partials[t], (size_t) np[t]);
   return fsumRound(merged, nm);
}


/* -------------------------------------------------------------------- */
/* parallelSums      Section 2: the 3 parallel sums for 1 .. 12 threads */
/*                                                                      */
/* Input : x  array of n doubles, exact  exact sum                      */
/* Output: prints s - exact and the time (best of NB_RUNS) per method   */
/* -------------------------------------------------------------------- */
static void parallelSums(const double x[], size_t n, double exact)
{
   double res[3][NB_THREADS], t, best[3];

   printf("%7s %14s %14s %14s %10s %10s %10s\n", "threads",
          "naive(+)", "Kahan/thread", "exact", "t naive", "t Kahan",
          "t exact");
   printf("%7s %14s %14s %14s %10s %10s %10s\n", "",
          "sum - exact", "sum - exact", "sum - exact", "ms", "ms", "ms");
   for (int k = 0; k < NB_THREADS; k++)
   {
      int nt = THREADS[k];

      best[0] = best[1] = best[2] = 1e300;
      for (int r = 0; r < NB_RUNS; r++)
      {
         t = omp_get_wtime();
         res[0][k] = ompNaive(x, n, nt, 0);
         t = omp_get_wtime() - t;
         if (t < best[0]) best[0] = t;

         t = omp_get_wtime();
         res[1][k] = ompKahan(x, n, nt);
         t = omp_get_wtime() - t;
         if (t < best[1]) best[1] = t;

         t = omp_get_wtime();
         res[2][k] = ompExact(x, n, nt);
         t = omp_get_wtime() - t;
         if (t < best[2]) best[2] = t;
      }
      printf("%7d %14.3e %14.3e %14.3e %10.2f %10.2f %10.2f\n", nt,
             res[0][k] - exact, res[1][k] - exact, res[2][k] - exact,
             1e3 * best[0], 1e3 * best[1], 1e3 * best[2]);
   }
   printf("%7s %14d %14d %14d\n", "distinct", countDistinct(res[0],
          NB_THREADS), countDistinct(res[1], NB_THREADS),
          countDistinct(res[2], NB_THREADS));
}


/* -------------------------------------------------------------------- */
/* repeatedRuns      Section 3: NB_REPEAT runs with 12 threads          */
/*                                                                      */
/* Input : x  array of n doubles, exact  exact sum                      */
/* Output: prints s - exact for every run                               */
/* -------------------------------------------------------------------- */
static void repeatedRuns(const double x[], size_t n, double exact)
{
   const int nt = 12;
   double    res[4][NB_REPEAT];

   for (int r = 0; r < NB_REPEAT; r++)
   {
      res[0][r] = ompNaive(x, n, nt, 0);
      res[1][r] = ompNaive(x, n, nt, 1);
      res[2][r] = ompKahan(x, n, nt);
      res[3][r] = ompExact(x, n, nt);
   }

   printf("%4s %16s %16s %16s %16s\n", "run", "naive static",
          "naive dynamic", "Kahan/thread", "exact");
   for (int r = 0; r < NB_REPEAT; r++)
      printf("%4d %16.3e %16.3e %16.3e %16.3e\n", r + 1,
             res[0][r] - exact, res[1][r] - exact, res[2][r] - exact,
             res[3][r] - exact);
   printf("%4s %16d %16d %16d %16d\n", "dist", countDistinct(res[0],
          NB_REPEAT), countDistinct(res[1], NB_REPEAT),
          countDistinct(res[2], NB_REPEAT), countDistinct(res[3],
          NB_REPEAT));
}


int main(int argc, char *argv[])
{
   size_t    n = DEFAULT_N;
   double  * x, * work, exact;
   mtState_t g;
   uint32_t  key[5] = { 0x123, 0x234, 0x345, 0x456, 2 };

   if (argc > 1)
      n = (size_t) strtod(argv[1], NULL);      /* 1e7 accepted          */
   if (n < MIN_N)
   {
      fprintf(stderr, "Usage: %s [N >= %d]   (default 1e7)\n", argv[0],
              MIN_N);
      return EXIT_FAILURE;
   }
   x    = malloc(n * sizeof *x);
   work = malloc(n * sizeof *work);
   if (x == NULL || work == NULL)
   {
      perror("malloc");
      return EXIT_FAILURE;
   }

   init_by_array(&g, key, 5);           /* same stream as U[-1,1] in    */
   for (size_t i = 0; i < n; i++)       /* sumBench                     */
      x[i] = 2.0 * genrand_res53(&g) - 1.0;
   exact = sumExact(x, n);

   printf("=== reproSum: N = %zu values U[-1,1], exact sum = %.17g ===\n",
          n, exact);
   printf("(ulp of the exact sum = %.3e; OpenMP max threads = %d)\n",
          nextafter(fabs(exact), INFINITY) - fabs(exact),
          omp_get_max_threads());

   printf("\n=== 1. Sequential sums, same values in 6 orders ===\n");
   sequentialOrders(x, work, n, exact);

   printf("\n=== 2. Parallel sums (OpenMP) versus number of threads ===\n");
   parallelSums(x, n, exact);

   printf("\n=== 3. Same binary, same data, 12 threads, %d runs ===\n",
          NB_REPEAT);
   repeatedRuns(x, n, exact);

   free(x);
   free(work);
   return EXIT_SUCCESS;
}
