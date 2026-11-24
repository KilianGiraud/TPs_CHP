/* -------------------------------------------------------------------- */
/* summation.c       Summation algorithms (see summation.h)             */
/*                                                                      */
/* Every function sums x[0] .. x[n-1] in this order (except the         */
/* pairwise and the 8-accumulator versions, whose order is given in     */
/* their comment). Compiled normally (-O2, no -ffast-math) in sumBench  */
/* and reproSum, with -O3 -ffast-math in sumBench_fast to show how the  */
/* compiler destroys the compensations.                                 */
/* -------------------------------------------------------------------- */
#include <math.h>
#include "summation.h"

#ifdef __FAST_MATH__
#warning "summation.c built with -ffast-math: compensations will be lost"
#endif


/* -------------------------------------------------------------------- */
/* summationBuild    How this file was compiled                         */
/*                                                                      */
/* Input : none                                                         */
/* Output: description of the floating point mode of the build          */
/* -------------------------------------------------------------------- */
const char *summationBuild(void)
{
#ifdef __FAST_MATH__
   return "-ffast-math (the compiler may reorder and simplify)";
#else
   return "IEEE 754 strict (no -ffast-math)";
#endif
}


/* ==================================================================== */
/*  Double precision                                                    */
/* ==================================================================== */

/* -------------------------------------------------------------------- */
/* sumNaive          Recursive summation s = s + x[i]                   */
/*                                                                      */
/* Input : x  array of n doubles                                        */
/* Output: the sum; error <= (n - 1) eps sum|x_i| (+ O(eps^2))          */
/* -------------------------------------------------------------------- */
double sumNaive(const double x[], size_t n)
{
   double s = 0.0;

   for (size_t i = 0; i < n; i++)
      s += x[i];
   return s;
}


/* -------------------------------------------------------------------- */
/* sumNaive8         Naive sum with NB_LANES independent accumulators   */
/*                   (accumulator j sums the x[i] with i mod 8 = j),    */
/*                   then the accumulators are added                    */
/*                                                                      */
/* Input : x  array of n doubles                                        */
/* Output: the sum (another order than sumNaive, so another rounding)   */
/* -------------------------------------------------------------------- */
double sumNaive8(const double x[], size_t n)
{
   double s[NB_LANES] = { 0.0 };
   double total = 0.0;
   size_t nv = n - n % NB_LANES, i;

   for (i = 0; i < nv; i += NB_LANES)
      for (int j = 0; j < NB_LANES; j++)
         s[j] += x[i + j];
   for (; i < n; i++)                         /* last n mod 8 values   */
      s[0] += x[i];

   for (int j = 0; j < NB_LANES; j++)
      total += s[j];
   return total;
}


/* -------------------------------------------------------------------- */
/* sumPairwise       Pairwise (cascade) summation: sum of the first     */
/*                   half + sum of the second half, recursively; naive  */
/*                   loop below PAIRWISE_BLOCK values (as in NumPy)     */
/*                                                                      */
/* Input : x  array of n doubles                                        */
/* Output: the sum; error ~ (PAIRWISE_BLOCK + log2 n) eps sum|x_i|      */
/* -------------------------------------------------------------------- */
double sumPairwise(const double x[], size_t n)
{
   size_t half;

   if (n <= PAIRWISE_BLOCK)
      return sumNaive(x, n);
   half = n / 2;
   return sumPairwise(x, half) + sumPairwise(x + half, n - half);
}


/* -------------------------------------------------------------------- */
/* sumKahan          Kahan compensated summation (1965), Theorem 8 of   */
/*                   Goldberg (1991)                                    */
/*                                                                      */
/*  t = s + y loses the low-order bits of y; (t - s) is the part of y   */
/*  that was really added, so c = (t - s) - y = -(lost bits), exactly   */
/*  when |s| >= |y|. The lost bits are given back to the next term.     */
/*                                                                      */
/* Input : x  array of n doubles                                        */
/* Output: the sum; error <= (2 eps + O(n eps^2)) sum|x_i|              */
/* -------------------------------------------------------------------- */
double sumKahan(const double x[], size_t n)
{
   double s = 0.0,                   /* running sum                     */
          c = 0.0,                   /* compensation: -(lost bits)      */
          y, t;

   for (size_t i = 0; i < n; i++)
   {
      y = x[i] - c;                  /* term corrected by the lost bits */
      t = s + y;                     /* low-order bits of y lost here   */
      c = (t - s) - y;               /* -(lost bits), computed exactly  */
      s = t;
   }
   return s;
}


/* -------------------------------------------------------------------- */
/* sumKahan8         Kahan summation with NB_LANES independent (s, c)   */
/*                   pairs (lane j sums the x[i] with i mod 8 = j); the */
/*                   16 partial results are then added with Neumaier    */
/*                                                                      */
/* Input : x  array of n doubles                                        */
/* Output: the sum, Kahan accuracy, without the chain of 4 dependent    */
/*         operations per term that slows down sumKahan                 */
/* -------------------------------------------------------------------- */
double sumKahan8(const double x[], size_t n)
{
   double s[NB_LANES] = { 0.0 },
          c[NB_LANES] = { 0.0 },
          rest[2 * NB_LANES],
          y, t;
   size_t nv = n - n % NB_LANES, i;

   for (i = 0; i < nv; i += NB_LANES)
      for (int j = 0; j < NB_LANES; j++)
      {
         y    = x[i + j] - c[j];
         t    = s[j] + y;
         c[j] = (t - s[j]) - y;
         s[j] = t;
      }
   for (; i < n; i++)                         /* last n mod 8 values   */
   {
      y    = x[i] - c[0];
      t    = s[0] + y;
      c[0] = (t - s[0]) - y;
      s[0] = t;
   }

   for (int j = 0; j < NB_LANES; j++)          /* s_j - c_j, 16 values */
   {
      rest[2 * j]     = s[j];
      rest[2 * j + 1] = -c[j];
   }
   return sumNeumaier(rest, 2 * NB_LANES);
}


/* -------------------------------------------------------------------- */
/* sumNeumaier       Kahan-Babuska-Neumaier summation (1974): the error */
/*                   of s + x[i] is computed from the LARGER operand,   */
/*                   so it stays exact when |x[i]| > |s| (case where    */
/*                   Kahan fails, e.g. 1 + 1e100 + 1 - 1e100)           */
/*                                                                      */
/* Input : x  array of n doubles                                        */
/* Output: s + c. Each error added to c is EXACT (Fast2Sum applied to   */
/*         the larger operand), as in sumTwoSum: both functions add the */
/*         same errors in the same order and return the same result     */
/* -------------------------------------------------------------------- */
double sumNeumaier(const double x[], size_t n)
{
   double s = 0.0,                   /* running sum                     */
          c = 0.0,                   /* sum of the exact errors         */
          t;

   for (size_t i = 0; i < n; i++)
   {
      t = s + x[i];
      if (fabs(s) >= fabs(x[i]))
         c += (s - t) + x[i];        /* low-order bits of x[i] lost     */
      else
         c += (x[i] - t) + s;        /* low-order bits of s lost        */
      s = t;
   }
   return s + c;
}


/* -------------------------------------------------------------------- */
/* sumTwoSum         Sum2 of Ogita, Rump and Oishi (2005): each         */
/*                   addition s + x[i] is followed by Knuth's TwoSum,   */
/*                   which gives its EXACT error e without any test     */
/*                   (6 operations); the errors are summed in c         */
/*                                                                      */
/* Input : x  array of n doubles                                        */
/* Output: s + c, as accurate as the naive sum computed with twice the  */
/*         precision then rounded: error <= eps|sum| + g^2 sum|x_i|,    */
/*         g = n eps / (1 - n eps)                                      */
/* -------------------------------------------------------------------- */
double sumTwoSum(const double x[], size_t n)
{
   double s = 0.0,
          c = 0.0,
          t, z, e;

   for (size_t i = 0; i < n; i++)
   {
      t  = s + x[i];                 /* rounded sum                     */
      z  = t - s;                    /* part of x[i] really added       */
      e  = (s - (t - z)) + (x[i] - z);   /* exact error: s + x = t + e  */
      s  = t;
      c += e;
   }
   return s + c;
}


/* -------------------------------------------------------------------- */
/* sumKlein          Second-order iterative Kahan-Babuska algorithm     */
/*                   (A. Klein, Computing 76, 2006): the errors of the  */
/*                   Neumaier compensation are themselves compensated   */
/*                                                                      */
/* Input : x  array of n doubles                                        */
/* Output: s + (cs + ccs)                                               */
/* -------------------------------------------------------------------- */
double sumKlein(const double x[], size_t n)
{
   double s   = 0.0,                 /* running sum                     */
          cs  = 0.0,                 /* first-order compensation        */
          ccs = 0.0,                 /* second-order compensation       */
          t, c, cc;

   for (size_t i = 0; i < n; i++)
   {
      t = s + x[i];
      if (fabs(s) >= fabs(x[i])) c = (s - t) + x[i];
      else                       c = (x[i] - t) + s;
      s = t;

      t = cs + c;
      if (fabs(cs) >= fabs(c))   cc = (cs - t) + c;
      else                       cc = (c - t) + cs;
      cs   = t;
      ccs += cc;
   }
   return s + (cs + ccs);
}


/* -------------------------------------------------------------------- */
/* sumLongDouble     Naive sum in a wider accumulator: long double =    */
/*                   x87 80-bit format on x86-64 (64-bit mantissa,      */
/*                   eps = 2^-64)                                       */
/*                                                                      */
/* Input : x  array of n doubles                                        */
/* Output: the sum rounded to double                                    */
/* -------------------------------------------------------------------- */
double sumLongDouble(const double x[], size_t n)
{
   long double s = 0.0L;

   for (size_t i = 0; i < n; i++)
      s += x[i];
   return (double) s;
}


/* ==================================================================== */
/*  Single precision data                                               */
/* ==================================================================== */

/* -------------------------------------------------------------------- */
/* sumNaiveF         Naive sum of floats in a float                     */
/*                                                                      */
/* Input : x  array of n floats                                         */
/* Output: the float sum (eps = 2^-24 ~ 6e-8)                           */
/* -------------------------------------------------------------------- */
float sumNaiveF(const float x[], size_t n)
{
   float s = 0.0f;

   for (size_t i = 0; i < n; i++)
      s += x[i];
   return s;
}


/* -------------------------------------------------------------------- */
/* sumPairwiseF      Pairwise summation of floats in float              */
/*                                                                      */
/* Input : x  array of n floats                                         */
/* Output: the float sum                                                */
/* -------------------------------------------------------------------- */
float sumPairwiseF(const float x[], size_t n)
{
   size_t half;

   if (n <= PAIRWISE_BLOCK)
      return sumNaiveF(x, n);
   half = n / 2;
   return sumPairwiseF(x, half) + sumPairwiseF(x + half, n - half);
}


/* -------------------------------------------------------------------- */
/* sumKahanF         Kahan summation of floats in float                 */
/*                                                                      */
/* Input : x  array of n floats                                         */
/* Output: the float sum, error <= (2 eps + O(n eps^2)) sum|x_i|        */
/* -------------------------------------------------------------------- */
float sumKahanF(const float x[], size_t n)
{
   float s = 0.0f, c = 0.0f, y, t;

   for (size_t i = 0; i < n; i++)
   {
      y = x[i] - c;
      t = s + y;
      c = (t - s) - y;
      s = t;
   }
   return s;
}


/* -------------------------------------------------------------------- */
/* sumNaiveFD        Naive sum of floats in a DOUBLE accumulator (the   */
/*                   "simple way" of Goldberg, section 3.2.3)           */
/*                                                                      */
/* Input : x  array of n floats                                         */
/* Output: the double sum (each float is converted exactly)             */
/* -------------------------------------------------------------------- */
double sumNaiveFD(const float x[], size_t n)
{
   double s = 0.0;

   for (size_t i = 0; i < n; i++)
      s += x[i];
   return s;
}
