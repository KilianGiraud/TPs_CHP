/* -------------------------------------------------------------------- */
/* summation.h       Summation algorithms for arrays of floating point  */
/*                   numbers ("To go a little further", Q2 and Q3)      */
/*                                                                      */
/*   naive        s += x[i]                       error <= (n-1) eps S  */
/*   pairwise     recursive halves (as in NumPy)  error ~ log2(n) eps S */
/*   Kahan        compensated sum (1965)          error <= 2 eps S + .. */
/*   Neumaier     Kahan-Babuska-Neumaier (1974): also right when a term */
/*                is larger than the running sum                        */
/*   TwoSum       Sum2 of Ogita, Rump, Oishi (2005): error-free TwoSum  */
/*                of Knuth, as if computed with twice the precision     */
/*   Klein        second-order Kahan-Babuska (2006): ~ 3x precision     */
/*   long double  80-bit x87 accumulator (64-bit mantissa)              */
/*   exact        exact sum rounded once (Shewchuk 1997 / Python fsum)  */
/*   ...8         same with 8 independent accumulators (in the order    */
/*                i mod 8): no chain of dependent additions, faster     */
/*   (S = sum of |x_i|, eps = 2^-53 for double, 2^-24 for float)        */
/*                                                                      */
/* Declarations only (ISIMA style guide): the code is in summation.c    */
/* and exactSum.c.                                                      */
/*                                                                      */
/* WARNING: compensated sums rely on the exact IEEE 754 semantics of    */
/* every operation. They must NOT be compiled with -ffast-math or       */
/* -Ofast: the compiler would simplify the compensation                 */
/* c = ((s + y) - s) - y to 0 (Goldberg 1991, section 3.2.3).           */
/* summation.c is built with -ffast-math only in sumBench_fast, on      */
/* purpose, to show this trap.                                          */
/* -------------------------------------------------------------------- */
#ifndef SUMMATION_H
#define SUMMATION_H

#include <stddef.h>

#define PAIRWISE_BLOCK      128   /* below: naive loop (as in NumPy)    */
#define NB_LANES              8   /* accumulators of the ...8 versions  */
#define FSUM_MAX_PARTIALS  2100   /* partials do not overlap: at most   */
                                  /* one per bit position of a double   */
                                  /* (2^-1074 .. 2^1023)                */

/* ------------------------------ double ------------------------------ */
double sumNaive     (const double x[], size_t n);
double sumNaive8    (const double x[], size_t n);
double sumPairwise  (const double x[], size_t n);
double sumKahan     (const double x[], size_t n);
double sumKahan8    (const double x[], size_t n);
double sumNeumaier  (const double x[], size_t n);
double sumTwoSum    (const double x[], size_t n);
double sumKlein     (const double x[], size_t n);
double sumLongDouble(const double x[], size_t n);

/* ---------------------- exact sum (exactSum.c) ---------------------- */
int    fsumAdd      (double p[], int np, const double x[], size_t n);
double fsumRound    (const double p[], int np);
double sumExact     (const double x[], size_t n);

/* ------------------------- float data ------------------------------- */
float  sumNaiveF    (const float x[], size_t n);
float  sumPairwiseF (const float x[], size_t n);
float  sumKahanF    (const float x[], size_t n);
double sumNaiveFD   (const float x[], size_t n);  /* double accumulator */

/* ------------------------- build information ------------------------ */
const char * summationBuild(void);

#endif /* SUMMATION_H */
