/* -------------------------------------------------------------------- */
/* exactSum.c        Exact sum of an array of doubles, rounded ONCE to  */
/*                   the nearest double (ties to even): algorithm of    */
/*                   Shewchuk (1997) as written in Python's math.fsum   */
/*                   (R. Hettinger, CPython Modules/mathmodule.c)       */
/*                                                                      */
/* The running sum is kept EXACTLY as a list of "partials": doubles     */
/* that do not overlap (the lowest set bit of one is above the highest  */
/* bit of the previous one), by increasing magnitude, whose exact       */
/* mathematical sum is the exact sum of the values added so far. A new  */
/* value is merged with each partial by Fast2Sum (hi = x + y and        */
/* lo = y - (hi - x) are exact when |x| >= |y|); zero errors are        */
/* dropped. The final rounding sums the partials from the largest and   */
/* corrects the half-way cases.                                         */
/*                                                                      */
/* The partials of several sub-arrays can be merged exactly (fsumAdd on */
/* the partials): reproSum uses this for an exact parallel sum.         */
/*                                                                      */
/* REFERENCE of sumBench and reproSum: this file must never be compiled */
/* with -ffast-math (checked below). It is validated at the start of    */
/* sumBench with the test cases of CPython (test_math.py, test_fsum).   */
/* Assumption: finite values whose partial sums do not overflow.        */
/* -------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "summation.h"

#ifdef __FAST_MATH__
#error "exactSum.c must not be built with -ffast-math (wrong reference)"
#endif


/* -------------------------------------------------------------------- */
/* fsumAdd           Adds n values to an exact list of partials         */
/*                                                                      */
/* Input : p   partials (capacity FSUM_MAX_PARTIALS), np of them used   */
/*         x   array of n doubles to add                                */
/* Output: p updated; returns the new number of partials                */
/* -------------------------------------------------------------------- */
int fsumAdd(double p[], int np, const double x[], size_t n)
{
   double xi, y, t, hi, lo;
   int    i;

   for (size_t k = 0; k < n; k++)
   {
      xi = x[k];
      i  = 0;
      for (int j = 0; j < np; j++)             /* for y in partials     */
      {
         y = p[j];
         if (fabs(xi) < fabs(y))               /* Fast2Sum: |xi| >= |y| */
         {
            t = xi; xi = y; y = t;
         }
         hi = xi + y;
         lo = y - (hi - xi);                   /* exact error of xi + y */
         if (lo != 0.0)
            p[i++] = lo;
         xi = hi;
      }
      np = i;
      if (xi != 0.0)
      {
         if (np >= FSUM_MAX_PARTIALS || !isfinite(xi))
         {
            fprintf(stderr, "fsumAdd: overflow or too many partials\n");
            exit(EXIT_FAILURE);
         }
         p[np++] = xi;
      }
   }
   return np;
}


/* -------------------------------------------------------------------- */
/* fsumRound         Nearest double to the exact sum of the partials    */
/*                                                                      */
/* Input : p   np partials built by fsumAdd                             */
/* Output: the exact sum rounded to nearest, ties to even               */
/* -------------------------------------------------------------------- */
double fsumRound(const double p[], int np)
{
   double hi = 0.0, lo = 0.0, x, y, yr;
   int    n  = np;

   if (n > 0)
   {
      hi = p[--n];
      while (n > 0)       /* sum from the top while the sum stays exact */
      {
         x  = hi;
         y  = p[--n];
         hi = x + y;
         yr = hi - x;
         lo = y - yr;
         if (lo != 0.0)
            break;
      }
      /* half-way case: if the remaining partials push in the same      */
      /* direction as lo, the true sum is above the tie: round up       */
      /* (needed for e.g. 1e16 + 1 + 1e-16 = 10000000000000002)         */
      if (n > 0 && ((lo < 0.0 && p[n - 1] < 0.0)
                    || (lo > 0.0 && p[n - 1] > 0.0)))
      {
         y  = lo * 2.0;
         x  = hi + y;
         yr = x - hi;
         if (y == yr)
            hi = x;
      }
   }
   return hi;
}


/* -------------------------------------------------------------------- */
/* sumExact          Exact sum of an array, correctly rounded           */
/*                                                                      */
/* Input : x  array of n doubles                                        */
/* Output: the double nearest to the exact sum                          */
/* -------------------------------------------------------------------- */
double sumExact(const double x[], size_t n)
{
   double p[FSUM_MAX_PARTIALS];
   int    np = fsumAdd(p, 0, x, n);

   return fsumRound(p, np);
}
