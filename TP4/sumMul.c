/* -------------------------------------------------------------------- */
/* sumMul.c          Section "Worry about Addition and Subtraction      */
/*                   (more than Multiplication and Division)",          */
/*                   first code: 10 additions of 0.1f versus 0.1f * 10  */
/*                                                                      */
/*  1. Code of the subject.                                             */
/*  2. f * 10 (float product) versus f * 10.0 (double product).         */
/*  3. The same experiment with n = 10 .. 10^8 additions: the error of  */
/*     the float sum grows with n, and the sum finally stops growing    */
/*     (absorption: 0.1 is smaller than half the spacing between two    */
/*     consecutive floats around the sum). A product is rounded only    */
/*     once. The Kahan compensated sum (see sumBench) fixes the sum.    */
/*                                                                      */
/* Compile : gcc -O2 -Wall -Wextra -o sumMul sumMul.c                   */
/* Run     : ./sumMul                                                   */
/* -------------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#define MAX_POWER  8           /* n = 10^1 .. 10^MAX_POWER additions   */


/* -------------------------------------------------------------------- */
/* subjectCode       Code of the subject                                */
/*                                                                      */
/* Input : none                                                         */
/* Output: prints the float sum and the float products                  */
/* -------------------------------------------------------------------- */
static void subjectCode(void)
{
   float f   = 0.1f;
   float sum = 0;

   for (int i = 0; i < 10; ++i) sum += f;

   float product = f * 10;

   printf("sum = %1.15f, mul = %1.15f, mul2 = %1.15f\n", sum, product,
          f * 10);
}


/* -------------------------------------------------------------------- */
/* sumFloat          n additions of x in a float accumulator            */
/*                                                                      */
/* Input : x  term, n  number of additions                              */
/* Output: the float sum                                                */
/* -------------------------------------------------------------------- */
static float sumFloat(float x, long n)
{
   float s = 0.0f;

   for (long i = 0; i < n; i++)
      s += x;
   return s;
}


/* -------------------------------------------------------------------- */
/* sumKahanFloat     n additions of x with Kahan's compensated sum      */
/*                   (float accumulator, see summation.c)               */
/*                                                                      */
/* Input : x  term, n  number of additions                              */
/* Output: the compensated float sum                                    */
/* -------------------------------------------------------------------- */
static float sumKahanFloat(float x, long n)
{
   float s = 0.0f, c = 0.0f, y, t;

   for (long i = 0; i < n; i++)
   {
      y = x - c;                      /* term corrected by the error   */
      t = s + y;                      /* low bits of y are lost here   */
      c = (t - s) - y;                /* -(lost bits), recovered       */
      s = t;
   }
   return s;
}


/* -------------------------------------------------------------------- */
/* sumDouble         n additions of x in a double accumulator           */
/*                                                                      */
/* Input : x  term, n  number of additions                              */
/* Output: the double sum                                               */
/* -------------------------------------------------------------------- */
static double sumDouble(double x, long n)
{
   double s = 0.0;

   for (long i = 0; i < n; i++)
      s += x;
   return s;
}


/* -------------------------------------------------------------------- */
/* growingN          Sum of n terms 0.1f versus n * 0.1f                */
/*                                                                      */
/* Input : none                                                         */
/* Output: prints, for n = 10 .. 10^MAX_POWER, the float sum, the float */
/*         product, the Kahan float sum and their relative errors       */
/*         against the exact value n * 0.1f; then the same sum of 0.1   */
/*         in double                                                    */
/* -------------------------------------------------------------------- */
static void growingN(void)
{
   const float f = 0.1f;
   long        n = 1;

   printf("%10s %14s %8s %12s %8s %12s %8s\n", "n", "float sum",
          "rel.err", "n * 0.1f", "rel.err", "Kahan float", "rel.err");
   for (int k = 1; k <= MAX_POWER; k++)
   {
      n *= 10;
      /* n * 0.1f in double is EXACT: 24 + 27 bits < 53 bits          */
      double exact = (double) n * (double) f;
      float  s     = sumFloat(f, n);
      float  p     = (float) n * f;
      float  k2    = sumKahanFloat(f, n);

      printf("%10ld %14.4f %8.1e %12.1f %8.1e %12.1f %8.1e\n", n,
             s, fabs(s - exact) / exact, p, fabs(p - exact) / exact,
             k2, fabs(k2 - exact) / exact);
   }

   printf("\nSame sum of 0.1 in DOUBLE (53-bit mantissa):\n");
   printf("%10s %22s %8s\n", "n", "double sum", "rel.err");
   n = 1;
   for (int k = 1; k <= MAX_POWER; k++)
   {
      n *= 10;
      long double exact = (long double) n * 0.1;   /* error ~ 1e-19   */
      double      s     = sumDouble(0.1, n);

      printf("%10ld %22.10f %8.1e\n", n, s,
             (double) (fabsl(s - exact) / exact));
   }
}


int main(void)
{
   const float f = 0.1f;

   printf("=== 1. Code of the subject ===\n");
   subjectCode();

   printf("\n=== 2. Float product versus double product ===\n");
   printf("f * 10    (float)  = %1.15f\n", f * 10);
   printf("f * 10.0  (double) = %1.15f   (exact value of 0.1f * 10)\n",
          f * 10.0);

   printf("\n=== 3. n additions of 0.1f ===\n");
   growingN();
   return 0;
}
