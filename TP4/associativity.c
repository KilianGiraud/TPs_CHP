/* -------------------------------------------------------------------- */
/* associativity.c   Section "Out of order execution": floating point   */
/*                   addition and multiplication are NOT associative,   */
/*                   so changing the order of the operations changes    */
/*                   the result.                                        */
/*                                                                      */
/*  1. Examples of the subject: (1e-3 + 1) - 1 versus 1e-3 + (1 - 1)    */
/*     and (2^60 - 2^60) + 1 versus (2^60 + 1) - 2^60.                  */
/*  2. Other additions and multiplications, overflow, and the fused     */
/*     multiply-add (FMA: a * b + c with ONE rounding). When FMA        */
/*     instructions are enabled (-mfma, -march=native), gcc fuses       */
/*     p = a * b; p - c into one FMA by default (-ffp-contract=fast in  */
/*     GNU C mode): the SAME source gives another result. Build         */
/*     associativity_fma (Makefile) to see it; -ffp-contract=off or     */
/*     -std=c17 forbid this contraction.                                */
/*  3. A long sum in two orders: the harmonic series 1 + 1/2 + ... +    */
/*     1/n in float, forward and backward. Forward, the sum stops       */
/*     growing when 1/i is smaller than half the gap between two floats */
/*     around the sum.                                                  */
/*                                                                      */
/* The operands are volatile: the operations are done at run time.      */
/*                                                                      */
/* Compile : gcc -O2 -Wall -Wextra -o associativity associativity.c -lm */
/* Run     : ./associativity                                            */
/* -------------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#define N_HARMONIC  10000000L  /* number of terms of the harmonic sum  */
/* exact value of H(10^7), Python math.fsum (correctly rounded)         */
#define H_EXACT     16.695311365859855


/* -------------------------------------------------------------------- */
/* subjectExamples   The two examples of the subject                    */
/*                                                                      */
/* Input : none                                                         */
/* Output: prints both orders with %.18f                                */
/* -------------------------------------------------------------------- */
static void subjectExamples(void)
{
   volatile double one = 1.0, milli = pow(10.0, -3.0);
   volatile double big = pow(2.0, 60.0);

   printf("(pow(10,-3) + 1) - 1 = %.18f\n", (milli + one) - one);
   printf("pow(10,-3) + (1 - 1) = %.18f\n", milli + (one - one));
   printf("(2^60 - 2^60) + 1    = %.1f\n", (big - big) + one);
   printf("(2^60 + 1) - 2^60    = %.1f\n", (big + one) - big);
}


/* -------------------------------------------------------------------- */
/* otherExamples     More additions, multiplications and an FMA         */
/*                                                                      */
/* Input : none                                                         */
/* Output: prints both orders with 17 significant digits                */
/* -------------------------------------------------------------------- */
static void otherExamples(void)
{
   volatile double a = 0.1, b = 0.2, c = 0.3;
   volatile double huge = 1e308, ten = 10.0, tenth = 0.1, one = 1.0;
   double          p;

   printf("0.1 + 0.2            = %.17g  (== 0.3 : %s)\n", a + b,
          (a + b == c) ? "true" : "false");
   printf("(0.1 + 0.2) + 0.3    = %.17g\n", (a + b) + c);
   printf("0.1 + (0.2 + 0.3)    = %.17g\n", a + (b + c));
   printf("(0.1 * 0.2) * 0.3    = %.17g\n", (a * b) * c);
   printf("0.1 * (0.2 * 0.3)    = %.17g\n", a * (b * c));
   printf("(1e308 * 10) * 0.1   = %g   (overflow of the intermediate)\n",
          (huge * ten) * tenth);
   printf("1e308 * (10 * 0.1)   = %g\n", huge * (ten * tenth));

   p = tenth * ten;                    /* rounded product: exactly 1   */
   printf("p = 0.1 * 10; p - 1  = %.17g   (2 roundings, unless "
          "fused)\n", p - one);
   printf("fma(0.1, 10, -1)     = %.17g   (1 rounding: exact error of "
          "0.1 * 10)\n", fma(tenth, ten, -one));
}


/* -------------------------------------------------------------------- */
/* harmonic          Harmonic sum in float in two orders                */
/*                                                                      */
/* Input : none                                                         */
/* Output: prints the forward and backward float sums, the double sum   */
/*         and the iteration where the forward float sum stops growing  */
/* -------------------------------------------------------------------- */
static void harmonic(void)
{
   float  forward  = 0.0f,
          backward = 0.0f;
   double fwdDbl   = 0.0;
   long   stall    = 0;               /* first i where s + 1/i == s     */

   for (long i = 1; i <= N_HARMONIC; i++)
   {
      float next = forward + 1.0f / (float) i;

      if (next == forward && stall == 0)
         stall = i;
      forward = next;
      fwdDbl += 1.0 / (double) i;
   }
   for (long i = N_HARMONIC; i >= 1; i--)
      backward += 1.0f / (float) i;

   printf("n = %ld terms, exact H(n) = %.15f\n", N_HARMONIC, H_EXACT);
   printf("float  forward  (1 -> n) = %.7f   rel. error %.1e   "
          "(stops growing at i = %ld)\n", forward,
          fabs(forward - H_EXACT) / H_EXACT, stall);
   printf("float  backward (n -> 1) = %.7f   rel. error %.1e\n", backward,
          fabs(backward - H_EXACT) / H_EXACT);
   printf("double forward  (1 -> n) = %.15f   rel. error %.1e\n", fwdDbl,
          fabs(fwdDbl - H_EXACT) / H_EXACT);
}


int main(void)
{
#ifdef __FMA__
   printf("*** FMA instructions enabled: gcc may fuse a * b + c ***\n");
#endif
   printf("=== 1. Examples of the subject ===\n");
   subjectExamples();

   printf("\n=== 2. Other non-associative operations ===\n");
   otherExamples();

   printf("\n=== 3. Harmonic sum: order of the additions ===\n");
   harmonic();
   return 0;
}
