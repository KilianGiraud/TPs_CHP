/* -------------------------------------------------------------------- */
/* softmax.c         "To go a little further - For AI students":        */
/*                   why a naive softmax fails, and the stable version  */
/*                                                                      */
/*   softmax(z)_i = exp(z_i) / sum_j exp(z_j)                           */
/*                                                                      */
/*  - naive: exp(z) overflows to inf for z > 88.7 in float (z > 709.8   */
/*    in double) => inf / inf = NaN; it underflows to 0 for very        */
/*    negative logits => 0 / 0 = NaN.                                   */
/*  - stable: softmax is unchanged when the same constant is subtracted */
/*    from every logit (exp(z_i - m) = exp(z_i) / exp(m)). With         */
/*    m = max z_j, every exponent is <= 0 (no overflow) and the         */
/*    denominator is >= 1 (no division by 0).                           */
/*  - log-softmax (log-sum-exp trick), used by the cross-entropy loss:  */
/*    log softmax(z)_i = z_i - m - log(sum_j exp(z_j - m)), finite even */
/*    when softmax(z)_i underflows to 0 (where log(0) = -inf).          */
/*                                                                      */
/* Computed in float (32 bits), as in deep learning frameworks.         */
/*                                                                      */
/* Compile : gcc -O2 -Wall -Wextra -o softmax softmax.c -lm             */
/* Run     : ./softmax                                                  */
/* -------------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#define K        4                   /* number of classes               */
#define NB_SETS  5                   /* number of logit vectors tested  */


/* -------------------------------------------------------------------- */
/* softmaxNaive      Direct formula, float                              */
/*                                                                      */
/* Input : z  K logits                                                  */
/* Output: p  K probabilities (NaN when exp overflows or underflows)    */
/* -------------------------------------------------------------------- */
static void softmaxNaive(const float z[], float p[])
{
   float sum = 0.0f;

   for (int i = 0; i < K; i++)
      sum += expf(z[i]);
   for (int i = 0; i < K; i++)
      p[i] = expf(z[i]) / sum;
}


/* -------------------------------------------------------------------- */
/* softmaxNaiveDouble Direct formula in double (overflow at 709.8)      */
/*                                                                      */
/* Input : z  K logits                                                  */
/* Output: p  K probabilities, rounded to float                         */
/* -------------------------------------------------------------------- */
static void softmaxNaiveDouble(const float z[], float p[])
{
   double sum = 0.0;

   for (int i = 0; i < K; i++)
      sum += exp((double) z[i]);
   for (int i = 0; i < K; i++)
      p[i] = (float) (exp((double) z[i]) / sum);
}


/* -------------------------------------------------------------------- */
/* maxOf             Largest logit                                      */
/*                                                                      */
/* Input : z  K logits                                                  */
/* Output: max z_i                                                      */
/* -------------------------------------------------------------------- */
static float maxOf(const float z[])
{
   float m = z[0];

   for (int i = 1; i < K; i++)
      if (z[i] > m) m = z[i];
   return m;
}


/* -------------------------------------------------------------------- */
/* softmaxStable     exp(z_i - max) / sum exp(z_j - max), float         */
/*                                                                      */
/* Input : z  K logits                                                  */
/* Output: p  K probabilities (sum = 1 up to rounding)                  */
/* -------------------------------------------------------------------- */
static void softmaxStable(const float z[], float p[])
{
   float m   = maxOf(z);
   float sum = 0.0f;

   for (int i = 0; i < K; i++)
   {
      p[i] = expf(z[i] - m);          /* in ]0, 1], at least one is 1   */
      sum += p[i];
   }
   for (int i = 0; i < K; i++)
      p[i] /= sum;
}


/* -------------------------------------------------------------------- */
/* logSoftmax        z_i - max - log(sum exp(z_j - max)), float         */
/*                                                                      */
/* Input : z  K logits                                                  */
/* Output: lp K log-probabilities (always finite)                       */
/* -------------------------------------------------------------------- */
static void logSoftmax(const float z[], float lp[])
{
   float m   = maxOf(z);
   float sum = 0.0f;

   for (int i = 0; i < K; i++)
      sum += expf(z[i] - m);
   for (int i = 0; i < K; i++)
      lp[i] = z[i] - m - logf(sum);
}


/* -------------------------------------------------------------------- */
/* printVector       One line: label and K values                       */
/*                                                                      */
/* Input : label, v  K values                                           */
/* Output: prints the line                                              */
/* -------------------------------------------------------------------- */
static void printVector(const char *label, const float v[])
{
   printf("  %-15s", label);
   for (int i = 0; i < K; i++)
      printf(" %12.6g", v[i]);
   printf("\n");
}


int main(void)
{
   static const float LOGITS[NB_SETS][K] = {
      {     1.0f,     2.0f,     3.0f,     4.0f },
      {  1000.0f,  1001.0f,  1002.0f,  1003.0f },
      { -1000.0f, -1001.0f, -1002.0f, -1003.0f },
      {   100.0f,    50.0f,     0.0f,   -50.0f },
      {     0.0f,     0.0f,     0.0f,  -200.0f } };
   float p[K], lp[K], sum;

   printf("=== Softmax of %d logits (float) ===\n", K);
   for (int s = 0; s < NB_SETS; s++)
   {
      printVector("logits", LOGITS[s]);
      softmaxNaive(LOGITS[s], p);
      printVector("naive float", p);
      softmaxNaiveDouble(LOGITS[s], p);
      printVector("naive double", p);
      softmaxStable(LOGITS[s], p);
      printVector("stable", p);
      sum = 0.0f;
      for (int i = 0; i < K; i++) sum += p[i];
      printf("  %-15s %12.9f\n\n", "sum (stable)", sum);
   }

   printf("=== Cross-entropy loss -log p[3] for the last logits ===\n");
   softmaxNaive(LOGITS[NB_SETS - 1], p);
   logSoftmax(LOGITS[NB_SETS - 1], lp);
   printf("  naive        -log(softmax) = %g   (p[3] = %g underflows)\n",
          -logf(p[K - 1]), p[K - 1]);
   printf("  log-softmax  -logSoftmax   = %g   (exact: 200 + log 3 = "
          "%.4f)\n", -lp[K - 1], 200.0 + log(3.0));
   return 0;
}
