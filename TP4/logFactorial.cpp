// -------------------------------------------------------------------- //
// logFactorial.cpp  Section "Use Logarithms to Avoid Overflow and      //
//                   Underflow"                                         //
//                                                                      //
//  1. N! overflows a double for N > 170 (DBL_MAX ~ 1.8e308).           //
//  2. 200! / (190! 10!) computed directly (inf / inf = NaN), with the  //
//     logFactorial function of the subject, with lgamma (log of the    //
//     Gamma function, C99 / C++11) and exactly with 64-bit integers.   //
//  3. A binomial probability where the direct formula overflows AND    //
//     underflows: P(X = 1000) for X ~ B(2000, 1/2)                     //
//     = C(2000, 1000) / 2^2000 = (2e600) * (8.7e-603).                 //
//                                                                      //
// Compile : g++ -O2 -Wall -Wextra -o logFactorial logFactorial.cpp     //
// Run     : ./logFactorial                                             //
// -------------------------------------------------------------------- //
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <cstdio>

const int    N_BIG      = 200;           // 200! / (190! 10!)
const int    K_SMALL    = 10;
const int    N_BINOM    = 2000;          // X ~ B(2000, 1/2)
const int    K_BINOM    = 1000;
// exact P(X = 1000), computed with rational numbers (Python fractions)
const double EXACT_PROB = 0.01783901114585432;


// -------------------------------------------------------------------- //
// logFactorial      Function of the subject: log(n!) = sum of log(i)   //
//                   (simple but inefficient: n calls to log)           //
//                                                                      //
// Input : inN  integer >= 0                                            //
// Output: log(inN!)                                                    //
// -------------------------------------------------------------------- //
static double logFactorial(int inN)
{
   double sum = 0.0;

   for (int i = 2; i <= inN; ++i)
   {
      sum += std::log((double) i);
   }

   return sum;
}


// -------------------------------------------------------------------- //
// factorial         n! in double precision                             //
//                                                                      //
// Input : inN  integer >= 0                                            //
// Output: inN! (inf when inN > 170)                                    //
// -------------------------------------------------------------------- //
static double factorial(int inN)
{
   double f = 1.0;

   for (int i = 2; i <= inN; ++i)
      f *= i;
   return f;
}


// -------------------------------------------------------------------- //
// binomialExact     C(n, k) with 64-bit integers, exact as long as the //
//                   result times n fits in 64 bits                     //
//                                                                      //
// Input : inN, inK  0 <= inK <= inN                                    //
// Output: C(inN, inK); at step i, c = C(n - k + i, i) is an integer    //
// -------------------------------------------------------------------- //
static unsigned long long binomialExact(int inN, int inK)
{
   unsigned long long c = 1;

   for (int i = 1; i <= inK; ++i)
      c = c * (unsigned long long) (inN - inK + i) / (unsigned long long) i;
   return c;
}


// -------------------------------------------------------------------- //
// ratioOfFactorials 200! / (190! 10!) with four methods                //
//                                                                      //
// Input : none                                                         //
// Output: prints the results and their relative errors (the product    //
//         of the ratios avoids both the overflow and the cancellation  //
//         of the logs)                                                 //
// -------------------------------------------------------------------- //
static void ratioOfFactorials()
{
   unsigned long long exactInt = binomialExact(N_BIG, K_SMALL);
   double             exact    = (double) exactInt;
   double             direct   = factorial(N_BIG)
                               / (factorial(N_BIG - K_SMALL)
                                  * factorial(K_SMALL));
   double             lf200    = logFactorial(N_BIG);
   double             lf190    = logFactorial(N_BIG - K_SMALL);
   double             lf10     = logFactorial(K_SMALL);
   double             byLogs   = std::exp(lf200 - lf190 - lf10);
   double             byGamma  = std::exp(std::lgamma(N_BIG + 1.0)
                                          - std::lgamma(N_BIG - K_SMALL + 1.0)
                                          - std::lgamma(K_SMALL + 1.0));
   double             byRatios = 1.0;     // product of (n - k + i) / i

   for (int i = 1; i <= K_SMALL; ++i)
      byRatios = byRatios * (N_BIG - K_SMALL + i) / i;

   std::printf("170! = %.6e    171! = %g    DBL_MAX = %.6e\n",
               factorial(170), factorial(171), DBL_MAX);
   std::printf("200! = %g    190! = %g\n\n", factorial(N_BIG),
               factorial(N_BIG - K_SMALL));

   std::printf("exact (64-bit integers)     C(200,10) = %llu\n", exactInt);
   std::printf("direct 200!/(190! 10!)              = %g\n", direct);
   std::printf("logFactorial: log 200! = %.12f, log 190! = %.12f, "
               "log 10! = %.12f\n", lf200, lf190, lf10);
   std::printf("exp(logFact(200) - logFact(190) - logFact(10)) = %.17g"
               "   rel. error %.1e\n", byLogs,
               std::fabs(byLogs - exact) / exact);
   std::printf("exp(lgamma(201) - lgamma(191) - lgamma(11))    = %.17g"
               "   rel. error %.1e\n", byGamma,
               std::fabs(byGamma - exact) / exact);
   std::printf("product of (190 + i) / i, i = 1..10 (double)  = %.17g"
               "   rel. error %.1e\n", byRatios,
               std::fabs(byRatios - exact) / exact);
}


// -------------------------------------------------------------------- //
// binomialProbability  P(X = k) for X ~ B(n, 1/2)                      //
//                                                                      //
// Input : none                                                         //
// Output: prints the direct computation (NaN) and the log computation  //
// -------------------------------------------------------------------- //
static void binomialProbability()
{
   double c = 1.0;                            // C(n, k) in double

   for (int i = 1; i <= K_BINOM; ++i)
      c = c * (N_BINOM - K_BINOM + i) / i;
   double p      = std::pow(0.5, N_BINOM);
   double direct = c * p;
   double logP   = std::lgamma(N_BINOM + 1.0) - std::lgamma(K_BINOM + 1.0)
                 - std::lgamma(N_BINOM - K_BINOM + 1.0)
                 + N_BINOM * std::log(0.5);
   double byLogs = std::exp(logP);

   std::printf("C(2000,1000) = %g (overflow)   0.5^2000 = %g (underflow)\n",
               c, p);
   std::printf("direct  C(2000,1000) * 0.5^2000 = %g\n", direct);
   std::printf("logs    exp(log C(2000,1000) + 2000 log 0.5) = %.17g\n",
               byLogs);
   std::printf("exact   %.17g    rel. error %.1e\n", EXACT_PROB,
               std::fabs(byLogs - EXACT_PROB) / EXACT_PROB);
}


int main()
{
   std::cout << "=== 1-2. Ratio of factorials 200! / (190! 10!) ===\n";
   ratioOfFactorials();

   std::cout << "\n=== 3. Binomial probability P(X = 1000), "
                "X ~ B(2000, 1/2) ===\n";
   binomialProbability();
   return 0;
}
