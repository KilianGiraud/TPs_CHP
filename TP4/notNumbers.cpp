// -------------------------------------------------------------------- //
// notNumbers.cpp    Section "Numeric Operations don't always return    //
//                   Numbers": infinities and NaNs                      //
//                                                                      //
//  1. Codes of the subject: 2 * DBL_MAX, sqrt(-2), the test            //
//     if (x - x == 0) and the test x == x.                             //
//  2. Special results of usual operations and the IEEE 754 exception   //
//     flags they raise (<cfenv>): overflow, underflow, division by     //
//     zero, invalid, inexact.                                          //
//  3. Comparisons with a NaN, and functions that do NOT propagate a    //
//     NaN (pow(NaN, 0) = 1, fmax(NaN, 1) = 1, std::max depends on the  //
//     order of its arguments).                                         //
//                                                                      //
// The operands are volatile so that the compiler cannot compute the    //
// results at compile time: the operations and flags happen at run time.//
// Built with -ffast-math (notNumbers_fast), the NaN tests are removed  //
// by the compiler: compare the two outputs.                            //
//                                                                      //
// Compile : g++ -O2 -Wall -Wextra -o notNumbers notNumbers.cpp         //
// Run     : ./notNumbers                                               //
// -------------------------------------------------------------------- //
#include <iostream>
#include <cmath>
#include <cfloat>
#include <cfenv>
#include <cstdio>
#include <string>
#include <algorithm>

static volatile double vMax  = DBL_MAX;    // volatile: no constant folding
static volatile double vTwo  = 2.0;
static volatile double vZero = 0.0;
static volatile double vOne  = 1.0;
static volatile double vInf  = INFINITY;   // created before the tests


// -------------------------------------------------------------------- //
// subjectCode       Codes of the subject                               //
//                                                                      //
// Input : none                                                         //
// Output: prints 2 * DBL_MAX, sqrt(-2) and the result of the tests     //
//         x - x == 0 and x == x for x = 1, inf and NaN                 //
// -------------------------------------------------------------------- //
static void subjectCode()
{
   double x = vMax;

   std::cout << 2 * x << "\n";
   std::cout << std::sqrt(-vTwo) << "\n";

   const double VALUES[] = { 1.0, 2 * x, std::sqrt(-vTwo) };
   const char * NAMES[]  = { "1", "inf", "NaN" };

   for (int i = 0; i < 3; i++)
   {
      volatile double v = VALUES[i];

      std::printf("x = %-4s  x - x == 0 : %-5s   x == x : %-5s   "
                  "std::isnan(x) : %s\n", NAMES[i],
                  (v - v == 0) ? "true" : "false",
                  (v == v) ? "true" : "false",
                  std::isnan(v) ? "true" : "false");
   }
}


// -------------------------------------------------------------------- //
// flagsString       IEEE 754 exception flags currently raised          //
//                                                                      //
// Input : none                                                         //
// Output: names of the raised flags (fetestexcept)                     //
// -------------------------------------------------------------------- //
static std::string flagsString()
{
   const int    FLAGS[] = { FE_OVERFLOW, FE_UNDERFLOW, FE_DIVBYZERO,
                            FE_INVALID, FE_INEXACT };
   const char * NAMES[] = { "OVERFLOW", "UNDERFLOW", "DIVBYZERO",
                            "INVALID", "INEXACT" };
   std::string  s;

   for (int i = 0; i < 5; i++)
      if (std::fetestexcept(FLAGS[i]))
         s += (s.empty() ? "" : " ") + std::string(NAMES[i]);
   return s.empty() ? "-" : s;
}


// -------------------------------------------------------------------- //
// show              Prints one operation, its result and its flags     //
//                                                                      //
// Input : inLabel  the operation                                       //
//         inValue  its result (computed after feclearexcept)           //
// Output: one line of the table                                        //
// -------------------------------------------------------------------- //
static void show(const char *inLabel, double inValue)
{
   std::printf("%-26s = %-13g %s\n", inLabel, inValue,
               flagsString().c_str());
}


// -------------------------------------------------------------------- //
// specialResults    Usual operations giving inf, NaN, 0 or subnormals  //
//                                                                      //
// Input : none                                                         //
// Output: table operation | result | flags                             //
// -------------------------------------------------------------------- //
static void specialResults()
{
   double r;

#define OPERATION(label, expr) \
   std::feclearexcept(FE_ALL_EXCEPT); r = (expr); show(label, r)

   OPERATION("2 * DBL_MAX",              vTwo * vMax);
   OPERATION("-2 * DBL_MAX",             -vTwo * vMax);
   OPERATION("1 / 0",                    vOne / vZero);
   OPERATION("-1 / 0",                   -vOne / vZero);
   OPERATION("0 / 0",                    vZero / vZero);
   OPERATION("sqrt(-2)",                 std::sqrt(-vTwo));
   OPERATION("log(0)",                   std::log(vZero));
   OPERATION("log(-1)",                  std::log(-vOne));
   OPERATION("inf - inf",                vInf - vInf);
   OPERATION("0 * inf",                  vZero * vInf);
   OPERATION("inf + 1",                  vInf + vOne);
   OPERATION("exp(1000)",                std::exp(1000 * vOne));
   OPERATION("exp(-1000)",               std::exp(-1000 * vOne));
   OPERATION("DBL_MIN / 2  (subnormal)", DBL_MIN * vOne / vTwo);
   OPERATION("DBL_MIN / 3  (subnormal)", DBL_MIN * vOne / (vTwo + vOne));
   OPERATION("DBL_TRUE_MIN / 2",         DBL_TRUE_MIN * vOne / vTwo);
   OPERATION("0.1 + 0.2",                0.1 * vOne + 0.2);
   OPERATION("0.5 + 0.25",               0.5 * vOne + 0.25);
#undef OPERATION
}


// -------------------------------------------------------------------- //
// nanBehaviour      Comparisons with a NaN and NaN propagation         //
//                                                                      //
// Input : none                                                         //
// Output: prints the results                                           //
// -------------------------------------------------------------------- //
static void nanBehaviour()
{
   volatile double n = std::sqrt(-vTwo);           // a NaN
   double          one = vOne;

   std::printf("NaN == NaN : %d   NaN != NaN : %d   NaN < 1 : %d   "
               "NaN > 1 : %d   NaN >= 1 : %d   !(NaN < 1) : %d\n",
               n == n, n != n, n < one, n > one, n >= one, !(n < one));
   std::printf("NaN + 1 = %g   NaN * 0 = %g   (a NaN propagates...)\n",
               n + one, n * 0.0);
   std::printf("pow(NaN, 0) = %g   fmax(NaN, 1) = %g   hypot(inf, NaN) = "
               "%g   (... but not always)\n", std::pow(n, 0.0),
               std::fmax(n, one), std::hypot((double) vInf, n));
   std::printf("std::max(NaN, 1.0) = %g   std::max(1.0, NaN) = %g   "
               "(depends on the order)\n", std::max((double) n, one),
               std::max(one, (double) n));
}


int main()
{
#ifdef __FAST_MATH__
   std::cout << "*** compiled with -ffast-math ***\n";
#endif
   std::cout << "=== 1. Codes of the subject ===\n";
   subjectCode();

   std::cout << "\n=== 2. Special results and IEEE 754 exception flags ===\n";
   specialResults();

   std::cout << "\n=== 3. NaN behaviour ===\n";
   nanBehaviour();
   return 0;
}
