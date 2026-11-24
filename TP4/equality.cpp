// -------------------------------------------------------------------- //
// equality.cpp      Sections "Don't test for equality" and             //
//                   "Relative errors"                                  //
//                                                                      //
//  1. Code of the subject: sqrt(10) squared is not exactly 10.         //
//  2. The difference measured in ulps (unit in the last place).        //
//  3. FLT_EPSILON and DBL_EPSILON (float.h).                           //
//  4. Absolute tolerance versus relative tolerance on several pairs:   //
//     the absolute test is wrong for large AND for tiny numbers.       //
//  5. A comparison function with a relative tolerance and an absolute  //
//     floor (needed near 0).                                           //
//                                                                      //
// Compile : g++ -O2 -Wall -Wextra -o equality equality.cpp             //
// Run     : ./equality                                                 //
// -------------------------------------------------------------------- //
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <cstdio>

const double ABS_TOL = 1e-12;      // "1e-12 for double precision"
const double REL_TOL = 1e-12;


// -------------------------------------------------------------------- //
// subjectCode       Code of the subject (using namespace std replaced  //
//                   by explicit std::)                                 //
//                                                                      //
// Input : none                                                         //
// Output: prints x - y when sqrt(x)^2 differs from x                   //
// -------------------------------------------------------------------- //
static void subjectCode()
{
   double x = 10;
   double y = std::sqrt(x);

   y *= y;

   if (x == y)
      std::cout << "Square root is exact\n";
   else
      std::cout << x - y << "\n";
}


// -------------------------------------------------------------------- //
// ulpsBetween       Distance between two doubles in ulps of inX        //
//                                                                      //
// Input : inX, inY  two doubles of the same binade                     //
// Output: (inY - inX) / ulp(inX)                                       //
// -------------------------------------------------------------------- //
static double ulpsBetween(double inX, double inY)
{
   return (inY - inX) / (std::nextafter(inX, INFINITY) - inX);
}


// -------------------------------------------------------------------- //
// relativeError     Relative error of the subject                      //
//                                                                      //
// Input : inResult, inExpected  computed and expected values           //
// Output: |inResult - inExpected| / |inExpected|  (inf if expected 0)  //
// -------------------------------------------------------------------- //
static double relativeError(double inResult, double inExpected)
{
   return std::fabs((inResult - inExpected) / inExpected);
}


// -------------------------------------------------------------------- //
// almostEqual       Recommended comparison: relative tolerance with an //
//                   absolute floor for values close to 0 (same rule as //
//                   Python's math.isclose)                             //
//                                                                      //
// Input : inA, inB   values to compare                                 //
//         inRelTol   relative tolerance (e.g. 1e-12)                   //
//         inAbsTol   absolute floor (0 = purely relative)              //
// Output: true when |a - b| <= max(relTol * max(|a|, |b|), absTol)     //
// -------------------------------------------------------------------- //
static bool almostEqual(double inA, double inB, double inRelTol,
                        double inAbsTol)
{
   double diff = std::fabs(inA - inB);
   double big  = std::fmax(std::fabs(inA), std::fabs(inB));

   return diff <= std::fmax(inRelTol * big, inAbsTol);
}


// -------------------------------------------------------------------- //
// printPair         One line of the tolerance table                    //
//                                                                      //
// Input : inLabel  description of the pair                             //
//         inX      expected value                                      //
//         inY      computed value                                      //
// Output: prints x - y, the absolute test, the relative error and the  //
//         relative test                                                //
// -------------------------------------------------------------------- //
static void printPair(const char *inLabel, double inX, double inY)
{
   bool absEqual = std::fabs(inX - inY) < ABS_TOL;
   bool relEqual = almostEqual(inX, inY, REL_TOL, 0.0);

   std::printf("%-26s %10.2e  %-9s %10.2e  %s\n", inLabel, inX - inY,
               absEqual ? "equal" : "DIFFERENT", relativeError(inY, inX),
               relEqual ? "equal" : "DIFFERENT");
}


int main()
{
   std::cout << "=== 1. Code of the subject ===\n";
   subjectCode();

   double x = 10;
   double y = std::sqrt(x) * std::sqrt(x);

   std::cout << std::setprecision(17)
             << "\n=== 2. Details ===\n"
             << "sqrt(10)        = " << std::sqrt(x) << "\n"
             << "sqrt(10)^2      = " << y << "\n"
             << "x - y           = " << x - y << "\n"
             << "y - x in ulps   = " << ulpsBetween(x, y)
             << "   (ulp(10) = 2^-49 = "
             << std::nextafter(x, INFINITY) - x << ")\n";

   std::cout << "\n=== 3. Ready to use constants (float.h) ===\n"
             << "FLT_EPSILON = " << std::setprecision(9) << FLT_EPSILON
             << "  = 2^-23\n"
             << "DBL_EPSILON = " << std::setprecision(17) << DBL_EPSILON
             << "  = 2^-52\n";

   std::printf("\n=== 4. Absolute tolerance (%.0e) vs relative "
               "tolerance (%.0e) ===\n", ABS_TOL, REL_TOL);
   std::printf("%-26s %10s  %-9s %10s  %s\n", "pair (expected, computed)",
               "x - y", "abs test", "rel. err", "rel test");
   printPair("10, sqrt(10)^2", 10.0, std::sqrt(10.0) * std::sqrt(10.0));
   printPair("2e10, sqrt(2e10)^2", 2e10, std::sqrt(2e10) * std::sqrt(2e10));
   printPair("0.3, 0.1 + 0.2", 0.3, 0.1 + 0.2);
   printPair("1e-15, 2e-15", 1e-15, 2e-15);
   printPair("0, sin(pi)", 0.0, std::sin(M_PI));

   std::cout << "\n=== 5. Relative errors ===\n" << std::setprecision(6)
             << "result 99.5, expected 100 : relative error = "
             << relativeError(99.5, 100.0) << "\n"
             << "near 0 a relative test is useless: almostEqual(0, sin(pi),"
             << " 1e-12, 0)     = "
             << std::boolalpha << almostEqual(0.0, std::sin(M_PI), 1e-12, 0)
             << "\n"
             << "an absolute floor is needed:       almostEqual(0, sin(pi),"
             << " 1e-12, 1e-15) = "
             << almostEqual(0.0, std::sin(M_PI), 1e-12, 1e-15) << "\n";
   return 0;
}
