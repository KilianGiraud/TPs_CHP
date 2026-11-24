// -------------------------------------------------------------------- //
// ranges.cpp        Section "Floating Point Numbers have Finite        //
//                   Ranges"                                            //
//                                                                      //
//  1. Codes of the subject: 2^24 + 1 in float, (2^53 + 1) - 2^53 in    //
//     double, 1 + DBL_EPSILON / 2 == 1 ("Sorry!").                     //
//  2. Details: f + 1 (float) versus f + 1.0 (double), what DBL_EPSILON //
//     really is (1 + eps/2 is a tie rounded to even; 1 - eps/2 != 1).  //
//  3. Spacing between two consecutive floats / doubles (ulp) for       //
//     numbers of different magnitudes: the absolute precision          //
//     decreases when the magnitude increases.                          //
//  4. Limits of float, double and long double (float.h).               //
//                                                                      //
// Compile : g++ -O2 -Wall -Wextra -o ranges ranges.cpp                 //
// Run     : ./ranges                                                   //
// -------------------------------------------------------------------- //
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <cstdio>


// -------------------------------------------------------------------- //
// subjectCode       The three codes of the subject                     //
//                                                                      //
// Input : none                                                         //
// Output: prints f and f + 1, (x + 1) - x, and "Sorry!"                //
// -------------------------------------------------------------------- //
static void subjectCode()
{
   float  f = 16777216;
   double x, y;

   std::cout << f << " " << f + 1 << "\n";

   x = 9007199254740992;                     // 2^53
   std::cout << ((x + 1) - x) << "\n";

   x = 1.0;
   y = x + 0.5 * DBL_EPSILON;

   if (x == y) std::cout << "Sorry!\n";
}


// -------------------------------------------------------------------- //
// details           What happens around 2^24, 2^53 and 1 + eps         //
//                                                                      //
// Input : none                                                         //
// Output: prints the values with enough digits                         //
// -------------------------------------------------------------------- //
static void details()
{
   float  f = 16777216;                      // 2^24
   double x = 9007199254740992;              // 2^53

   std::cout << std::setprecision(17)
             << "default precision (6) hides everything: " << std::setw(0)
             << std::setprecision(6) << f << " " << f + 1 << "\n"
             << std::setprecision(17)
             << "f           = " << f << "\n"
             << "f + 1       = " << f + 1 << "   (float + int -> float)\n"
             << "f + 1.0     = " << f + 1.0
             << "   (float + double -> double)\n"
             << "f + 2       = " << f + 2 << "   (2^24 + 2 is a float)\n"
             << "(x + 1) - x = " << (x + 1) - x
             << "   x + 1 = 2^53 + 1 is a tie, rounded to even 2^53\n"
             << "(x + 2) - x = " << (x + 2) - x << "\n"
             << "(x + 3) - x = " << (x + 3) - x
             << "   x + 3 is a tie, rounded to even 2^53 + 4\n\n";

   double eps = DBL_EPSILON;

   std::cout << "DBL_EPSILON           = " << eps << " = 2^-52\n"
             << "nextafter(1, 2) - 1   = " << std::nextafter(1.0, 2.0) - 1.0
             << "  (definition: gap after 1)\n" << std::boolalpha
             << "1 + 0.5  * eps == 1 : " << (1.0 + 0.5 * eps == 1.0)
             << "   tie between 1 and 1 + eps -> even -> 1\n"
             << "1 + 0.51 * eps == 1 : " << (1.0 + 0.51 * eps == 1.0)
             << "  above the tie -> rounded up to 1 + eps\n"
             << "1 - 0.5  * eps == 1 : " << (1.0 - 0.5 * eps == 1.0)
             << "  below 1 the gap is eps/2\n";
}


// -------------------------------------------------------------------- //
// spacing           Gap between a number and the next representable    //
//                   one, in float and in double                        //
//                                                                      //
// Input : none                                                         //
// Output: prints ulp(v) for several magnitudes                         //
// -------------------------------------------------------------------- //
static void spacing()
{
   const double VALUES[] = { 1.0, 1000.0, 1e6, 16777216.0, 1e9,
                             9007199254740992.0, 1e16, 1e30, 1e300 };
   const int    NB = sizeof VALUES / sizeof VALUES[0];

   std::printf("%22s %14s %14s\n", "value", "float ulp", "double ulp");
   for (int i = 0; i < NB; i++)
   {
      double v  = VALUES[i];
      double ud = std::nextafter(v, INFINITY) - v;

      if (v <= FLT_MAX)
      {
         float vf = (float) v;
         float uf = std::nextafter(vf, INFINITY) - vf;
         std::printf("%22.16g %14.7g %14.7g\n", v, uf, ud);
      }
      else
         std::printf("%22.16g %14s %14.7g\n", v, "(> FLT_MAX)", ud);
   }
}


// -------------------------------------------------------------------- //
// limits            Characteristics of the three types (float.h)       //
//                                                                      //
// Input : none                                                         //
// Output: prints mantissa bits, digits, epsilon, min, max              //
// -------------------------------------------------------------------- //
static void limits()
{
   std::printf("%-30s %14s %14s %14s\n", "", "float", "double",
               "long double");
   std::printf("%-30s %14d %14d %14d\n", "bits of the significand",
               FLT_MANT_DIG, DBL_MANT_DIG, LDBL_MANT_DIG);
   std::printf("%-30s %14d %14d %14d\n", "decimal digits always kept",
               FLT_DIG, DBL_DIG, LDBL_DIG);
   std::printf("%-30s %14d %14d %14d\n", "digits for a round trip",
               FLT_DECIMAL_DIG, DBL_DECIMAL_DIG, LDBL_DECIMAL_DIG);
   std::printf("%-30s %14.3e %14.3e %14.3Le\n", "EPSILON",
               FLT_EPSILON, DBL_EPSILON, LDBL_EPSILON);
   std::printf("%-30s %14.3e %14.3e %14.3Le\n", "smallest normal (MIN)",
               FLT_MIN, DBL_MIN, LDBL_MIN);
   std::printf("%-30s %14.3e %14.3e %14.3Le\n", "smallest subnormal",
               FLT_TRUE_MIN, DBL_TRUE_MIN, LDBL_TRUE_MIN);
   std::printf("%-30s %14.3e %14.3e %14.3Le\n", "largest (MAX)",
               FLT_MAX, DBL_MAX, LDBL_MAX);
   std::printf("%-30s %14.7g %14.7g %14.7Lg\n",
               "consecutive integers up to", std::ldexp(1.0, FLT_MANT_DIG),
               std::ldexp(1.0, DBL_MANT_DIG), ldexpl(1.0L, LDBL_MANT_DIG));
}


int main()
{
   std::cout << "=== 1. Codes of the subject ===\n";
   subjectCode();

   std::cout << "\n=== 2. Details ===\n";
   details();

   std::cout << "\n=== 3. Spacing between consecutive numbers (ulp) ===\n";
   spacing();

   std::cout << "\n=== 4. Limits (float.h) ===\n";
   limits();
   return 0;
}
