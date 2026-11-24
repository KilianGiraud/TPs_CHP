// -------------------------------------------------------------------- //
// derivative.cpp    Section "Worry about Addition and Subtraction",    //
//                   second code: derivative of sin at x = 1 by         //
//                   (f(x+h) - f(x)) / h for h = 10^-1 .. 10^-19        //
//                                                                      //
//  1. Code of the subject (setprecision(15)).                          //
//  2. Same values with the digits after the first wrong digit replaced //
//     by '.' (presentation of the subject) and the error, compared     //
//     with two "tricks that change the problem":                       //
//       - central difference (f(x+h) - f(x-h)) / 2h : truncation error //
//         in h^2 instead of h, best h ~ eps^(1/3) instead of eps^(1/2) //
//       - complex step Im(f(x + ih)) / h : NO subtraction at all, so   //
//         no cancellation, full precision for any small h              //
//  3. Why "gcc on Linux" gave -0.017, -0.17, -1.7, -17 instead of 0    //
//     in the subject: x87 extended precision (80 bits) mixed with      //
//     double precision (64 bits), reproduced here with long double.    //
//                                                                      //
// Compile : g++ -O2 -Wall -Wextra -o derivative derivative.cpp         //
// Run     : ./derivative                                               //
// -------------------------------------------------------------------- //
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cfloat>                   // FLT_EVAL_METHOD
#include <complex>
#include <cstdio>
#include <cstring>
#include <string>

const int    MAX_I  = 19;           // h = 10^-1 .. 10^-MAX_I
const int    DIGITS = 15;           // decimals compared
const double X0     = 1.0;          // derivative of sin at x = 1


// -------------------------------------------------------------------- //
// subjectCode       Code of the subject (std:: instead of using        //
//                   namespace std)                                     //
//                                                                      //
// Input : none                                                         //
// Output: prints (sin(1+h) - sin(1)) / h for h = 10^-1 .. 10^-19       //
// -------------------------------------------------------------------- //
static void subjectCode()
{
   std::cout << std::setprecision(15);   // C++ way to set the precision

   for (int i = 1; i < 20; ++i)
   {
      double h = std::pow(10.0, -i);
      std::cout << (std::sin(1.0 + h) - std::sin(1.0)) / h << "\n";
   }

   std::cout << "True result: " << std::cos(1.0) << "\n";
}


// -------------------------------------------------------------------- //
// markDigits        Presentation of the subject: the digits after the  //
//                   first incorrect digit are replaced by '.'          //
//                                                                      //
// Input : inValue  computed value                                      //
//         inExact  exact value                                         //
// Output: string of DIGITS decimals with the wrong digits hidden       //
// -------------------------------------------------------------------- //
static std::string markDigits(double inValue, double inExact)
{
   char   v[64], e[64];
   size_t i;

   std::snprintf(v, sizeof v, "%.*f", DIGITS, inValue);
   std::snprintf(e, sizeof e, "%.*f", DIGITS, inExact);

   for (i = 0; v[i] != '\0' && v[i] == e[i]; i++)   // common prefix
      ;
   if (v[i] == '\0' || std::strlen(v) != std::strlen(e))
      return v;                     // all digits right, or other format
   for (size_t j = i + 1; v[j] != '\0'; j++)   // keep 1st wrong digit
      v[j] = '.';
   return v;
}


// -------------------------------------------------------------------- //
// compareMethods    Forward, central and complex-step differences      //
//                                                                      //
// Input : none                                                         //
// Output: prints, for each h, the three approximations of cos(1) and   //
//         their absolute errors                                        //
// -------------------------------------------------------------------- //
static void compareMethods()
{
   const double exact = std::cos(X0);

   std::printf("%6s  %-18s %8s  %-18s %8s  %-18s %8s\n", "h",
               "forward", "error", "central", "error", "complex step",
               "error");
   for (int i = 1; i <= MAX_I; i++)
   {
      double h   = std::pow(10.0, -i);
      double fwd = (std::sin(X0 + h) - std::sin(X0)) / h;
      double ctr = (std::sin(X0 + h) - std::sin(X0 - h)) / (2.0 * h);
      double cpx = std::imag(std::sin(std::complex<double>(X0, h))) / h;

      std::printf("1e-%-3d  %-18s %8.1e  %-18s %8.1e  %-18s %8.1e\n", i,
                  markDigits(fwd, exact).c_str(), std::fabs(fwd - exact),
                  markDigits(ctr, exact).c_str(), std::fabs(ctr - exact),
                  markDigits(cpx, exact).c_str(), std::fabs(cpx - exact));
   }
   std::printf("exact   %.*f\n", DIGITS, exact);
}


// -------------------------------------------------------------------- //
// x87Mix            Reproduces the values given by "gcc on Linux" in   //
//                   the subject: sin(1.0 + h) is rounded to a double   //
//                   (64 bits), but sin(1.0) is kept in an 80-bit x87   //
//                   register (FLT_EVAL_METHOD = 2 on 32-bit x86)       //
//                                                                      //
// Input : none                                                         //
// Output: prints (double(sin(1 + h)) - sinl(1)) / h for small h        //
// -------------------------------------------------------------------- //
static void x87Mix()
{
   long double sin1Ext = sinl(1.0L);             // 64-bit mantissa
   long double sin1Dbl = std::sin(1.0);          // 53-bit mantissa

   std::printf("FLT_EVAL_METHOD here = %d (0: every operation rounded to"
               " its type, SSE2)\n", (int) FLT_EVAL_METHOD);
   std::printf("double(sin 1) - extended(sin 1) = %.6Le\n",
               sin1Dbl - sin1Ext);
   for (int i = 14; i <= MAX_I; i++)
   {
      double      h = std::pow(10.0, -i);
      long double a = std::sin(1.0 + h);         // rounded to double
      std::printf("h = 1e-%d : %12.6Lf\n", i, (a - sin1Ext) / h);
   }
}


int main()
{
   std::cout << "=== 1. Code of the subject ===\n";
   subjectCode();

   std::cout << "\n=== 2. Correct digits and errors, three methods ===\n";
   compareMethods();

   std::cout << "\n=== 3. Mixing 80-bit and 64-bit (old gcc, x87) ===\n";
   x87Mix();
   return 0;
}
