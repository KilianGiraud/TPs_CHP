/* -------------------------------------------------------------------- */
/* compare.c         Section "How does it impact my computing?"         */
/*                                                                      */
/*  1. The two codes of the subject: the FLOATS a = 0.7 and b = 0.5     */
/*     are compared with the constants 0.7 and 0.5, which are DOUBLES.  */
/*  2. The "naive expected behaviour": 0.7f / 0.5f constants, or        */
/*     double variables.                                                */
/*  3. What is really stored for 0.1, 0.5, 0.6 and 0.7, in float and    */
/*     in double: IEEE 754 fields (sign | exponent | mantissa, as on    */
/*     binaryconvert.com, Figure 1 of the subject) and exact decimal    */
/*     value of the stored number.                                      */
/*  4. How many digits must be printed to see the stored value.         */
/*                                                                      */
/* Compile : gcc -O2 -Wall -Wextra -o compare compare.c                 */
/* Run     : ./compare                                                  */
/* -------------------------------------------------------------------- */
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define FLT_EXP_BITS     8      /* IEEE 754 binary32: 1 + 8 + 23 bits  */
#define FLT_MANT_BITS   23
#define DBL_EXP_BITS    11      /* IEEE 754 binary64: 1 + 11 + 52 bits */
#define DBL_MANT_BITS   52
#define EXACT_DECIMALS  70      /* enough for the exact decimal value  */
#define BUF_SIZE       128


/* -------------------------------------------------------------------- */
/* verdict           Message printed by the codes of the subject        */
/*                                                                      */
/* Input : condA  result of the test a < 0.7                            */
/*         condB  result of the test b < 0.5                            */
/*                                                                      */
/* Output: same nested logic as the subject (condB is only looked at    */
/*         when condA is true)                                          */
/* -------------------------------------------------------------------- */
static const char *verdict(int condA, int condB)
{
   if (condA)
   {
      if (condB) return "2 conditions are right";
      else       return "only 1 is right";
   }
   return "0 conditions are right";
}


/* -------------------------------------------------------------------- */
/* subjectCode1      First code of the subject, unchanged               */
/*                                                                      */
/* Input : none                                                         */
/* Output: prints how many conditions are "right"                       */
/* -------------------------------------------------------------------- */
static void subjectCode1(void)
{
   float a = 0.7;
   float b = 0.5;

   if (a < 0.7)
   {
      if (b < 0.5) printf("2 conditions are right");
      else         printf("only 1 is right");
   }
   else printf("0 conditions are right");
   printf("\n");
}


/* -------------------------------------------------------------------- */
/* subjectCode2      Second code of the subject (curly quotes of the    */
/*                   PDF replaced by "", newlines and labels added)     */
/*                                                                      */
/* Input : none                                                         */
/* Output: prints the float and double values of 0.7 and the result of  */
/*         the same tests as subjectCode1                               */
/* -------------------------------------------------------------------- */
static void subjectCode2(void)
{
   float  a = 0.7, b = 0.5;                   /* These are FLOATS      */
   double c = (double) a;
   double d = 0.7;

   printf("a  %%0.10f : %0.10f\n", a);
   printf("c  %%0.15f : %0.15f\n", c);

   printf("d  %%0.15f : %0.15f\n", d);        /* an approximation      */
   printf("c  %%0.18f : %0.18f\n", c);
   printf("d  %%0.18f : %0.18f\n", d);        /* "all the details"     */

   if (a < .7)                                /* test done in DOUBLE   */
   {
      if (b < .5) printf("2 are right");      /* in DOUBLE too         */
      else        printf("1 is right");
   }
   else printf("0 are right");
   printf("\n");
}


/* -------------------------------------------------------------------- */
/* expectedBehaviour The two corrections proposed by the subject        */
/*                                                                      */
/* Input : none                                                         */
/* Output: prints the verdict with float constants (.7f, .5f) and with  */
/*         double variables                                             */
/* -------------------------------------------------------------------- */
static void expectedBehaviour(void)
{
   float  af = 0.7f, bf = 0.5f;
   double ad = 0.7,  bd = 0.5;

   printf("float a, b compared with .7  and .5  (doubles) : %s\n",
          verdict(af < .7, bf < .5));
   printf("float a, b compared with .7f and .5f (floats)  : %s\n",
          verdict(af < .7f, bf < .5f));
   printf("double a, b compared with .7 and .5  (doubles) : %s\n",
          verdict(ad < .7, bd < .5));
}


/* -------------------------------------------------------------------- */
/* bitString         Binary string of the nbBits lowest bits of v       */
/*                                                                      */
/* Input : v       value                                                */
/*         nbBits  number of bits to write (most significant first)     */
/*         out     buffer of at least nbBits + 1 characters             */
/*                                                                      */
/* Output: out contains a string of '0' and '1'                         */
/* -------------------------------------------------------------------- */
static void bitString(uint64_t v, int nbBits, char out[])
{
   for (int i = 0; i < nbBits; i++)
      out[i] = ((v >> (nbBits - 1 - i)) & 1U) ? '1' : '0';
   out[nbBits] = '\0';
}


/* -------------------------------------------------------------------- */
/* exactDecimal      Exact decimal value of a double (every double is a */
/*                   finite binary fraction, so it has a finite decimal */
/*                   expansion; glibc's printf prints it exactly)       */
/*                                                                      */
/* Input : v    value (a float converted to double is exact)            */
/*         out  buffer of at least BUF_SIZE characters                  */
/*                                                                      */
/* Output: out contains the decimal expansion without trailing zeros    */
/* -------------------------------------------------------------------- */
static void exactDecimal(double v, char out[])
{
   size_t len;

   snprintf(out, BUF_SIZE, "%.*f", EXACT_DECIMALS, v);
   len = strlen(out);
   while (len > 2 && out[len - 1] == '0' && out[len - 2] != '.')
      out[--len] = '\0';
}


/* -------------------------------------------------------------------- */
/* printFields       IEEE 754 fields of the float and double nearest    */
/*                   to a decimal constant                              */
/*                                                                      */
/* Input : label  the decimal constant as written in the source         */
/*         f      the constant stored in a float                        */
/*         d      the constant stored in a double                       */
/*                                                                      */
/* Output: prints hexadecimal, sign, exponent and mantissa bits         */
/* -------------------------------------------------------------------- */
static void printFields(const char *label, float f, double d)
{
   uint32_t bf;
   uint64_t bd;
   char     exp[BUF_SIZE],
            mant[BUF_SIZE];

   memcpy(&bf, &f, sizeof bf);                 /* bits of the float    */
   memcpy(&bd, &d, sizeof bd);                 /* bits of the double   */

   bitString(bf >> FLT_MANT_BITS, FLT_EXP_BITS, exp);
   bitString(bf, FLT_MANT_BITS, mant);
   printf("%-4s float   0x%08X           %u  %-11s  %s\n", label,
          (unsigned) bf, (unsigned) (bf >> 31), exp, mant);

   bitString(bd >> DBL_MANT_BITS, DBL_EXP_BITS, exp);
   bitString(bd, DBL_MANT_BITS, mant);
   printf("%-4s double  0x%016llX   %u  %-11s  %s\n", label,
          (unsigned long long) bd, (unsigned) (bd >> 63), exp, mant);
}


/* -------------------------------------------------------------------- */
/* representation    Stored values of 0.1, 0.5, 0.6, 0.7                */
/*                                                                      */
/* Input : none                                                         */
/* Output: prints the IEEE fields, the exact decimal values and how the */
/*         float compares with the double (the test of the subject)     */
/* -------------------------------------------------------------------- */
static void representation(void)
{
   static const char  * LABELS[] = { "0.1", "0.5", "0.6", "0.7" };
   static const float   FLOATS[] = { 0.1f, 0.5f, 0.6f, 0.7f };
   static const double  DOUBLES[] = { 0.1, 0.5, 0.6, 0.7 };
   const int            nb = 4;
   char                 bufF[BUF_SIZE],
                        bufD[BUF_SIZE];

   printf("val  type    hexadecimal          s  exponent     mantissa\n");
   for (int i = 0; i < nb; i++)
      printFields(LABELS[i], FLOATS[i], DOUBLES[i]);

   printf("\nExact decimal value of the stored numbers:\n");
   for (int i = 0; i < nb; i++)
   {
      double fd = (double) FLOATS[i];          /* exact conversion     */

      exactDecimal(fd, bufF);
      exactDecimal(DOUBLES[i], bufD);
      printf("  %sf = %s\n", LABELS[i], bufF);
      printf("  %s  = %s\n", LABELS[i], bufD);
      printf("  => (double) %sf %s %s  so  a < %s is %s\n\n", LABELS[i],
             fd < DOUBLES[i] ? "< " : (fd > DOUBLES[i] ? "> " : "=="),
             LABELS[i], LABELS[i], fd < DOUBLES[i] ? "TRUE" : "false");
   }
}


/* -------------------------------------------------------------------- */
/* digits            Formats needed to see the stored value of 0.7      */
/*                                                                      */
/* Input : none                                                         */
/* Output: prints 0.7 (double) and 0.7f (float) with several formats    */
/* -------------------------------------------------------------------- */
static void digits(void)
{
   double d = 0.7;
   float  f = 0.7f;
   char   buf[BUF_SIZE];

   exactDecimal(d, buf);
   printf("double 0.7  %%.15f : %.15f\n", d);
   printf("            %%.18f : %.18f\n", d);
   printf("            %%.17g : %.17g   (17 digits: unique double)\n", d);
   printf("            %%a    : %a   (exact, hexadecimal)\n", d);
   printf("            exact : %s  (%d decimals)\n", buf,
          (int) strlen(buf) - 2);

   exactDecimal(f, buf);
   printf("float  0.7f %%.9g  : %.9g   (9 digits: unique float)\n", f);
   printf("            %%a    : %a\n", f);
   printf("            exact : %s  (%d decimals)\n", buf,
          (int) strlen(buf) - 2);
}


int main(void)
{
   printf("=== 1. First code of the subject ===\n");
   subjectCode1();

   printf("\n=== 2. Second code of the subject ===\n");
   subjectCode2();

   printf("\n=== 3. Expected behaviour (float constants or doubles) ===\n");
   expectedBehaviour();

   printf("\n=== 4. What is really stored in memory ===\n");
   representation();

   printf("=== 5. How many digits to print? ===\n");
   digits();
   return 0;
}
