/* ------------------------------------------------------------------ */
/* simPi.c - Monte Carlo estimation of pi                              */
/*                                                                    */
/* Draw points (xr, yr) uniformly in the unit square [0,1]^2 and count */
/* those falling inside the quarter disc x^2 + y^2 <= 1.               */
/* Area ratio = (pi/4) / 1  =>  pi ~= 4 * inside / n                   */
/*                                                                    */
/* NOTE: rand() is used here only as a first approach, as requested.  */
/* It is statistically weak; use Mersenne Twister / MRG32k3a for real  */
/* scientific work.                                                    */
/*                                                                    */
/* Compile : gcc -O2 -Wall -o simPi simPi.c -lm                        */
/* Run     : ./simPi            (default test sequence)                */
/*           ./simPi 1000000    (single run with n points)             */
/* ------------------------------------------------------------------ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Uniform random number in [0, 1] */
static inline double uniform(void)
{
    return (double) rand() / (double) RAND_MAX;
}

/* ------------------------------------------------------------------ */
/* simPi : estimate pi using nbPoints random points                    */
/* ------------------------------------------------------------------ */
double simPi(long long nbPoints)
{
    long long inside = 0;
    long long i;
    double xr, yr;

    for (i = 0; i < nbPoints; i++)
    {
        xr = uniform();
        yr = uniform();
        if (xr * xr + yr * yr <= 1.0)
            inside++;
    }
    return 4.0 * (double) inside / (double) nbPoints;
}

/* Run one estimation and print result, absolute error and time */
static void runTest(long long n)
{
    clock_t start = clock();
    double  piEst = simPi(n);
    double  secs  = (double)(clock() - start) / CLOCKS_PER_SEC;
    double  err   = fabs(piEst - M_PI);
    /* theoretical standard deviation of the estimator */
    double  sigma = sqrt(M_PI * (4.0 - M_PI) / (double) n);

    printf("n = %14lld | pi ~ %.8f | |error| = %.2e | sigma_th = %.2e | %7.2f s\n",
           n, piEst, err, sigma, secs);
}

int main(int argc, char *argv[])
{
    srand(12345);   /* fixed seed => reproducible experiment */

    if (argc > 1)
    {
        runTest(atoll(argv[1]));
        return 0;
    }

    printf("--- Required tests ---\n");
    runTest(1000LL);
    runTest(1000000LL);
    runTest(1000000000LL);

    printf("\n--- Precision study (powers of 10) ---\n");
    for (long long n = 10; n <= 1000000000LL; n *= 10)
        runTest(n);

    return 0;
}