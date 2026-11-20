/* ------------------------------------------------------------------ */
/* sphereCI.c - Q2: Monte Carlo volume of a sphere (3D) and estimation */
/*              of PI with a confidence interval around the mean.     */
/*              Replicates computed in parallel with OpenMP.          */
/*                                                                    */
/* One replicate (experiment): nbPoints points drawn uniformly in the */
/* cube [0,R]^3 (one octant, by symmetry). The fraction f of points   */
/* inside the sphere x^2 + y^2 + z^2 <= R^2 converges to pi/6, so     */
/*      V  = 4/3 pi R^3  ~=  8 R^3 f         and   pi ~= 6 f          */
/* With dim = 2 it is the quarter disc of TP2 / openMPcomputePI.c:    */
/*      V  = pi R^2      ~=  4 R^2 f         and   pi ~= 4 f          */
/*                                                                    */
/* Each replicate owns its generator (no shared state between threads):*/
/*   mt  : MT19937, init_by_array with Makoto's original key          */
/*         {0x123, 0x234, 0x345, 0x456} extended with the replicate # */
/*   mrg : MRG32k3a, replicate r uses stream r (jump of r x 2^127)    */
/* The stream only depends on the replicate number, not on the thread */
/* => same results whatever the number of threads.                    */
/*                                                                    */
/* Confidence interval: method of the appendix of the subject         */
/*   X(n) mean, S^2(n) unbiased variance, S(n)/sqrt(n) standard error */
/*   R = t(n-1, 1-alpha/2) * sqrt(S^2(n) / n)  ->  [X(n) - R, X(n) + R]*/
/*                                                                    */
/* Compile : gcc -O2 -Wall -fopenmp -o sphereCI sphereCI.c -lm        */
/* Run     : ./sphereCI <mt|mrg> <nbReplicates> [nbPoints] [dim]      */
/*           ./sphereCI mrg 10              (10^9 points, 3D)         */
/*           ./sphereCI mt 40 1e6 2         (10^6 points, 2D)         */
/*           OMP_NUM_THREADS=4 ./sphereCI mt 10                       */
/* ------------------------------------------------------------------ */

/* Make math.h expose M_PI (POSIX / MSVC) so it comes from <math.h>    */
#define _DEFAULT_SOURCE
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "mt19937ar.h"
#include "mrg32k3a.h"

#define RADIUS          1.0     /* radius of the sphere              */
#define MAX_REPLICATES  50      /* Student table used up to nu = 40  */

/* ================================================================== */
/*  Student law: t(nu, 1 - alpha/2), nu = n - 1 degrees of freedom    */
/*  T975: Table 1 of the appendix (alpha = 0.05, 95 % interval)       */
/*  T995: same table for alpha = 0.01 (99 % interval)                 */
/*  Index = nu (the first column of Table 1 is nu, not n:             */
/*  12.706 is the value for 1 degree of freedom, i.e. n = 2).         */
/* ================================================================== */
static const double T975[31] = { 0.0,  /* nu = 0: no interval */
    12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.308, 2.262, 2.228,
     2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
     2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042 };
static const double T995[31] = { 0.0,
    63.657, 9.925, 5.841, 4.604, 4.032, 3.707, 3.499, 3.355, 3.250, 3.169,
     3.106, 3.055, 3.012, 2.977, 2.947, 2.921, 2.898, 2.878, 2.861, 2.845,
     2.831, 2.819, 2.807, 2.797, 2.787, 2.779, 2.771, 2.763, 2.756, 2.750 };
#define T975_NU40  2.021
#define T995_NU40  2.704

/* t(nu) read in the table; t40 = value of the table for nu = 40      */
static double studentT(int nu, const double tab[31], double t40)
{
    if (nu <= 30)
        return tab[nu];
    if (nu >= 40)           /* nu = 40 .. 49: value of nu = 40,        */
        return t40;         /* slightly larger => conservative interval */

    /* 30 < nu < 40: linear interpolation in 1/nu between the two rows */
    double w = (1.0 / nu - 1.0 / 40.0) / (1.0 / 30.0 - 1.0 / 40.0);
    return t40 + w * (tab[30] - t40);
}

/* ================================================================== */
/*  Statistics of a sample of n replicates (appendix of the subject)  */
/* ================================================================== */
typedef struct
{
    double mean;            /* X(n)                                  */
    double s2;              /* S^2(n), unbiased estimate of sigma^2  */
    double s;               /* S(n)                                  */
    double se;              /* standard error of the mean S/sqrt(n)  */
} stats_t;

/* Two passes: the deviations are taken from the mean, which avoids  */
/* the cancellation of the formula E[X^2] - E[X]^2 (values ~ 3.14,    */
/* variance ~ 1e-9).                                                  */
static stats_t computeStats(const double x[], int n)
{
    stats_t st;
    double  sum = 0.0, sq = 0.0;

    for (int i = 0; i < n; i++)
        sum += x[i];
    st.mean = sum / n;

    for (int i = 0; i < n; i++)
        sq += (x[i] - st.mean) * (x[i] - st.mean);
    st.s2 = sq / (n - 1);
    st.s  = sqrt(st.s2);
    st.se = st.s / sqrt((double) n);
    return st;
}

/* ================================================================== */
/*  One replicate: fraction of points inside the sphere (or disc)     */
/*  Same code for both generators, only the drawing function changes */
/* ================================================================== */
static double fractionMT(mt_state *g, long long nbPoints, int dim)
{
    long long inside = 0;
    double    x, y, z, r2;

    for (long long i = 0; i < nbPoints; i++)
    {
        x  = RADIUS * genrand_real1(g);
        y  = RADIUS * genrand_real1(g);
        r2 = x * x + y * y;
        if (dim == 3)
        {
            z   = RADIUS * genrand_real1(g);
            r2 += z * z;
        }
        if (r2 <= RADIUS * RADIUS)
            inside++;
    }
    return (double) inside / (double) nbPoints;
}

static double fractionMRG(mrg_state *g, long long nbPoints, int dim)
{
    long long inside = 0;
    double    x, y, z, r2;

    for (long long i = 0; i < nbPoints; i++)
    {
        x  = RADIUS * MRG32k3a(g);
        y  = RADIUS * MRG32k3a(g);
        r2 = x * x + y * y;
        if (dim == 3)
        {
            z   = RADIUS * MRG32k3a(g);
            r2 += z * z;
        }
        if (r2 <= RADIUS * RADIUS)
            inside++;
    }
    return (double) inside / (double) nbPoints;
}

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s <mt|mrg> <nbReplicates> [nbPoints] [dim]\n"
            "  nbReplicates : 2 .. %d\n"
            "  nbPoints     : per replicate, default 1e9\n"
            "  dim          : 3 (sphere, default) or 2 (disc)\n",
            prog, MAX_REPLICATES);
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    int       useMT, nbReplicates, dim = 3;
    long long nbPoints = 1000000000LL;

    /* ---------------- arguments ---------------- */
    if (argc < 3)
        usage(argv[0]);
    if      (strcmp(argv[1], "mt")  == 0) useMT = 1;
    else if (strcmp(argv[1], "mrg") == 0) useMT = 0;
    else    usage(argv[0]);

    nbReplicates = atoi(argv[2]);
    if (argc > 3) nbPoints = (long long) strtod(argv[3], NULL); /* 1e9 ok */
    if (argc > 4) dim      = atoi(argv[4]);

    if (nbReplicates < 2 || nbReplicates > MAX_REPLICATES
        || nbPoints <= 0 || (dim != 2 && dim != 3))
        usage(argv[0]);

    /* pi = k f with k = 4 (2D) or 6 (3D); exact volume of the ball    */
    double k        = (dim == 3) ? 6.0 : 4.0;
    double cube     = pow(2.0 * RADIUS, dim);    /* volume of [-R,R]^d */
    double exactVol = (dim == 3) ? 4.0 / 3.0 * M_PI * pow(RADIUS, 3)
                                 : M_PI * RADIUS * RADIUS;

    double *piRep  = malloc(nbReplicates * sizeof(double));
    double *volRep = malloc(nbReplicates * sizeof(double));
    if (piRep == NULL || volRep == NULL) { perror("malloc"); return EXIT_FAILURE; }

    /* ---------------- parallel replicates ---------------- */
    double start = omp_get_wtime();

    /* dynamic: a thread takes the next replicate as soon as it is free */
    #pragma omp parallel for schedule(dynamic, 1)
    for (int r = 0; r < nbReplicates; r++)
    {
        double f;

        if (useMT)
        {
            uint32_t key[5] = { 0x123, 0x234, 0x345, 0x456, (uint32_t) r };
            mt_state g;
            init_by_array(&g, key, 5);
            f = fractionMT(&g, nbPoints, dim);
        }
        else
        {
            mrg_state g;
            mrg_init_stream(&g, r);
            f = fractionMRG(&g, nbPoints, dim);
        }
        volRep[r] = cube * f;
        piRep[r]  = k * f;
    }

    double secs = omp_get_wtime() - start;

    /* ---------------- confidence intervals ---------------- */
    int     nu  = nbReplicates - 1;
    stats_t sp  = computeStats(piRep, nbReplicates);
    stats_t sv  = computeStats(volRep, nbReplicates);
    double  t95 = studentT(nu, T975, T975_NU40);
    double  t99 = studentT(nu, T995, T995_NU40);
    double  r95 = t95 * sp.se;          /* = t * sqrt(S^2(n) / n)     */
    double  r99 = t99 * sp.se;
    double  rv  = t95 * sv.se;
    char    label[32];

    /* theoretical standard deviation of one replicate: binomial      */
    double p       = M_PI / k;
    double sigmaTh = k * sqrt(p * (1.0 - p) / (double) nbPoints);

    printf("=== %s | %s (R = %g) | %d replicates x %lld points | %d threads ===\n",
           useMT ? "MT19937" : "MRG32k3a", dim == 3 ? "sphere 3D" : "disc 2D",
           RADIUS, nbReplicates, nbPoints, omp_get_max_threads());
    for (int r = 0; r < nbReplicates; r++)
        printf("  replicate %2d : V = %.10f | pi = %.10f\n",
               r + 1, volRep[r], piRep[r]);

    printf("\n  Volume (exact %s = %.10f)\n",
           dim == 3 ? "4/3 pi R^3" : "pi R^2", exactVol);
    printf("    mean V              = %.10f\n", sv.mean);
    printf("    95%% CI on V         = [%.10f ; %.10f]  R = %.3e\n",
           sv.mean - rv, sv.mean + rv, rv);

    printf("\n  PI\n");
    printf("    mean pi             = %.10f\n", sp.mean);
    printf("    S^2(n)              = %.3e\n", sp.s2);
    printf("    S(n)                = %.3e  (theory for 1 replicate: %.3e)\n",
           sp.s, sigmaTh);
    printf("    std error S/sqrt(n) = %.3e\n", sp.se);
    snprintf(label, sizeof label, "t(%d ; 0.975)", nu);
    printf("    %-20s= %.3f\n", label, t95);
    printf("    95%% CI radius R     = %.3e\n", r95);
    printf("    95%% CI              = [%.10f ; %.10f]  contains M_PI: %s\n",
           sp.mean - r95, sp.mean + r95,
           fabs(sp.mean - M_PI) <= r95 ? "yes" : "NO");
    snprintf(label, sizeof label, "t(%d ; 0.995)", nu);
    printf("    %-20s= %.3f\n", label, t99);
    printf("    99%% CI              = [%.10f ; %.10f]  contains M_PI: %s\n",
           sp.mean - r99, sp.mean + r99,
           fabs(sp.mean - M_PI) <= r99 ? "yes" : "NO");
    printf("    absolute error      = %.3e\n", fabs(sp.mean - M_PI));
    printf("    relative error      = %.3e\n", fabs(sp.mean - M_PI) / M_PI);
    printf("  time (wall)           = %.2f s\n", secs);

    free(piRep);
    free(volRep);
    return EXIT_SUCCESS;
}
