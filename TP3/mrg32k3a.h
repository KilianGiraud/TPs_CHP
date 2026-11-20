/* ------------------------------------------------------------------ */
/* mrg32k3a.h - Combined Multiple Recursive Generator MRG32k3a        */
/*              (P. L'Ecuyer, Operations Research 47(1), 1999)        */
/*                                                                    */
/* MRG32k3a() is the code given on Moodle (combmrg2.c), unchanged      */
/* except that the 6 state variables s10 ... s22 are stored in a       */
/* structure given as parameter (reentrant => usable by threads) and   */
/* that the constants are prefixed (MRG_) to avoid name clashes.       */
/*                                                                    */
/* Independent streams: mrg_init_stream(g, k) starts the generator    */
/* k x 2^127 steps after the default seed, using the jump-ahead        */
/* matrices of L'Ecuyer's RngStreams package. Two streams cannot       */
/* overlap unless one of them draws more than 2^127 numbers.           */
/* ------------------------------------------------------------------ */
#ifndef MRG32K3A_H
#define MRG32K3A_H

#include <stdint.h>

#define MRG_NORM  2.328306549295728e-10     /* 1 / (m1 + 1)         */
#define MRG_M1    4294967087.0              /* 2^32 - 209           */
#define MRG_M2    4294944443.0              /* 2^32 - 22853         */
#define MRG_A12      1403580.0
#define MRG_A13N      810728.0
#define MRG_A21       527612.0
#define MRG_A23N     1370589.0

#define MRG_DEFAULT_SEED  12345.0           /* RngStreams default   */

typedef struct
{
    double s10, s11, s12;                   /* component 1          */
    double s20, s21, s22;                   /* component 2          */
} mrg_state;

/* Uniform double in ]0, 1[  (code from combmrg2.c)                   */
static inline double MRG32k3a(mrg_state *g)
{
    long   k;
    double p1, p2;
    /* Component 1 */
    p1 = MRG_A12 * g->s11 - MRG_A13N * g->s10;
    k = p1 / MRG_M1;   p1 -= k * MRG_M1;   if (p1 < 0.0) p1 += MRG_M1;
    g->s10 = g->s11;   g->s11 = g->s12;   g->s12 = p1;
    /* Component 2 */
    p2 = MRG_A21 * g->s22 - MRG_A23N * g->s20;
    k  = p2 / MRG_M2;  p2 -= k * MRG_M2;   if (p2 < 0.0) p2 += MRG_M2;
    g->s20 = g->s21;   g->s21 = g->s22;   g->s22 = p2;
    /* Combination */
    if (p1 <= p2) return ((p1 - p2 + MRG_M1) * MRG_NORM);
    else return ((p1 - p2) * MRG_NORM);
}

/* v = A.v mod m, computed exactly with 64-bit integers:               */
/* A[i][j] < 2^32 and v[j] < 2^32, so each product fits in 64 bits.    */
static inline void mrg_matvec(const uint64_t A[3][3], double v[3], uint64_t m)
{
    uint64_t x[3] = { (uint64_t) v[0], (uint64_t) v[1], (uint64_t) v[2] };

    for (int i = 0; i < 3; i++)
    {
        uint64_t s = 0;
        for (int j = 0; j < 3; j++)
            s += (A[i][j] * x[j]) % m;
        v[i] = (double) (s % m);
    }
}

/* Jump 2^127 steps ahead: start of the next stream                   */
static inline void mrg_jump_stream(mrg_state *g)
{
    /* A1^(2^127) mod m1 and A2^(2^127) mod m2 (RngStreams.c)         */
    static const uint64_t A1p127[3][3] = {
        { 2427906178U, 3580155704U,  949770784U },
        {  226153695U, 1230515664U, 3580155704U },
        { 1988835001U,  986791581U, 1230515664U } };
    static const uint64_t A2p127[3][3] = {
        { 1464411153U,  277697599U, 1610723613U },
        {   32183930U, 1464411153U, 1022607788U },
        { 2824425944U,   32183930U, 2093834863U } };

    double v1[3] = { g->s10, g->s11, g->s12 };
    double v2[3] = { g->s20, g->s21, g->s22 };

    mrg_matvec(A1p127, v1, (uint64_t) MRG_M1);
    mrg_matvec(A2p127, v2, (uint64_t) MRG_M2);

    g->s10 = v1[0];  g->s11 = v1[1];  g->s12 = v1[2];
    g->s20 = v2[0];  g->s21 = v2[1];  g->s22 = v2[2];
}

/* Default seed {12345 x 6}, then jump to the beginning of stream k    */
static inline void mrg_init_stream(mrg_state *g, int k)
{
    g->s10 = g->s11 = g->s12 = MRG_DEFAULT_SEED;
    g->s20 = g->s21 = g->s22 = MRG_DEFAULT_SEED;
    for (int i = 0; i < k; i++)
        mrg_jump_stream(g);
}

#endif /* MRG32K3A_H */
