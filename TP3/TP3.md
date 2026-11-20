# TP3 – Fine Random Number Generators and Confidence Intervals

Two high-quality pseudo-random number generators, Mersenne Twister MT19937 and
L'Ecuyer's MRG32k3a, are checked against reference values and timed. They are
then used to estimate the volume of a sphere, and π, with independent Monte
Carlo replicates computed in parallel (OpenMP), and a confidence interval
around the mean (Student law).

## Principle

### Q1 – Generators

| Generator | Author | Output | Period | State |
|---|---|---|---|---|
| MT19937 | Matsumoto & Nishimura (1998) | `genrand_real1`: [0, 1] | 2¹⁹⁹³⁷ − 1 | 624 × 32 bits |
| MRG32k3a | L'Ecuyer (1999) | ]0, 1[ | ≈ 2¹⁹¹ | 6 × 32 bits |

Both are written in reentrant form: the state is a structure given as
parameter instead of global variables, so that each replicate (thread) owns
its own generator. The algorithms themselves are unchanged and are checked
against reference values at the start of `benchRNG` (see Results).

### Q2 – Volume of the sphere

Points `(x, y, z)` are drawn uniformly in the cube `[0,R]³` (one octant, by
symmetry). The fraction `f` of points with `x² + y² + z² ≤ R²` converges to
`(4/3 π R³ / 8) / R³ = π/6`, so

```
V = 4/3 π R³ ≈ 8 R³ f          π ≈ 6 f
```

With `dim = 2` the program computes the quarter disc of TP2 (`π ≈ 4 f`).
Standard deviation of one estimate of π with `N` points (binomial law):

```
3D : σ = 6 √(p(1−p)/N),  p = π/6   ≈ 3.00 / √N
2D : σ = 4 √(p(1−p)/N),  p = π/4   ≈ 1.64 / √N
```

The 3D estimator of π is 1.8 times noisier than the 2D one and needs 3 random
numbers per point instead of 2.

### Q2 – Confidence interval (appendix of the subject)

For `n` independent replicates `X1 … Xn` (each one with its own random stream):

```
X(n)   = Σ Xi / n                              mean
S²(n)  = Σ (Xi − X(n))² / (n − 1)              unbiased variance
S(n)/√n                                        standard error of the mean
R      = t(n−1, 1−α/2) × √(S²(n) / n)          confidence radius
CI     = [X(n) − R ; X(n) + R]                 at 1 − α = 95 %
```

`t(n−1, 0.975)` is read in Table 1 of the subject with **ν = n − 1 degrees of
freedom** (see Notes). A 99 % interval (α = 0.01, `t(n−1, 0.995)`) is printed
as well.

## Files

| File | Question | Description |
|---|---|---|
| `mt19937ar.h` | Q1-A | Mersenne Twister MT19937, reentrant version of the original `mt19937ar.c` |
| `mrg32k3a.h` | Q1-B | MRG32k3a from `combmrg2.c` (reentrant) + jump-ahead of 2¹²⁷ steps for independent streams |
| `benchRNG.c` | Q1-C | Reference checks of both generators, then timing of 10⁹ calls (+ `rand()` for comparison) |
| `RunBench` | Q1-C | Runs the benchmark compiled with `-O0`, `-O2` and `-O3 -march=native` |
| `sphereCI.c` | Q2 | Volume of the sphere and π, replicates in parallel (OpenMP), 95 % and 99 % confidence intervals |
| `RunCI` | Q2 | 10 and 40 replicates × 10⁹ points with both generators + two small sensitivity studies |
| `results/` | Q1, Q2 | Outputs of `RunBench` and `RunCI` |
| `openMPcomputePI.c` | – | Provided: 2D OpenMP example (xorshift generator) |
| `combmrg2.c` | – | Provided: MRG32k3a, MRG32k5a and MRG63k3a (reference only, not compilable as is) |

## Requirements

- Linux (or macOS / WSL on Windows)
- `gcc` with OpenMP support (`libgomp`, installed with gcc on most distributions)
- Optional: `htop` to watch the threads

## Build

```bash
gcc -O0 -Wall -o benchRNG_O0 benchRNG.c
gcc -O2 -Wall -o benchRNG_O2 benchRNG.c
gcc -O3 -march=native -Wall -o benchRNG_O3 benchRNG.c
gcc -O2 -Wall -fopenmp -o sphereCI sphereCI.c -lm
chmod +x RunBench RunCI
```

The generators are header files (`static inline` functions): there is nothing
else to compile or link. The binaries are listed in `.gitignore`.

## Usage

### Q1 – Generators: checks and performance

```bash
./benchRNG_O2          # reference checks + 10^9 calls per generator
./benchRNG_O2 1e8      # <nbCalls>
./RunBench             # -O0, -O2, -O3 one after the other -> results/bench_O*.txt
```

The three runs of `RunBench` are sequential on purpose: two benchmarks running
at the same time would share the CPU and distort the timings.

### Q2 – Confidence intervals

```bash
./sphereCI mt 10                   # MT19937, 10 replicates, 10^9 points, 3D
./sphereCI mrg 40                  # MRG32k3a, 40 replicates
./sphereCI mt 10 1e6               # <mt|mrg> <nbReplicates> [nbPoints] [dim]
./sphereCI mt 10 1e6 2             # 2D (quarter disc, as in TP2)
OMP_NUM_THREADS=4 ./sphereCI mt 10 # choose the number of threads
./RunCI                            # all the runs of the Results -> results/ (~8 min)
```

- `nbReplicates`: from 2 to 50 (the Student table is used up to ν = 40).
- `nbPoints` accepts the scientific notation (`1e9`).
- Each replicate owns its generator, whose stream only depends on the
  replicate number:
  - `mt`: `init_by_array` with Makoto Matsumoto's original key
    `{0x123, 0x234, 0x345, 0x456}`, extended with the replicate number as a
    5th word;
  - `mrg`: replicate *r* uses stream *r*, i.e. the default seed
    `{12345 × 6}` advanced by *r* × 2¹²⁷ steps.
- The results are therefore **identical whatever the number of threads**
  (checked with 1, 5 and 12 threads). The loop uses `schedule(dynamic, 1)`:
  a thread takes the next replicate as soon as it is free.

## Results

Machine: AMD Ryzen 5 5600H (6 cores / 12 threads), gcc 16.2.1.

### Q1 – Reference checks

| Check | Expected (reference) | Result |
|---|---|---|
| MT `init_by_array({0x123, 0x234, 0x345, 0x456})`, 5 first outputs | `1067595299 955945823 477289528 4107218783 4228976476` (`mt19937ar.out`) | OK |
| MT `init_genrand(5489)`, 10000th output | `4123659995` (C++11 standard, `std::mt19937`) | OK |
| MRG32k3a seed `{12345 × 6}`, first output | `0.1270111220` (exact integer computation) | OK |
| MRG32k3a jump of 2¹²⁷ steps | 2nd stream of RngStreams: `{3692455944, 1366884236, 2968912127, 335948734, 4161675175, 475798818}` | OK |

The checks pass at the three optimization levels, and the means of the draws
are identical (at the printed precision) between `-O0`, `-O2` and `-O3`:
optimization changes the speed, not the random sequences.

### Q1 – Performance, 10⁹ calls (ns per call = total time in s)

| Generator | `-O0` | `-O2` | `-O3 -march=native` |
|---|---|---|---|
| `rand()` (TP2, reference) | 7.43 | 7.40 | 7.40 |
| MT19937 `genrand_real1` | 9.15 | 4.30 | 2.05 |
| MRG32k3a | 48.74 | 28.71 | 27.14 |

- **MT19937** gains ×2.1 with `-O2` and ×4.5 with `-O3 -march=native`: once
  optimized it is by far the fastest, as announced by the subject.
- **`rand()`** does not change: it is in the C library, already compiled, so
  our options do not apply to it.
- **MRG32k3a** only gains ×1.7 and stays 5 to 13 times slower than MT: each
  call chains 2 divisions and 3 data-dependent branches (see Notes).

### Q2 – Confidence intervals, 10⁹ points per replicate (3D sphere)

σ of one replicate (theory): 9.48·10⁻⁵. Exact π = 3.1415926536.

| Generator | n | mean π | S(n) | t(n−1 ; 0.975) | 95 % radius R | 95 % CI | π ∈ CI | time |
|---|---|---|---|---|---|---|---|---|
| MT19937 | 10 | 3.1416330228 | 5.99·10⁻⁵ | 2.262 | 4.29·10⁻⁵ | [3.1415901724 ; 3.1416758732] | yes | 15.1 s |
| MT19937 | 40 | 3.1415992318 | 7.92·10⁻⁵ | 2.023 | 2.53·10⁻⁵ | [3.1415738960 ; 3.1416245677] | yes | 57.5 s |
| MRG32k3a | 10 | 3.1415737494 | 8.99·10⁻⁵ | 2.262 | 6.43·10⁻⁵ | [3.1415094622 ; 3.1416380366] | yes | 82.9 s |
| MRG32k3a | 40 | 3.1415963326 | 1.07·10⁻⁴ | 2.023 | 3.41·10⁻⁵ | [3.1415622124 ; 3.1416304529] | yes | 329.9 s |

Volume of the sphere with 40 replicates (exact 4/3 π = 4.1887902048):
MT 4.1887989758 ± 3.38·10⁻⁵, MRG32k3a 4.1887951102 ± 4.55·10⁻⁵ (95 %).

- **More replicates ⇒ smaller confidence radius.** From 10 to 40 replicates the
  radius is divided by 1.7 (MT) and 1.9 (MRG). With σ known, the expected
  factor is (2.262 / 2.023) × √(40/10) = 2.24 (theoretical radii 6.78·10⁻⁵
  and 3.03·10⁻⁵). The gap comes from S(n): with 10 replicates it is a noisy
  estimate of σ (5.99·10⁻⁵ for MT, low by chance); with 40 it gets closer to
  9.48·10⁻⁵. This uncertainty on S is exactly what the Student law accounts
  for (t = 2.262 instead of 1.960).
- The 4 intervals contain π. With 40 replicates (4·10¹⁰ points) the error on π
  is 6.6·10⁻⁶ (MT) and 3.7·10⁻⁶ (MRG32k3a).
- Both generators give statistically equivalent results; MRG32k3a costs 5.5
  to 5.7 times more time.

### Q2 – Effect of the number of replicates (MT, 10⁶ points each)

| n | t(n−1 ; 0.975) | 95 % radius R | R with σ theory | 95 % CI | π ∈ CI |
|---|---|---|---|---|---|
| 5 | 2.776 | 2.47·10⁻³ | 3.72·10⁻³ | [3.1383042930 ; 3.1432469070] | yes |
| 10 | 2.262 | 2.50·10⁻³ | 2.14·10⁻³ | [3.1397010252 ; 3.1446965748] | yes |
| 20 | 2.093 | 1.30·10⁻³ | 1.40·10⁻³ | [3.1405810524 ; 3.1431841476] | yes |
| 30 | 2.045 | 1.10·10⁻³ | 1.12·10⁻³ | [3.1404858257 ; 3.1426841743] | yes |
| 40 | 2.023 | 8.94·10⁻⁴ | 9.59·10⁻⁴ | [3.1404780850 ; 3.1422656150] | yes |
| 50 | 2.021 | 7.47·10⁻⁴ | 8.56·10⁻⁴ | [3.1405881746 ; 3.1420816654] | yes |

The radius decreases like `t(n−1) / √n`. With few replicates it fluctuates
(5 → 10: the S of the 5 first replicates is low by chance): below ~10
replicates, both `t` and the uncertainty on S are large. The replicates 1 … n
are the same streams for every n, so these intervals are nested samples, not
independent experiments.

### Q2 – Effect of the number of points per replicate (MT, 10 replicates)

| N | 95 % radius R | R with σ theory | 95 % CI | π ∈ CI |
|---|---|---|---|---|
| 10³ | 7.59·10⁻² | 6.78·10⁻² | [3.0614978493 ; 3.2133021507] | yes |
| 10⁴ | 2.24·10⁻² | 2.14·10⁻² | [3.1130336956 ; 3.1578063044] | yes |
| 10⁵ | 6.95·10⁻³ | 6.78·10⁻³ | [3.1329965676 ; 3.1468914324] | yes |
| 10⁶ | 2.50·10⁻³ | 2.14·10⁻³ | [3.1397010252 ; 3.1446965748] | yes |
| 10⁷ | 6.61·10⁻⁴ | 6.78·10⁻⁴ | [3.1407424139 ; 3.1420633861] | yes |
| 10⁸ | 1.57·10⁻⁴ | 2.14·10⁻⁴ | [3.1412230140 ; 3.1415368140] | **no** |
| 10⁹ | 4.29·10⁻⁵ | 6.78·10⁻⁵ | [3.1415901724 ; 3.1416758732] | yes |

- The radius is divided by about √10 ≈ 3.16 for every factor 10 on N: the slow
  `1/√N` convergence of Monte Carlo, one more correct digit costs 100 times
  more points.
- **The 10⁸ interval does not contain π, and this is expected**: a 95 %
  interval misses the true value in 5 % of the cases on average (1 time out
  of 20). Over the 15 distinct intervals of this TP, 1 misses. Using more
  points does not give "more certainty", only a narrower interval at the same
  95 % confidence.

### Q2 – Parallel speed-up (12 replicates × 5·10⁷ points, 3D)

```bash
for gen in mt mrg; do
    for t in 1 6 12; do OMP_NUM_THREADS=$t ./sphereCI $gen 12 5e7 | grep time; done
done
```

| Threads | MT19937 | speed-up | MRG32k3a | speed-up |
|---|---|---|---|---|
| 1 | 6.26 s | 1.0 | 50.24 s | 1.0 |
| 6 | 1.05 s | 6.0 | 8.36 s | 6.0 |
| 12 | 0.82 s | 7.6 | 4.18 s | 12.0 |

With one thread per physical core (6), the speed-up is ideal for both
generators: the replicates are fully independent. The 6 extra hardware
threads (SMT) bring little to MT, which already keeps the execution units
busy, but double the speed of MRG32k3a: its code mostly waits (divisions,
mispredicted branches), and the second hardware thread of each core uses the
idle units.

## Notes

- **Reading Table 1 of the subject**: the first column is the number of
  degrees of freedom ν = n − 1, not the number of replicates n (12.706 is t
  for ν = 1, i.e. 2 replicates; there is no interval with 1 replicate). With
  n = 10 replicates the right value is t(9) = 2.262, not 2.228.
  Between ν = 30 and ν = 40 the table is interpolated linearly in 1/ν
  (t(39) = 2.023, exact 2.0227); from ν = 40 to 49 the value of ν = 40
  (2.021) is used, slightly above the exact value, so the interval is
  slightly conservative. The number of replicates is limited to 50. The table
  value for ν = 8 (2.308) is kept as given, the exact value being 2.306.
- **α = 0.01**: the table for α = 0.01 asked in the appendix is in
  `sphereCI.c` (`T995`, from 63.657 for ν = 1 to 2.704 for ν = 40); the 99 %
  interval is printed after the 95 % one.
- **Independent streams**: for MRG32k3a, the streams are 2¹²⁷ draws apart
  (jump-ahead matrices of L'Ecuyer's RngStreams) and a replicate uses
  3·10⁹ ≈ 2³¹·⁵ draws: the streams cannot overlap, independence is guaranteed
  by construction. For MT19937 the replicates use different initialization
  keys: practically independent (period 2¹⁹⁹³⁷ − 1, `init_by_array` mixes the
  whole 624-word state), but non-overlapping is not formally guaranteed. The
  rigorous technique for MT is jump-ahead or pre-computed saved states
  (sequence splitting).
- **double, not float**: all the computations use `double` and the counters
  are `long long`. A `float` counter stops growing at 2²⁴ = 16 777 216
  (adding 1 no longer changes it), far below the 10⁹ points of a replicate,
  and with ~7 significant digits a `float` is at the limit of the 10⁻⁶
  precision reached on π.
- **Variance in two passes**: S²(n) is computed from the deviations to the
  mean, as in the appendix, and not with `Σ Xi²/n − X(n)²`, which subtracts
  two numbers ≈ 9.87 to get a variance ≈ 10⁻⁸: about half of the significant
  digits of a `double` would be lost (and all of them with a `float`).
- **Why MRG32k3a is slow here**: each call chains 2 double divisions and
  `double ↔ long` conversions, plus 3 branches that depend on the random
  values (`if (p < 0)` twice and the final combination), very often
  mispredicted. A test variant where these branches are replaced by masks
  (same recurrence, outputs identical bit for bit) runs at 13.6 ns per call
  instead of 28.1 ns at `-O2`. The provided code is kept unchanged in this TP.
- **`-march=native`**: the `-O3` binary uses the instructions of the build CPU
  (AVX2, FMA…) and may not run on another machine; `-O2` is the portable
  choice.
- **`M_PI`**: `sphereCI.c` defines `_DEFAULT_SOURCE` / `_USE_MATH_DEFINES`
  before including `<math.h>`, as in TP2.

## Troubleshooting

| Problem | Fix |
|---|---|
| `undefined reference to 'omp_get_wtime'` / `omp.h: No such file` | Compile with `-fopenmp` (install `libgomp` if needed) |
| `undefined reference to 'sqrt'` / `'pow'` | Add `-lm` at the end of the gcc command |
| `Permission non accordée` / `Permission denied` | `chmod +x RunBench RunCI` |
| `./benchRNG_O0: No such file or directory` | Build the binaries first (section Build), in the `TP3` directory |
| `Illegal instruction` with `benchRNG_O3` | Built with `-march=native` on another CPU: rebuild it on this machine |
| `Usage: ./sphereCI ...` | Generator `mt` or `mrg`, 2 to 50 replicates, dim 2 or 3 |
| `Reference check FAILED` | A generator was modified: compare with `mt19937ar.h` / `mrg32k3a.h` |
| `/bin/bash^M: bad interpreter` | Windows line endings: `sed -i 's/\r$//' RunBench RunCI` |
