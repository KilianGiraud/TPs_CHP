# TP4 – Beware of floating point arithmetic

Every code example of the subject is run and explained (one program per
section, in C or C++ as in the subject), with additional experiments to
understand the results. The "To go a little further" part is then done: a
stable softmax for the AI students and, for HPC, Kahan summation and other
compensated sums measured on huge random vectors (results and computing time),
and the reproducibility of a parallel reduction with OpenMP. This file is the
trace of the lab asked in question 4.

## Contents

| Section of the subject | Program | Main observation |
|---|---|---|
| How does it impact my computing? | `compare.c` | `0.7f` promoted to double is < 0.7 and 0.5 is exact: "only 1 is right"; the rounding direction depends on the value (0.6f > 0.6) |
| Don't test for equality, Relative errors | `equality.cpp` | `sqrt(10)²` = 10 + 1 ulp; an absolute tolerance is wrong for large and for tiny numbers |
| Worry about addition and subtraction (1) | `sumMul.c` | 10⁷ additions of 0.1f give 1 087 937 (+8.8 %), 10⁸ give 2 097 152 (−79 %); a product is rounded once |
| Worry about addition and subtraction (2) | `derivative.cpp` | output of the subject reproduced digit for digit; gcc's negative values explained (x87); complex step: exact |
| Floating point numbers have finite ranges | `ranges.cpp` | 2²⁴ + 1 and 2⁵³ + 1 are ties rounded to even; DBL_EPSILON is the gap after 1 |
| Use logarithms to avoid overflow and underflow | `logFactorial.cpp` | 200!/(190! 10!) = NaN directly, 10⁻¹³ relative error with the logs, exact with a product of ratios |
| Numeric operations don't always return numbers | `notNumbers.cpp` | inf, NaN, IEEE flags; `x - x == 0` also fails for inf; a NaN does not always propagate |
| Out of order execution | `associativity.c` | (a + b) + c ≠ a + (b + c); gcc fuses a·b + c into an FMA with `-mfma`; order of a long sum |
| Going further – AI | `softmax.c` | naive softmax gives NaN for logits ±1000; subtract the max; log-softmax |
| Going further – HPC, Q1 to Q3 | `summation.c`, `exactSum.c`, `sumBench.c` | Kahan exact to the last bit on 10⁷–10⁸ values, 3.9× slower on one core, free with 8 accumulators or with threads; variants up to the exact sum |
| Reduction phase (HPC) | `reproSum.c` | an OpenMP reduction changes with the number of threads AND from one run to the next; compensated or exact reductions do not |

## Files

| File | Section / question | Description |
|---|---|---|
| `compare.c` | How does it impact… | The 2 codes, the 2 fixes, IEEE fields and exact values of 0.1, 0.5, 0.6, 0.7 |
| `equality.cpp` | Equality, Relative errors | sqrt(10)², ulps, absolute vs relative tolerance, `almostEqual` |
| `sumMul.c` | Addition (1) | 10 × 0.1f vs 0.1f × 10, then n = 10 … 10⁸ additions |
| `derivative.cpp` | Subtraction (2) | Derivative of sin: forward, central, complex step; x87 reproduction |
| `ranges.cpp` | Finite ranges | 2²⁴, 2⁵³, DBL_EPSILON, gap between consecutive numbers, `float.h` limits |
| `logFactorial.cpp` | Logarithms | 200!/(190! 10!) computed 4 ways, binomial probability B(2000, ½) |
| `notNumbers.cpp` | Not numbers | inf, NaN, IEEE exception flags, comparisons and propagation of NaN |
| `associativity.c` | Out of order execution | Non-associativity, FMA, harmonic sum forward and backward |
| `softmax.c` | Going further (AI) | Naive and stable softmax, log-softmax |
| `summation.h`, `summation.c` | Going further (HPC) Q2, Q3 | naive, pairwise, Kahan, Neumaier, Sum2, Klein, long double, 8-accumulator and float versions |
| `exactSum.c` | Q2, Q3 | Exact sum (Shewchuk / Python `math.fsum`): the reference |
| `mt19937ar.h`, `mt19937ar.c` | Q2 | MT19937 of TP3 + `genrand_res53`, split into `.h` / `.c` |
| `sumBench.c` | Q2, Q3 | Accuracy and time of every algorithm on 4 data sets (double, float), condition number study |
| `reproSum.c` | Reduction phase | Order of the additions, OpenMP reductions with 1 … 12 threads, repeated runs |
| `Makefile` | – | Builds everything, including the 3 demonstration binaries |
| `RunAll` | – | Runs everything, outputs saved in `results/` |
| `results/` | – | Outputs of `./RunAll big` |
| `*.pdf` | – | Provided: subject, Goldberg (1991), C / C++ style guides |

## Requirements

- Linux (or WSL / macOS with gcc), `gcc` and `g++`, `make`
- OpenMP (`libgomp`, installed with gcc on most distributions) for `reproSum`
- An x86-64 processor: the `long double` (x87) remarks are specific to x86;
  `associativity_fma` needs FMA instructions (Intel since 2013, AMD since 2012)
- About 200 MB of free memory (1.3 GB for `./sumBench 1e8`)

## Build

```bash
make            # every program
make clean      # removes the binaries and object files
```

Equivalent commands:

```bash
gcc -O2 -Wall -Wextra -o compare       compare.c       -lm
gcc -O2 -Wall -Wextra -o sumMul        sumMul.c        -lm
gcc -O2 -Wall -Wextra -o associativity associativity.c -lm
gcc -O2 -Wall -Wextra -o softmax       softmax.c       -lm
g++ -O2 -Wall -Wextra -o equality      equality.cpp
g++ -O2 -Wall -Wextra -o derivative    derivative.cpp
g++ -O2 -Wall -Wextra -o ranges        ranges.cpp
g++ -O2 -Wall -Wextra -o logFactorial  logFactorial.cpp
g++ -O2 -Wall -Wextra -o notNumbers    notNumbers.cpp
gcc -O2 -Wall -Wextra -c summation.c exactSum.c mt19937ar.c sumBench.c
gcc -o sumBench sumBench.o summation.o exactSum.o mt19937ar.o -lm
gcc -O2 -Wall -Wextra -fopenmp -o reproSum reproSum.c summation.o exactSum.o mt19937ar.o -lm

# demonstration binaries, built on purpose with risky options
g++ -O2 -Wall -Wextra -ffast-math -o notNumbers_fast notNumbers.cpp
gcc -O2 -Wall -Wextra -mfma -o associativity_fma associativity.c -lm
gcc -O3 -ffast-math -Wall -Wextra -c -o summation_fast.o summation.c
gcc -o sumBench_fast sumBench.o summation_fast.o exactSum.o mt19937ar.o -lm
```

Everything is built with `-O2`, without `-ffast-math` and without
`-march=native`: strict IEEE 754 arithmetic, same results on any x86-64
machine. The `#warning` printed while compiling `summation_fast.o` is expected.

## Usage

```bash
./compare  ./equality  ./sumMul  ./derivative  ./ranges        # part 1: no argument
./logFactorial  ./notNumbers  ./associativity  ./softmax
./notNumbers_fast  ./associativity_fma          # same sources, risky options

./sumBench           # N = 1e7 values per data set (~8 s)
./sumBench 1e8       # N = 1e8 values (~70 s, 1.3 GB)
./sumBench_fast      # same, summation.c compiled with -O3 -ffast-math
./reproSum           # N = 1e7 values, OpenMP 1 .. 12 threads (~6 s)
./reproSum 1e8
./RunAll             # everything -> results/*.txt (~25 s)
./RunAll big         # + sumBench 1e8 (~1 min 35 s)
```

## Results and observations

Machine: AMD Ryzen 5 5600H (6 cores / 12 threads, Zen 3, AVX2 + FMA), 7 GB of
RAM, Fedora Linux, gcc / g++ 16.2.1. All the outputs are in `results/`.

### 1. How does it impact my computing? – `compare.c`

Both codes print **`only 1 is right`**.

| Constant | stored in a float (32 bits) | stored in a double (64 bits) |
|---|---|---|
| 0.7 | 0.699999988079071044921875 | 0.6999999999999999555910790149937383830547332763671875 |
| 0.5 | 0.5 | 0.5 |

- `a < 0.7`: `0.7` is a **double** constant, so `a` is promoted to double.
  The float is a coarser rounding of 0.7 (24 bits instead of 53) and it is
  rounded **down**: 0.69999998808 < 0.69999999999999995559 → true.
- `b < 0.5`: 0.5 = 2⁻¹ is exact in both formats → false.
- With `.7f` / `.5f`, or with double variables: `0 conditions are right`, the
  naive expected behavior.
- Figure 1 of the subject is checked: 0.1f = `0x3DCCCCCD`, 0.1 =
  `0x3FB999999999999A`, with the same sign / exponent / mantissa bits.
- **The result depends on the value.** 0.7f is rounded down, but 0.1f
  (0.100000001490…) and 0.6f (0.600000023841…) are rounded **up**: with 0.6
  instead of 0.7, `a < 0.6` would be false. The rule is not "a float is smaller
  than the double", but "the float and the double are two different roundings
  of the decimal constant".
- **How many digits to print?** `%0.18f` prints 0.699999999999999956, which is
  still not the stored value (it has 52 decimals). 17 significant digits
  identify a double uniquely (`%.17g`, `DBL_DECIMAL_DIG`), 9 a float (`%.9g`),
  and `%a` prints the exact value in hexadecimal (`0x1.6666666666666p-1`).

### 2. Don't test for equality, Relative errors – `equality.cpp`

- The code of the subject prints `-1.77636e-15` (default precision of `cout`:
  6 digits). The exact difference is −2⁻⁴⁹ = −1.7763568394002505·10⁻¹⁵:
  `sqrt(10)²` = 10.000000000000002 is **the next double after 10** (1 ulp).
  The value of the subject, `-1.778636e-015`, is probably a typo for
  `-1.776357e-015`.
- `FLT_EPSILON` = 2⁻²³ ≈ 1.19·10⁻⁷, `DBL_EPSILON` = 2⁻⁵² ≈ 2.22·10⁻¹⁶.

| Pair (expected, computed) | x − y | absolute test (< 10⁻¹²) | relative error | relative test (≤ 10⁻¹²) |
|---|---|---|---|---|
| 10, sqrt(10)² | −1.78·10⁻¹⁵ | equal | 1.8·10⁻¹⁶ | equal |
| 2·10¹⁰, sqrt(2·10¹⁰)² | −3.81·10⁻⁶ | **different** (wrong) | 1.9·10⁻¹⁶ | equal |
| 0.3, 0.1 + 0.2 | −5.55·10⁻¹⁷ | equal | 1.9·10⁻¹⁶ | equal |
| 10⁻¹⁵, 2·10⁻¹⁵ | −1.00·10⁻¹⁵ | **equal** (wrong: factor 2) | 1.0 | different |
| 0, sin(π) | −1.22·10⁻¹⁶ | equal | inf | **different** (wrong) |

- An absolute tolerance only works for numbers around 1 (too small for 10¹⁰,
  too large for 10⁻¹⁵). A relative tolerance works at any magnitude, except
  near 0 where the relative error is meaningless. A robust comparison combines
  both: `|a − b| ≤ max(rel × max(|a|, |b|), abs)` (`almostEqual`, the rule of
  Python's `math.isclose`). The values of `rel` and `abs` still depend on the
  expected accuracy of the computation, as the subject says.
- Relative error of the subject: result 99.5, expected 100 → 0.005.

### 3. Worry about addition (more than multiplication) – `sumMul.c`

Code of the subject:
`sum = 1.000000119209290, mul = 1.000000000000000, mul2 = 1.000000000000000`.

- The 10 additions of 0.1f make 10 roundings that accumulate: 1 + 2⁻²³ (1 ulp
  above 1). The product `f * 10` is rounded **once**: 1.0000000149 rounded to
  float is 1.
- `mul2` (`f * 10` passed to `printf`) is also computed in float on x86-64
  (10 is converted to float, `FLT_EVAL_METHOD = 0`). With `f * 10.0`, the
  product is a double: 1.000000014901161, the exact value of 0.1f × 10.

| n | float sum | rel. error | float product n × 0.1f | Kahan float sum |
|---|---|---|---|---|
| 10 | 1.0000001 | 1.0·10⁻⁷ | 1.0 | 1.0 |
| 10³ | 99.999 | 9.6·10⁻⁶ | 100.0 | 100.0 |
| 10⁵ | 9 998.557 | 1.4·10⁻⁴ | 10 000.0 | 10 000.0 |
| 10⁶ | 100 958.34 | 9.6·10⁻³ | 100 000.0 | 100 000.0 |
| 10⁷ | 1 087 937 | **8.8·10⁻²** | 1 000 000.0 | 1 000 000.0 |
| 10⁸ | 2 097 152 | **7.9·10⁻¹** | 10 000 000.0 | 10 000 000.0 |

- At 10⁸ additions, the float sum is **stuck at 2²¹ = 2 097 152**: there, two
  consecutive floats are 0.25 apart, 0.1 is less than half of this gap and
  `s + 0.1` is rounded back to `s` (absorption).
- The product and the Kahan sum (part 2) give the correctly rounded value for
  every n (the relative error 1.5·10⁻⁸ is only the rounding of the exact
  n × 0.1f to a float).
- In double, the same sum of 0.1 has a relative error of 1.9·10⁻⁹ for
  n = 10⁸, growing **linearly** with n: the term is always the same, so the
  rounding errors have the same sign and do not cancel out.

### 4. Worry about subtraction: derivative of sin – `derivative.cpp`

The program reproduces the output given in the subject digit for digit
(Visual C++ column), including the four `0` at the end:

| h | (sin(1+h) − sin(1))/h | error | central difference | error | complex step | error |
|---|---|---|---|---|---|---|
| 10⁻¹ | `0.4..............` | 4.3·10⁻² | `0.53.............` | 9.0·10⁻⁴ | `0.541............` | 9.0·10⁻⁴ |
| 10⁻⁵ | `0.5402...........` | 4.2·10⁻⁶ | `0.54030230585....` | **1.1·10⁻¹¹** | `0.54030230587....` | 9.0·10⁻¹² |
| 10⁻⁸ | `0.540302302......` | **3.0·10⁻⁹** | `0.540302308......` | 2.6·10⁻⁹ | `0.540302305868140` | **0** |
| 10⁻¹² | `0.54034..........` | 4.3·10⁻⁵ | `0.5402...........` | 1.2·10⁻⁵ | `0.540302305868140` | 0 |
| 10⁻¹⁵ | `0.55.............` | 1.5·10⁻² | `0.55.............` | 1.5·10⁻² | `0.540302305868140` | 0 |
| 10⁻¹⁶ … 10⁻¹⁹ | `0` | 0.54 | `0` (from 10⁻¹⁷) | 0.54 | `0.540302305868140` | 0 |

- Error ≈ h/2·|sin''| (truncation) + 2ε·|sin|/h (cancellation in the
  subtraction): minimum for h ≈ √ε ≈ 10⁻⁸, with an error of about 10⁻⁸, as
  observed. For h ≤ 10⁻¹⁶, 1 + h == 1 in double and the result is exactly 0.
- **Tricks that change the problem**, as advised by the subject: the central
  difference (f(x+h) − f(x−h))/2h has a truncation error in h², so the best h
  is ≈ ε^(1/3) ≈ 10⁻⁵ (error 1.1·10⁻¹¹). The **complex step**
  f'(x) ≈ Im f(x + ih)/h has **no subtraction at all**: no cancellation, and
  cos(1) is obtained to the last digit for every h ≤ 10⁻⁸ (`std::complex`).
- **Why "gcc on Linux" gave −0.017…, −0.17…, −1.7… and 17….** With gcc 16 on
  x86-64 we get 0, like Visual C++. The old values come from the **x87 FPU**
  used by 32-bit x86 compilers (`FLT_EVAL_METHOD = 2`): intermediate results
  stay in 80-bit registers (64-bit mantissa). When sin(1+h) is rounded to a
  double (stored in memory) while sin(1) stays in an 80-bit register, the
  numerator is the rounding error of sin(1) in double, −1.79·10⁻¹⁸, divided by
  h. Reproduced with `long double`: −0.0179, −0.179, −1.79, −17.9 (the 4th one
  is negative too). The same source gives different results depending on the
  compiler, its options and the hardware.

### 5. Floating point numbers have finite ranges – `ranges.cpp`

- `cout << f << " " << f+1` prints `1.67772e+07 1.67772e+07`: "16777216
  twice" needs `setprecision(8)` at least. Then `f + 1` = 16777216: `1` is
  converted to float and 2²⁴ + 1 is not a float (24-bit significand). It is a
  **tie** between 2²⁴ and 2²⁴ + 2, rounded to the even significand. `f + 1.0`
  is computed in double: 16777217.
- For x = 2⁵³: `(x + 1) − x = 0` (same tie in double), `(x + 2) − x = 2`,
  `(x + 3) − x = 4` (tie rounded to 2⁵³ + 4).
- `Sorry!` is printed. More precisely, `DBL_EPSILON` = 2⁻⁵² is the **gap
  between 1 and the next double**. 1 + ε/2 is exactly halfway between 1 and
  1 + ε and is rounded to even, i.e. 1. But 1 + 0.51·ε ≠ 1: the smallest e
  such that 1 + e ≠ 1 is just above ε/2, not ε. And 1 − ε/2 ≠ 1, because the
  gap below 1 is ε/2.

| Value | gap to the next float | gap to the next double |
|---|---|---|
| 1 | 1.19·10⁻⁷ | 2.22·10⁻¹⁶ |
| 10⁶ | 0.0625 | 1.16·10⁻¹⁰ |
| 2²⁴ = 16 777 216 | 2 | 3.73·10⁻⁹ |
| 10⁹ | 64 | 1.19·10⁻⁷ |
| 2⁵³ ≈ 9.0·10¹⁵ | 1.07·10⁹ | 2 |
| 10³⁰ | 7.6·10²² | 1.4·10¹⁴ |

The **absolute** precision decreases when the magnitude increases (the
relative precision is constant): a small number added to a large one is lost
(`sumMul.c`, harmonic sum of `associativity.c`, float sums of part 2).

| | float | double | long double (x87) |
|---|---|---|---|
| Bits of the significand | 24 | 53 | 64 |
| Digits always kept / needed for a round trip | 6 / 9 | 15 / 17 | 18 / 21 |
| EPSILON | 1.19·10⁻⁷ | 2.22·10⁻¹⁶ | 1.08·10⁻¹⁹ |
| Smallest normal / subnormal | 1.18·10⁻³⁸ / 1.40·10⁻⁴⁵ | 2.23·10⁻³⁰⁸ / 4.94·10⁻³²⁴ | 3.36·10⁻⁴⁹³² / 3.65·10⁻⁴⁹⁵¹ |
| Largest | 3.40·10³⁸ | 1.80·10³⁰⁸ | 1.19·10⁴⁹³² |
| Consecutive integers up to | 2²⁴ ≈ 1.7·10⁷ | 2⁵³ ≈ 9.0·10¹⁵ | 2⁶⁴ ≈ 1.8·10¹⁹ |

### 6. Use logarithms to avoid overflow and underflow – `logFactorial.cpp`

170! = 7.26·10³⁰⁶ but 171! = inf (`DBL_MAX` = 1.80·10³⁰⁸). Computed directly,
200!/(190! 10!) = inf/inf = **NaN** (printed `-nan`).

| Method | Result | Relative error |
|---|---|---|
| Exact (64-bit integers) | 22 451 004 309 013 280 | – |
| `exp(logFactorial(200) − logFactorial(190) − logFactorial(10))` | 22 451 004 309 011 012 | 1.0·10⁻¹³ |
| `exp(lgamma(201) − lgamma(191) − lgamma(11))` | 22 451 004 309 011 172 | 9.4·10⁻¹⁴ |
| Product of the ratios (190 + i)/i, i = 1..10 | 22 451 004 309 013 280 | 0 |

- The logs avoid the overflow but cost 3 digits: log 200! = 863.23 and
  log 190! = 810.48 have an absolute precision of about 10⁻¹³ (the ulp of
  863), the subtraction gives 37.65 with the same absolute error, and this
  error becomes the **relative** error of exp(37.65). Even `lgamma` cannot
  avoid it. When possible, a formula without huge intermediate values (product
  of the ratios) is better.
- Overflow **and** underflow in the same formula: P(X = 1000) for
  X ~ B(2000, ½) is C(2000, 1000) × 0.5²⁰⁰⁰ with C(2000, 1000) ≈ 2·10⁶⁰⁰ = inf
  and 0.5²⁰⁰⁰ ≈ 8.7·10⁻⁶⁰³ = 0, so inf × 0 = NaN. With `lgamma`:
  0.017839011145854369, relative error 2.7·10⁻¹⁵ (exact value
  0.01783901114585432, computed with Python fractions).

### 7. Numeric operations don't always return numbers – `notNumbers.cpp`

`2 * DBL_MAX` prints `inf`, `sqrt(-2)` prints `-nan` (glibc prints the sign:
the default NaN of x86 has its sign bit set).

| x | `x - x == 0` | `x == x` | `std::isnan(x)` |
|---|---|---|---|
| 1 | true | true | false |
| inf | **false** (inf − inf = NaN) | true | false |
| NaN | false | **false** | true |

So `x - x == 0` tests "x is finite", and `x == x` tests "x is not a NaN"
(infinities included), as stated by the subject.

IEEE 754 exception flags (`<cfenv>`), raised at run time:

| Operation | Result | Flags |
|---|---|---|
| 2 × DBL_MAX, exp(1000) | inf | OVERFLOW INEXACT |
| 1 / 0, log(0) | inf, −inf | DIVBYZERO |
| 0 / 0, sqrt(−2), log(−1), inf − inf, 0 × inf | NaN | INVALID |
| inf + 1 | inf | – |
| exp(−1000), DBL_TRUE_MIN / 2 | 0 | UNDERFLOW INEXACT |
| DBL_MIN / 2 | 1.11·10⁻³⁰⁸ (subnormal) | – (exact result: no underflow flag) |
| DBL_MIN / 3 | 7.42·10⁻³⁰⁹ (subnormal) | UNDERFLOW INEXACT |
| 0.1 + 0.2 | 0.30000000000000004 | INEXACT |
| 0.5 + 0.25 | 0.75 | – |

- Every comparison with a NaN is false (`NaN == NaN`, `NaN < 1`, `NaN >= 1`)
  except `!=`, so `!(x < 1)` is not the same as `x >= 1`.
- **A NaN does not always propagate**: `pow(NaN, 0) = 1`,
  `fmax(NaN, 1) = 1` (Goldberg §3.2.2: a max should ignore NaNs),
  `hypot(inf, NaN) = inf`, and `std::max(NaN, 1.0) = NaN` but
  `std::max(1.0, NaN) = 1`: the result depends on the order of the arguments.
- **With `-ffast-math` (`notNumbers_fast`)**: `x == x` is true for NaN,
  `std::isnan(NaN)` is false, `x - x == 0` is true for inf and NaN,
  `NaN * 0 = 0`, and the subnormals are flushed to 0 (DBL_MIN / 2 = 0,
  flush-to-zero set at program start). The compiler assumes there is no NaN
  and no inf: the tests of this section are removed.

### 8. Out of order execution: the operations are not associative – `associativity.c`

| Expression | Result |
|---|---|
| (10⁻³ + 1) − 1 | 0.000999999999999890 |
| 10⁻³ + (1 − 1) | 0.001000000000000000 |
| (2⁶⁰ − 2⁶⁰) + 1 | 1 |
| (2⁶⁰ + 1) − 2⁶⁰ | 0 |
| (0.1 + 0.2) + 0.3 | 0.6000000000000001 |
| 0.1 + (0.2 + 0.3) | 0.6 |
| (0.1 × 0.2) × 0.3 | 0.006000000000000001 |
| 0.1 × (0.2 × 0.3) | 0.006 |
| (10³⁰⁸ × 10) × 0.1 | inf |
| 10³⁰⁸ × (10 × 0.1) | 10³⁰⁸ |

- Products are not associative either, but the difference stays at 1 ulp,
  unless an intermediate result overflows.
- **FMA**: `fma(0.1, 10, −1)` = 5.55·10⁻¹⁷ (one rounding: the exact error of
  the product 0.1 × 10), while `p = 0.1 * 10; p − 1` = 0 (two roundings).
  Built with `-mfma` (`associativity_fma`) or `-march=native`, **gcc fuses
  `p = a * b; p - c` into an FMA by itself** (`-ffp-contract=fast` is the
  default in GNU C mode) and the same source prints 5.55·10⁻¹⁷.
  `-ffp-contract=off` or `-std=c17` forbid this contraction. FMA is more
  accurate, but the results change with the compilation options and the
  target processor.
- **Order of a long sum**: harmonic sum H(10⁷) in float. Forward (1 → n), the
  sum **stops growing at i = 2 097 152 = 2²¹**, when 1/i becomes smaller than
  half the gap between floats around 15.4: 15.4036827 (error 7.7 %). Backward
  (small terms first): 16.6860313 (5.6·10⁻⁴). Exact: 16.695311365859855.
- **About out-of-order execution in the processor.** The out-of-order engine
  of a core executes independent instructions in another order, but every
  instruction receives the same operands and the results are committed in
  program order: on one thread, a given binary computes bit-identical results
  at every run (checked: `sumBench` prints the same 17 digits at every run).
  It even makes some sums faster without changing them (`pairwise` is faster
  than `naive` because its blocks of 128 values are independent chains of
  additions). What really changes the order of the floating point operations
  is (1) the **compiler** (vectorization allowed by `-ffast-math`, FMA
  contraction), (2) the **parallel runtime** (partition of the loop depending
  on the number of threads, order in which the partial results are combined
  depending on the timing of the threads) and (3) libraries that choose their
  code path at run time (BLAS kernels per processor, for example). This
  dynamic non-reproducibility is measured with `reproSum` (part 2).

### 9. Going further – AI students: softmax – `softmax.c`

In float (32 bits), as in deep learning frameworks:

| Logits | Naive softmax | Stable softmax |
|---|---|---|
| (1, 2, 3, 4) | 0.032 0.087 0.237 0.644 | same |
| (1000, 1001, 1002, 1003) | NaN NaN NaN NaN (inf / inf) | 0.032 0.087 0.237 0.644 |
| (−1000, −1001, −1002, −1003) | NaN NaN NaN NaN (0 / 0) | 0.644 0.237 0.087 0.032 |
| (100, 50, 0, −50) | NaN 0 0 0 (exp(100) > FLT_MAX; fine in double) | 1 1.9·10⁻²² 3.8·10⁻⁴⁴ 0 |

- exp overflows for z > 88.7 in float (709.8 in double, only 11.1 in half
  precision, used for mixed-precision training) and underflows to 0 for very
  negative logits.
- Softmax does not change when the same constant is subtracted from every
  logit. With m = max zⱼ, every exp(zᵢ − m) ≤ 1 (no overflow) and the
  denominator is ≥ 1 (no 0/0): the logits (1000 … 1003) give exactly the
  result of (1 … 4).
- Cross-entropy loss −log pᵢ for the logits (0, 0, 0, −200) and the last
  class: naive p₃ = 0 (underflow), so the loss is **inf**; with the
  log-softmax zᵢ − m − log Σⱼ exp(zⱼ − m) (log-sum-exp trick) the loss is
  **201.099** (exact: 200 + ln 3 = 201.0986).

### 10. Going further – HPC: compensated sums

#### Q1 – Professor William Kahan

William Morton Kahan (born in 1933 in Toronto), professor emeritus of
mathematics and computer science at the University of California, Berkeley,
is the **main architect of the IEEE 754 standard (1985)**, hence his nickname
of "father of floating point". He received the **Turing Award in 1989** for
his fundamental contributions to numerical analysis. From 1976 he worked with
Intel on the 8087 floating point coprocessor, whose arithmetic was the basis
of the standard. He also wrote the "paranoia" program (1983), which tests the
floating point arithmetic of a machine, and numerically robust formulas
(area of needle-like triangles, roots of a quadratic equation, …).

What he proposed so that the rounding errors of the standard can be
compensated, or at least detected:

- **Compensated summation (1965)**: W. Kahan, "Further remarks on reducing
  truncation errors", Communications of the ACM 8(1). The error of every
  addition is computed exactly and given back to the next term (Theorem 8 of
  Goldberg).
- **Extended precision** in the standard (the 80-bit format of the x87, with
  a 64-bit significand), to compute intermediate results with more precision
  than the stored data.
- **Correct rounding** of +, −, ×, / and √, which makes error-free
  transformations such as TwoSum exact; **gradual underflow** (subnormal
  numbers: x − y = 0 ⇔ x = y); **inf, NaN and exception flags**, so that a
  computation can go on and the program can test what happened; **directed
  rounding modes**, also useful to estimate the sensitivity of a result to
  rounding errors by running it again with another rounding mode.

#### Q2 – Kahan summation vs regular summation: results and computing time – `sumBench.c`

Kahan's algorithm (Goldberg, Theorem 8), in `summation.c`:

```c
s = 0.0; c = 0.0;
for (i = 0; i < n; i++)
{
   y = x[i] - c;        /* term corrected by the bits lost before       */
   t = s + y;           /* the low-order bits of y are lost here        */
   c = (t - s) - y;     /* = -(lost bits), computed exactly             */
   s = t;
}
```

Error bounds (ε = 2⁻⁵³): naive ≤ (n − 1)·ε·Σ|xᵢ|, Kahan ≤ (2ε + O(nε²))·Σ|xᵢ|.
The error of Kahan's sum no longer grows with n.

**Protocol.** 4 data sets of N = 10⁷ doubles (80 MB, more than the 16 MB L3
cache) drawn with MT19937 (`genrand_res53`, 53 random bits per value), one
stream per data set: U[0,1]; 10^U[−10,10] (20 orders of magnitude); U[−1,1]
(cancellations); ±10^U[−10,10]. The reference is the exact sum rounded once
(`exactSum.c`, checked bit for bit on the 12 test cases of CPython's
`math.fsum`). Time = best of 5 runs. The computed values are identical at
every run.

Relative error, double, N = 10⁷ (0 = the correctly rounded sum, bit for bit):

| Algorithm | U[0,1] | 10^U[−10,10] | U[−1,1] (cond 5.8·10³) | ±10^U[−10,10] (cond 898) | ns / value |
|---|---|---|---|---|---|
| naive | 1.1·10⁻¹⁴ | **5.8·10⁻¹²** | 6.5·10⁻¹⁴ | 9.8·10⁻¹⁴ | 0.73 |
| naive, 8 accumulators | 5.4·10⁻¹⁵ | 7.1·10⁻¹³ | 3.5·10⁻¹⁴ | 9.3·10⁻¹⁵ | 0.40 |
| **Kahan** | **0** | **0** | **0** | **0** | **2.87** |
| Kahan, 8 accumulators | 0 | 0 | 0 | 0 | 0.67 |
| exact (fsum) | 0 | 0 | 0 | 0 | 5.3 – 16.4 |

With N = 10⁸ (`results/sumBench_1e8.txt`): naive 4.6·10⁻¹³ (U[0,1]) and
6.2·10⁻¹¹ (10^U[−10,10]), Kahan still exact, same times per value.

Single precision (the same values rounded to float), N = 10⁷:

| Algorithm | U[0,1] | 10^U[−10,10] | U[−1,1] | ±10^U[−10,10] | ns / value |
|---|---|---|---|---|---|
| naive (float accumulator) | 2.7·10⁻⁵ | **3.0·10⁻³** | 2.5·10⁻⁵ | 7.8·10⁻⁵ | 0.72 |
| Kahan (float accumulator) | 3.6·10⁻⁸ | 2.8·10⁻⁸ | 2.2·10⁻⁹ | 3.2·10⁻⁹ | 2.87 |
| naive in a double accumulator | 5.4·10⁻¹⁵ | 5.8·10⁻¹² | 2.6·10⁻¹⁶ | 5.2·10⁻¹⁴ | 0.71 |
| naive float, **N = 10⁸** | **16 777 216 instead of 5.0·10⁷ (−66 %)** | 3.3·10⁻² | 1.3·10⁻⁴ | 2.8·10⁻⁶ | 0.72 |

- **Results.** On every data set, Kahan gives the **correctly rounded sum**:
  0 error in double; in float, the float nearest to the exact sum (3.6·10⁻⁸ is
  only the rounding of 5 002 029.678 to 5 002 029.5). The naive sum loses 2 to
  4 digits in double and up to the whole result in float: with N = 10⁸ values
  U[0,1], the float sum stops at 2²⁴ = 16 777 216, where two floats are 2
  apart and every term smaller than 1 is absorbed (section "finite ranges" of
  the subject, inside a real reduction).
- **The range of the values matters more than N.** With 20 orders of
  magnitude, half of the values are smaller than 1 and are partly absorbed by
  a running sum of about 10¹⁵. The lost bits all have the same sign, so the
  error is systematic: 500 times larger than with U[0,1] (5.8·10⁻¹² vs
  1.1·10⁻¹⁴). The computed sum is below the exact one.
- **Computing time.** Kahan is **3.9 times slower** than the naive loop (2.87
  vs 0.73 ns per value), although it only does 4 additions instead of 1: the
  loop is limited by the **latency** of a floating point addition (3 cycles on
  Zen 3, i.e. ≈ 0.72 ns at 4.2 GHz). The naive loop is a chain of dependent
  additions, 1 per value; in Kahan's loop the compensation chains 4 dependent
  operations per value (y → t → t − s → c): 4 × 0.72 = 2.88 ns, exactly what
  is measured.
- **Break the chain, keep the compensation.** With 8 independent accumulators
  (accumulator j sums the values i ≡ j mod 8), **Kahan takes 0.67 ns per
  value: faster than the naive sequential loop**, with the same accuracy. The
  naive sum with 8 accumulators takes 0.40 ns (limited by the memory
  bandwidth of one core). In parallel, the threads play the same role (see
  `reproSum`: same time for the naive and Kahan sums with 12 threads).
- For float data, the cheapest fix is a **double accumulator**: free (0.71
  ns) and accurate to 10⁻¹⁴–10⁻¹⁶ (the "simple way" of Goldberg §3.2.3).
  Kahan in float gives the float nearest to the exact sum.

#### The `-ffast-math` trap – `sumBench_fast`

`summation.c` compiled with `-O3 -ffast-math` (this file only; the reference
`exactSum.c` refuses to compile with this option):

| Algorithm | U[0,1], `-O2` | U[0,1], `-O3 -ffast-math` | ns / value, `-O2` → fast-math |
|---|---|---|---|
| naive | 1.1·10⁻¹⁴ | 3.4·10⁻¹⁵ | 0.73 → 0.38 |
| Kahan | 0 | **3.4·10⁻¹⁵ (= naive)** | 2.87 → 0.38 |
| Neumaier, Sum2, Klein | 0 | **3.4·10⁻¹⁵ (= naive)** | 0.95 – 1.36 → 0.38 |
| Kahan float, 10^U[−10,10] | 2.8·10⁻⁸ | **7.6·10⁻⁴ (= naive float)** | 2.86 → 0.19 |

- With `-ffast-math` the compiler may use the algebra of real numbers:
  c = ((s + y) − s) − y = 0, the compensation disappears and **every
  compensated sum becomes the naive sum, bit for bit** (Kahan with 8
  accumulators becomes the naive sum with 8 accumulators; the same holds on
  every line of the condition number study). Goldberg already
  warned about it in 1991 (§3.2.3): "an optimizer that believed floating-point
  arithmetic obeyed the laws of algebra would conclude that C = 0, rendering
  the algorithm completely useless".
- The naive sum is vectorized (several partial sums in parallel): it is faster
  and here even slightly more accurate, but it is **another result** (another
  order); the tiny sum 1 + 10¹⁰⁰ + 1 − 10¹⁰⁰ even gives 2 instead of 0.
  `-Ofast` includes `-ffast-math`, and Intel compilers use a "fast" floating
  point model by default: compensated sums need `-fp-model=precise` there.

#### Q3 – Other variants of compensated sums

| Algorithm | Author(s), year | Idea | Error (S = Σ\|xᵢ\|) | Time here |
|---|---|---|---|---|
| naive | – | s += x | ≤ (n − 1)·ε·S | 0.73 ns |
| pairwise (cascade) | used by NumPy, Julia | sum of the two halves, recursively (naive below 128 values) | ≈ ε·log₂n·S | 0.44 ns |
| **Kahan** | Kahan, 1965 | the lost bits are given back to the next term | ≤ (2ε + O(nε²))·S | 2.87 ns |
| **Kahan–Babuška–Neumaier** | Babuška 1969, Neumaier 1974 | the error of s + x is computed from the larger operand, the errors are summed apart | ≤ ε·\|Σx\| + O(n²ε²)·S | 0.95 ns |
| **Sum2 (TwoSum)** | Møller 1965, Knuth 1969 (TwoSum); Ogita, Rump, Oishi 2005 | exact error of each addition without any test (6 operations) | as the naive sum computed with twice the precision: ≤ ε·\|Σx\| + γ²ₙ₋₁·S | 1.19 ns |
| **Klein** | Klein, 2006 | 2nd-order Kahan–Babuška: the errors of the compensation are compensated | observed ~10⁴ times more robust than Neumaier | 1.35 ns |
| SumK, AccSum | Ogita, Rump, Oishi 2005–2008 | K − 1 error-free passes; faithful rounding | as K times the precision | – |
| doubly compensated | Priest, 1992 | sort by decreasing \|x\|, then double compensation | ≤ 2ε·\|Σx\| | sort O(n log n) |
| wider accumulator | – | long double (x87), double for float data, double-double, `__float128` | depends on the format | 1.57 ns (long double) |
| **exact (fsum)** | Shewchuk 1997, Python `math.fsum` | exact list of non-overlapping partials, rounded once | correctly rounded | 5 – 16 ns |
| long accumulator, reproducible sums | Kulisch; ReproBLAS (Demmel, Nguyen 2013); ExBLAS (Collange, Defour, Graillat, Iakymchuk 2015) | fixed-point accumulator of ~4200 bits, or pre-rounding of all the values to a common grid | exact, or bit for bit identical whatever the order | – |

The algorithms in bold, the long double, the pairwise and the exact sums are
implemented in `summation.c` / `exactSum.c` and measured above.

**Tiny example (Neumaier's case)**: 1 + 10¹⁰⁰ + 1 − 10¹⁰⁰ = 2. naive,
pairwise, Kahan and long double give **0**; Neumaier, Sum2, Klein and exact
give **2**. Kahan fails because the new term (10¹⁰⁰) is larger than the
running sum: (t − s) − y is then no longer the exact error.

**Ill-conditioned sums.** N = 10⁶ values: pairs (v, −v) with |v| up to 2^E,
which cancel exactly, plus small values U[0,1); the exact sum is known.
cond = Σ|xᵢ| / |Σxᵢ|; relative error (≥ 1: no correct digit at all):

| cond | naive | pairwise | Kahan | Neumaier = Sum2 | Klein | long double | exact |
|---|---|---|---|---|---|---|---|
| 4 | 3.4·10⁻¹⁵ | 1.2·10⁻¹⁶ | 0 | 0 | 0 | 1.2·10⁻¹⁶ | 0 |
| 3.0·10⁵ | 5.5·10⁻¹² | 3.9·10⁻¹³ | 1.4·10⁻¹⁴ | 0 | 0 | 6.7·10⁻¹⁵ | 0 |
| 1.6·10¹¹ | 6.3·10⁻⁶ | 1.9·10⁻⁷ | 3.9·10⁻⁸ | 0 | 0 | 3.5·10⁻⁹ | 0 |
| 1.3·10¹⁴ | **0.90** | 1.7·10⁻⁴ | 1.9·10⁻⁵ | 2.9·10⁻¹⁴ | 0 | 5.1·10⁻⁶ | 0 |
| 1.1·10¹⁷ | 9.3 | 1.1 | **8.0·10⁻²** | 3.6·10⁻¹³ | 2.3·10⁻¹⁶ | 0.59 | 0 |
| 1.0·10²⁰ | 7.8·10³ | 2.7·10² | 73 | 1.1·10⁻¹⁰ | 3.2·10⁻¹³ | 7.4 | 0 |
| 8.3·10²⁵ | 1.3·10⁹ | 1.4·10⁸ | 4.3·10⁷ | 5.2·10⁻⁵ | 1.7·10⁻⁸ | 1.9·10⁶ | 0 |
| 7.3·10²⁸ | 1.0·10¹² | 3.6·10¹⁰ | 5.0·10¹⁰ | **0.98** | 4.3·10⁻⁶ | 6.8·10⁸ | 0 |
| 7.0·10³¹ | 7.6·10¹⁵ | 2.2·10¹⁴ | 2.1·10¹³ | 35 | **4.9·10⁻²** | 4.1·10¹² | 0 |
| 4.6·10⁵⁸ | 6.9·10⁴¹ | 2.1·10⁴¹ | 3.7·10⁴⁰ | 9.1·10²⁷ | 1.0 | 6.0·10³⁸ | **0** |

- The relative error of a sum is governed by **ε × cond**, much more than by
  N. Every algorithm has a breaking point: the naive sum and the long double
  around 10¹⁴–10¹⁷; Kahan around 10¹⁶ (its bound is relative to Σ|xᵢ|: it is
  a first-order method); Neumaier and Sum2 around 10²⁸ (twice the working
  precision: the error grows like ε²·cond); Klein around 10³²; only the exact
  sum always works.
- **Neumaier and Sum2 give the same results bit for bit** (on every data
  set): both compute the exact error of every addition (Fast2Sum on the larger
  operand for Neumaier, TwoSum for Sum2) and add the same errors in the same
  order. Sum2 replaces the test by 2 additions (no branch).
- Kahan is enough for well-conditioned sums (all the data sets of the first
  table). Neumaier and Sum2 are cheaper (0.95 and 1.19 ns vs 2.87 ns: their
  compensation is outside the chain of dependent additions) and more robust:
  for a sequential sum they are the best compromise.

#### Reduction phase and reproducibility – `reproSum.c`

The same 10⁷ values U[−1,1] (exact sum −859.10864839042449, 1 ulp =
1.1·10⁻¹³).

**1. Sequential sums, the same values in 6 orders** (forward, backward,
increasing |x|, decreasing |x|, 2 random permutations): the naive sum gives
**6 different results** (from 2.3·10⁻¹³ with increasing |x| to 2.9·10⁻¹⁰ with
decreasing |x|); Kahan, Neumaier and exact give **the same correctly rounded
result** in every order.

**2. OpenMP, 1 to 12 threads** (sum − exact, then time):

| Threads | naive `reduction(+)` | Kahan per thread, exact combination | exact (partials per thread, merged) | time naive | time Kahan | time exact |
|---|---|---|---|---|---|---|
| 1 | 5.6·10⁻¹¹ | 0 | 0 | 7.3 ms | 28.6 ms | 53.4 ms |
| 2 | 7.5·10⁻¹¹ | 0 | 0 | 3.8 ms | 14.5 ms | 27.2 ms |
| 3 | −1.2·10⁻¹⁰ | 0 | 0 | 2.6 ms | 9.7 ms | 18.7 ms |
| 4 | 2.9·10⁻¹¹ | 0 | 0 | 2.0 ms | 7.3 ms | 14.1 ms |
| 6 | 3.0·10⁻¹² | 0 | 0 | 2.0 ms | 4.9 ms | 9.8 ms |
| 8 | 7.1·10⁻¹¹ | 0 | 0 | 2.1 ms | 3.7 ms | 7.6 ms |
| 12 | −2.7·10⁻¹¹ | 0 | 0 | 2.6 ms | **2.5 ms** | 5.3 ms |
| distinct results | **7** | 1 | 1 | | | |

**3. Same binary, same data, 12 threads, 10 runs**: the naive
`reduction(+)` with `schedule(static)` gives **2 to 5 different results** over
10 runs (it depends on the run); with `schedule(dynamic)`, **10 different
results out of 10**; Kahan per thread and the exact sum always give the same
result.

- The result of a naive parallel reduction depends on the number of threads
  (partition of the loop) **and on the run**: at the end of the loop, each
  thread adds its partial sum to the shared variable when it arrives, so the
  order of these last additions depends on the timing of the threads (and with
  `dynamic`, each thread even gets different blocks). This is the dynamic
  non-reproducibility described in the subject; it comes from the scheduling
  of the threads.
- Compensating inside each thread and combining the partial results exactly,
  in a fixed order (thread number), makes the result independent of the
  number of threads, here bit for bit. The exact parallel sum (exact partials
  per thread, merged exactly) is reproducible **by construction**, whatever
  the partition and the order.
- **Cost.** The naive reduction is limited by the memory bandwidth from 4
  threads on (80 MB in 2 ms ≈ 40 GB/s). Kahan, 4 times slower on 1 thread,
  catches up: **with 12 threads, Kahan takes the same time as the naive sum
  (2.5 ms)**. In HPC, where reductions are memory-bound, the compensation is
  almost free; the exact sum costs twice as much (5.3 ms).

## Notes

- **Codes of the subject**: the curly quotes `“ ”` of the PDF are replaced
  by `"`, the missing `#include`, declarations and `\n` are added,
  `using namespace std` is replaced by explicit `std::`; otherwise the codes
  are unchanged.
- **Probable typos in the subject**: `-1.778636e-015` (exact value
  −1.7763568·10⁻¹⁵); "prints 16777216 twice" needs `setprecision(8)` or more
  (`cout` prints `1.67772e+07` by default); the 4th "gcc" value is −17.9,
  negative like the others.
- **Compilation options**: everything is built with `-O2`, without
  `-ffast-math` (it breaks compensated sums and NaN tests) and without
  `-march=native` (FMA contraction, non-portable binaries). All the values of
  this file are those of strict IEEE 754 arithmetic, identical on any x86-64
  machine, except the times and the run-dependent results of `reproSum`
  (section 3).
- **Random values with 53 bits**: `genrand_res53` (MT19937 of TP3, split into
  `mt19937ar.h` / `mt19937ar.c`). With `genrand_real1` (32 random bits) the
  values would be multiples of 2⁻³²: sums of up to about 2²¹ of them would be
  exact in double, which would hide the rounding errors studied here.
- **The reference**: `exactSum.c` (Shewchuk's algorithm, as in Python's
  `math.fsum`) gives the exact sum rounded once. It is checked bit for bit on
  the 12 test cases of CPython (`test_fsum`) at the start of `sumBench`, and
  refuses to compile with `-ffast-math` (`#error`).
- **ISIMA style guide** (provided with the lab): every file starts with a
  comment block (name and role), every function is preceded by a block (role,
  input, output), user types end with `_t`, constants are in upper case
  (`const` in C++), declarations are aligned, the indentation is 3 spaces,
  lines are at most 80 characters, and the `.h` files only contain
  declarations (unlike TP3, where the generators were `static inline`
  functions in headers).
- **`long double`**: on x86-64 it is the x87 80-bit format (64-bit
  significand). It is twice as slow as double here and not portable (it is a
  plain double on Windows/MSVC and on Apple ARM, a software 128-bit format on
  Linux ARM64).

## Troubleshooting

| Problem | Fix |
|---|---|
| `make: command not found` | Install `make`, or use the commands of the section Build |
| `omp.h: No such file` / `undefined reference to 'omp_get_wtime'` | Compile `reproSum` with `-fopenmp` (install `libgomp` if needed) |
| `warning: #warning "summation.c built with -ffast-math…"` | Expected: it is `summation_fast.o`, built on purpose for `sumBench_fast` |
| `#error "exactSum.c must not be built with -ffast-math"` | `-ffast-math` or `-Ofast` was added to `CFLAGS`: remove it |
| `Permission denied` on `RunAll` | `chmod +x RunAll` |
| `malloc: Cannot allocate memory` | N too large for the memory: `./sumBench 1e7` (80 MB per array) |
| `Illegal instruction` with `associativity_fma` | Processor without FMA: skip this demonstration |
| Times different from this file | Processor frequency, other running programs: the times vary; the computed values must be identical, except `reproSum` section 3 (run-dependent on purpose) |
| `/bin/bash^M: bad interpreter` | Windows line endings: `sed -i 's/\r$//' RunAll` |
