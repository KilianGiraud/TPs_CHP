# TP2 – Monte Carlo Simulation and Parallel Replications

Estimation of π with the Monte Carlo method, first in a single run, then with
independent sequential replicates, and finally with replicates launched in
parallel on an SMP machine using Unix multitasking (SPMD model).

## Principle

Points `(xr, yr)` are drawn uniformly in the unit square `[0,1]²`. The fraction
of points falling inside the quarter disc `x² + y² ≤ 1` converges to `π/4`, so

```
π ≈ 4 × inDisk / nbPoints
```

Each point needs two random numbers. The standard deviation of one estimate is
`σ = √(π(4−π)/n) ≈ 1.64/√n`: every additional correct digit costs 100 times
more points.

## Files

| File | Question | Description |
|---|---|---|
| `simPi.c` | Q1 | Single estimation of π with `rand()`, tests with 10³, 10⁶, 10⁹ points and a precision study |
| `replicates.c` | Q2 | 10 independent replicates run sequentially with Mersenne Twister, mean, absolute and relative error |
| `simPI.c` | Q3 | One replication, seed given on the command line, result printed on stdout |
| `Sim10times` | Q3 | Shell script launching 10 `simPI` processes in parallel with different seeds |

## Requirements

- Linux (or macOS / WSL on Windows)
- `gcc`
- Optional: `htop` or `top` to watch the processes

## Build

```bash
gcc -O2 -Wall -o simPi      simPi.c      -lm
gcc -O2 -Wall -o replicates replicates.c -lm
gcc -O2 -Wall -o simPI      simPI.c
chmod +x Sim10times
```

`-O2` matters: without optimization the 10⁹-point runs are several times slower.

## Usage

### Q1 – Single estimation

```bash
./simPi            # tests with 1e3, 1e6, 1e9 points + precision study
./simPi 1000000    # single run with the given number of points
```

### Q2 – Sequential replicates

```bash
./replicates              # 10 replicates for 1e3, 1e6 and 1e9 points
./replicates 10 1000000   # <nbReplicates> <nbPoints>
```

The generator is initialized **once** with the original key from Makoto
Matsumoto's reference code (`{0x123, 0x234, 0x345, 0x456}`) and is never
re-seeded between replicates: each replicate uses the next part of the same
stream, which guarantees independence.

### Q3 – Parallel replications (SPMD)

```bash
./simPI 19937              # one replication, 1e9 points by default
./simPI 19937 1000000      # <seed> [nbPoints]
./simPI 19937 > PI1        # result redirected to a file

./Sim10times               # 10 replications in parallel -> PI1 … PI10
time ./Sim10times          # measure the parallel execution time
```

Each line of `Sim10times` runs `./simPI <seed> > PIi &` in the background; Linux
spreads the 10 processes over the available cores. A final `wait` makes the
script return only when all replications are finished.

Display the results and a quick mean:

```bash
for i in $(seq 1 10); do echo "PI$i: $(cat PI$i)"; done
cat PI* | awk '{s+=$1; n++} END {printf "mean = %.10f\n", s/n}'
```

## Results

### Q1 – Precision vs number of points (theory)

| Target precision | n (1σ, ~68 %) | n (2σ, ~95 %) |
|---|---|---|
| 10⁻² (3.14) | ≈ 2.7·10⁴ | ≈ 1.1·10⁵ |
| 10⁻³ | ≈ 2.7·10⁶ | ≈ 1.1·10⁷ |
| 10⁻⁴ (3.1415) | ≈ 2.7·10⁸ | ≈ 1.1·10⁹ |

### Q2 – 10 sequential replicates (Mersenne Twister)

| nbPoints | meanPI | absolute error | relative error |
|---|---|---|---|
| 10³ | 3.1480000 | 6.4·10⁻³ | 2.0·10⁻³ |
| 10⁶ | 3.1412544 | 3.4·10⁻⁴ | 1.1·10⁻⁴ |
| 10⁹ | 3.1415644 | 2.8·10⁻⁵ | 9.0·10⁻⁶ |

### Q3 – 10 parallel replications, 10⁹ points each

| Seed | File | Estimated π |
|---|---|---|
| 19937 | PI1 | 3.1415328240 |
| 52149 | PI2 | 3.1416400000 |
| 31415 | PI3 | 3.1416069640 |
| 27182 | PI4 | 3.1416070920 |
| 14142 | PI5 | 3.1416196400 |
| 17320 | PI6 | 3.1416013360 |
| 44257 | PI7 | 3.1415138560 |
| 70913 | PI8 | 3.1415850600 |
| 63113 | PI9 | 3.1416275640 |
| 86226 | PI10 | 3.1416601640 |

| Statistic | Value |
|---|---|
| mean π | 3.1415994500 |
| absolute error | 6.8·10⁻⁶ |
| relative error | 2.2·10⁻⁶ |
| standard deviation | ≈ 4.5·10⁻⁵ (theory: 5.2·10⁻⁵) |
| standard error of the mean | ≈ 1.4·10⁻⁵ |

The results are bit-for-bit reproducible: the same seeds give the same values
on any machine.

## Notes

- **`rand()` is only used in Q1** as a first approach. It is statistically weak
  and slow (glibc takes a lock on each call); Q2 and Q3 use Mersenne Twister
  MT19937 (period 2¹⁹⁹³⁷ − 1).
- **Parallel speed-up depends on the number of logical cores** (`nproc`). With
  10 or more cores, the 10 replications take about the time of one; with fewer,
  the total time is roughly `10 / nproc` times a single run.
- **Independence of the parallel streams**: different seeds with
  `init_genrand` give practically independent streams, but no overlap is not
  formally guaranteed. The rigorous approach is to use pre-computed, widely
  spaced generator states (jump-ahead), one per process.
- **`M_PI` undefined**: `M_PI` is POSIX, not ISO C. `simPi.c` defines a fallback,
  `replicates.c` defines `_DEFAULT_SOURCE` / `_USE_MATH_DEFINES` before
  including `<math.h>`.

## Troubleshooting

| Problem | Fix |
|---|---|
| `Permission non accordée` / `Permission denied` | `chmod +x Sim10times` |
| Still denied after `chmod` | The directory may be mounted `noexec` (check with `findmnt -T .`): copy the files to a local directory |
| `./simPI: No such file or directory` | Compile `simPI.c` in the same directory as the script |
| `/bin/bash^M: bad interpreter` | Windows line endings: `sed -i 's/\r$//' Sim10times` |