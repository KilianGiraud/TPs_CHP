# -------------------------------------------------------------------- #
# Makefile          TP4 - Beware of floating point arithmetic          #
#                                                                      #
#   make            builds every program                               #
#   make run        builds, then runs everything (RunAll) -> results/  #
#   make clean      removes the binaries and the object files          #
#                                                                      #
# Every program is built with -O2, WITHOUT -ffast-math and WITHOUT     #
# -march=native: strict IEEE 754 semantics, same results on any x86-64 #
# machine. Three binaries are built ON PURPOSE with risky options:     #
#   sumBench_fast      summation.c compiled with -O3 -ffast-math: the  #
#                      compensations of Kahan & co. are removed. The   #
#                      link is done without -ffast-math so that only   #
#                      summation.c is affected (no flush-to-zero).     #
#   notNumbers_fast    -ffast-math: NaN tests removed, subnormals = 0  #
#   associativity_fma  -mfma: gcc fuses p = a * b; p - c into an FMA   #
# -------------------------------------------------------------------- #

CC       = gcc
CXX      = g++
CFLAGS   = -O2 -Wall -Wextra
CXXFLAGS = -O2 -Wall -Wextra
FASTMATH = -O3 -ffast-math
OMPFLAGS = -fopenmp
LDLIBS   = -lm

C_PROGS   = compare sumMul associativity softmax
CXX_PROGS = equality derivative ranges logFactorial notNumbers
SUM_PROGS = sumBench sumBench_fast reproSum
DEMOS     = notNumbers_fast associativity_fma
SUM_OBJS  = exactSum.o mt19937ar.o

all: $(C_PROGS) $(CXX_PROGS) $(SUM_PROGS) $(DEMOS)

# ----------- part 1: one program per section of the subject ----------
$(C_PROGS): %: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

$(CXX_PROGS): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

notNumbers_fast: notNumbers.cpp
	$(CXX) $(CXXFLAGS) -ffast-math -o $@ $<

associativity_fma: associativity.c
	$(CC) $(CFLAGS) -mfma -o $@ $< $(LDLIBS)

# ----------- part 2: summation algorithms ----------------------------
summation.o: summation.c summation.h
	$(CC) $(CFLAGS) -c -o $@ $<

summation_fast.o: summation.c summation.h
	$(CC) $(FASTMATH) -Wall -Wextra -c -o $@ $<

exactSum.o: exactSum.c summation.h
	$(CC) $(CFLAGS) -c -o $@ $<

mt19937ar.o: mt19937ar.c mt19937ar.h
	$(CC) $(CFLAGS) -c -o $@ $<

sumBench.o: sumBench.c summation.h mt19937ar.h
	$(CC) $(CFLAGS) -c -o $@ $<

sumBench: sumBench.o summation.o $(SUM_OBJS)
	$(CC) -o $@ $^ $(LDLIBS)

sumBench_fast: sumBench.o summation_fast.o $(SUM_OBJS)
	$(CC) -o $@ $^ $(LDLIBS)

reproSum: reproSum.c summation.o $(SUM_OBJS) summation.h mt19937ar.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -o $@ $< summation.o $(SUM_OBJS) $(LDLIBS)

# ----------- utilities -----------------------------------------------
run: all
	./RunAll

clean:
	rm -f $(C_PROGS) $(CXX_PROGS) $(SUM_PROGS) $(DEMOS) *.o

.PHONY: all run clean
