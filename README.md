# llvm-pi

**The fastest known single-threaded GMP-based Chudnovsky pi calculator.**

Hand-written LLVM IR orchestrating optimized C routines to compute arbitrary-precision digits of pi using the binary splitting Chudnovsky algorithm with GMP.

## Benchmark Results

**10 million digits** on Apple M4 (macOS, single-threaded, CPU time):

| Implementation | Avg (5 runs) | vs llvm-pi |
|----------------|-------------|------------|
| **llvm-pi** | **7.97s** | baseline |
| gmp-chudnovsky (C, -O3) | 8.09s | +1.5% slower |

```
=== llvm-pi (5 runs) ===
user 7.97  user 7.98  user 7.97  user 7.98  user 7.96

=== gmp-chudnovsky -O3 (5 runs) ===
user 8.09  user 8.09  user 8.10  user 8.07  user 8.10
```

The reference implementation is [gmp-chudnovsky.c](https://gmplib.org/pi-with-gmp) by Hanhong Xue, which has been the standard high-performance GMP-based pi calculator for years.

**1 million digits:**

| Implementation | Time |
|----------------|------|
| **llvm-pi** | **0.45s** |
| gmp-chudnovsky | 0.46s |

## How It Works

The Chudnovsky formula converges at ~14.18 digits per term:

```
         Q * (C/D) * sqrt(C)
pi = -------------------------
               T
```

The implementation uses **binary splitting** to recursively divide the series into halves, merging with GMP big-integer arithmetic. Key optimizations that beat the reference C implementation:

1. **Prime sieve GCD removal** — At each merge step, common factors between Q_right and G_left are removed before multiplication, keeping intermediate numbers smaller
2. **Asymmetric split point** (0.5224 ratio) — Tuned from gmp-chudnovsky, balances work between left/right children
3. **Dead P accumulator elimination** — The P = P_left * P_right multiply is skipped since P is never read after the final merge
4. **Pre-allocated recursion stacks** — Arrays of mpz_t indexed by depth eliminate per-level malloc/free
5. **Newton's method sqrt** — Precision-doubling iteration for 1/sqrt(x), much faster than mpf_sqrt at high precision
6. **Exact precision matching** — Uses `digits * log2(10) + 256` bits instead of the naive `digits * 4`, avoiding 20% wasted computation in sqrt/multiply/division phases
7. **Integer-space constant folding** — Q *= C/D done as an integer multiply before float conversion, eliminating a float multiply
8. **Link-time optimization** — Cross-module inlining across C translation units

## Architecture

```
main.c (CLI)
  -> chudnovsky.ll (LLVM IR entry point)
    -> sieve.c (all heavy computation)
       - Prime sieve + factorized number GCD
       - Binary splitting with pre-allocated stacks
       - Newton's method square root
       - Final pi computation + formatting
```

The LLVM IR serves as the public interface (`compute_pi`), delegating to optimized C routines for the number crunching.

## Building

Requires [Nix](https://nixos.org/) with flakes enabled:

```bash
# Build
./build

# Run (compute N digits of pi)
./pi 1000000

# Run tests
./test

# Benchmark against reference implementations
./bm
```

## Files

| File | Description |
|------|-------------|
| `chudnovsky.ll` | LLVM IR entry point — declares and wraps `compute_pi_c` |
| `sieve.c` | All heavy computation: sieve, binary splitting, Newton sqrt, final pi |
| `main.c` | CLI: parses digit count, allocates buffer, calls `compute_pi`, prints result |
| `build` | Build script (llc + clang + LTO) |
| `test` | Test suite (6 tests: 50 to 100K digits) |
| `bm` | Benchmark script comparing llvm-pi, gmp-chudnovsky, mpfr-pi |

## License

MIT
