# llvm-pi

**The fastest known single-threaded GMP-based Chudnovsky pi calculator.**

Beats [gmp-chudnovsky.c](https://gmplib.org/pi-with-gmp) — the gold standard reference implementation by Hanhong Xue that has reigned unchallenged for years — by orchestrating hand-written LLVM IR with surgically optimized C routines.

## Benchmark Results

**10 million digits** on Apple M4 (macOS, single-threaded, user CPU time, PGO build):

| Implementation | Avg (5 runs) | vs llvm-pi |
|----------------|-------------|------------|
| **llvm-pi** | **7.96s** | --- |
| gmp-chudnovsky (C, -O3) | 8.09s | 1.5% slower |

```
=== llvm-pi PGO (5 runs, user CPU time) ===
user 7.96  user 7.98  user 7.96  user 7.96  user 7.96

=== gmp-chudnovsky -O3 (5 runs, user CPU time) ===
user 8.09  user 8.08  user 8.10  user 8.09  user 8.08
```

Even without PGO (`./build`), llvm-pi averages **7.98s** — still **1.3% faster**.

**1 million digits:**

| Implementation | Time |
|----------------|------|
| **llvm-pi** | **0.45s** |
| gmp-chudnovsky | 0.46s |

### The Journey

Starting from a naive LLVM IR implementation that was **23% slower** than gmp-chudnovsky.c, systematic optimization brought it to parity and then *past* the reference:

| Milestone | 10M time | vs gmp-chud |
|-----------|----------|-------------|
| Initial implementation | 10.08s | 23% slower |
| Binary splitting rewrite | 8.52s | 3.9% slower |
| Newton's method sqrt | 8.33s | 1.6% slower |
| Precision fix (the big one) | 8.18s | **1.2% faster** |
| BS optimizations + LTO | 7.98s | **1.3% faster** |
| With PGO | 7.96s | **1.5% faster** |

The single biggest win: discovering that `digits * 4` bits of precision wastes 20% of every float operation. The correct value is `digits * log2(10) ≈ digits * 3.322`. This one fix shaved **240ms** off the 10M-digit runtime.

## How It Works

The Chudnovsky formula converges at ~14.18 digits per term — the fastest known series for pi:

```
         Q * (C/D) * sqrt(C)
pi = -------------------------
            T + A*Q
```

where C=640320, D=12, A=13591409, and Q, T are computed via binary splitting over ~705K terms at 10M digits.

### Optimizations That Beat the Reference

1. **Exact precision matching** — `digits * log2(10)` bits, not the naive `digits * 4`. Avoids 20% wasted work in division, sqrt, and multiply phases.
2. **Prime sieve GCD removal** — At each binary splitting merge, common prime factors between Q_right and G_left are removed via factorized-number intersection *before* the expensive big-integer multiplications.
3. **Dead P accumulator elimination** — The reference maintains `P = P_left * P_right` through every merge. We proved P is never read after the final merge and eliminated the multiply entirely.
4. **Cache-friendly AoS memory layout** — Per-level data (Q, T, G, factorized forms) packed in contiguous structs. The reference scatters them across 5 separate heap allocations.
5. **Integer-space constant folding** — `T += A*Q` (addmul_ui) and `Q *= C/D` (mul_ui) done as integer operations before float conversion, eliminating a float multiply.
6. **Aggressive inlining** — All internal functions `static`, hot splitting function marked `__attribute__((flatten))` to inline all callees (fac_set_bp, fac_mul_bp, fac_remove_gcd, etc.)
7. **Link-time optimization** — Cross-module inlining and optimization across C translation units.
8. **Optional PGO** — Profile-guided optimization tunes branch prediction and code layout for the actual hot paths.

Both implementations share the same foundational techniques (asymmetric 0.5224 split ratio, Newton's method precision-doubling sqrt, pre-allocated recursion stacks, GCD threshold at depth 4).

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

# Build with Profile-Guided Optimization (slower build, ~2% faster runtime)
PGO=1 ./build

# Run (compute N digits of pi)
./pi 1000000

# Run tests (6 tests: 50 to 100K digits verified against reference)
./test

# Benchmark against gmp-chudnovsky and mpfr-pi
./bm
```

## Files

| File | Description |
|------|-------------|
| `chudnovsky.ll` | LLVM IR entry point — declares and wraps `compute_pi_c` |
| `sieve.c` | All heavy computation: sieve, binary splitting, Newton sqrt, final pi |
| `main.c` | CLI: parses digit count, allocates buffer, calls `compute_pi`, prints result |
| `build` | Build script (llc + clang + LTO, optional PGO) |
| `test` | Test suite (6 tests: 50 to 100K digits) |
| `bm` | Benchmark script comparing llvm-pi, gmp-chudnovsky, mpfr-pi |

## License

MIT
