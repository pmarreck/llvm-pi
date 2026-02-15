# llvm-pi Specification

Status: Draft v1.0 (2026-02-15 EST)
Owner: Peter Marreck
Target branch: `yolo`

## 1. Goal

Compute millions of digits of pi as fast as possible using hand-written LLVM IR (`.ll`) for the core algorithm, with C for I/O glue. A secondary goal is to evaluate `codescan` for navigating real-world LLVM IR.

## 2. Done Criteria (Release 0.1)

- `./build` produces a working binary from `.ll` + `.c` sources linked against GMP.
- `./run 1000` prints 1000 correct digits of pi.
- `./run 1000000` completes in reasonable time on Apple Silicon M4.
- `./test` validates correctness at 100, 1000, 10000, and 100000 digits.
- `codescan index` and `codescan search` work on the `.ll` files.
- No external dependencies beyond what `flake.nix` provides.

## 3. Algorithm

Binary splitting Chudnovsky formula.

The Chudnovsky series:
```
1/pi = 12 * sum_{k=0}^{N} (-1)^k * (6k)! * (13591409 + 545140134*k)
                           / ((3k)! * (k!)^3 * 640320^(3k+3/2))
```

Naive term-by-term summation is O(N^2) in digit operations. Binary splitting reduces this to O(N * log(N)^2 * M(N)) where M(N) is the cost of one N-digit multiplication.

Binary splitting computes three integer sequences P(a,b), Q(a,b), T(a,b) over a range [a,b):
- Base case (single term k=a):
  - P(a,a+1) = (6a-5)(2a-1)(6a-1)
  - Q(a,a+1) = a^3 * 640320^3 / 24  (i.e. k^3 * C3_OVER_24)
  - T(a,a+1) = P(a,a+1) * (13591409 + 545140134*a)  [negated if a is odd]
- Merge (split at midpoint m):
  - P(a,b) = P(a,m) * P(m,b)
  - Q(a,b) = Q(a,m) * Q(m,b)
  - T(a,b) = T(a,m) * Q(m,b) + P(a,m) * T(m,b)

Final computation:
  pi = (Q(0,N) * 426880 * sqrt(10005)) / T(0,N)

The sqrt is computed via GMP's mpf_sqrt at the required precision.

## 4. Architecture

```
main.c          CLI: parse digit count, call compute_pi, print result
    |
    v
chudnovsky.ll   Core: binary splitting loop, GMP mpz_* calls
    |
    v
libgmp          External: arbitrary-precision integer/float arithmetic
```

### Files
- `chudnovsky.ll` — hand-written LLVM IR implementing binary splitting Chudnovsky.
- `main.c` — CLI entry point, argument parsing, output formatting.
- `flake.nix` — provides gmp, clang, llvm tools.
- `build` — shell script: llc + clang + link.
- `test` — shell script: correctness checks against known pi digits.

### Build Pipeline
```
chudnovsky.ll --[llc]--> chudnovsky.o
main.c        --[clang]--> main.o
chudnovsky.o + main.o --[clang -lgmp]--> pi
```

## 5. LLVM IR Design

### External GMP Declarations
The `.ll` file declares GMP C functions as externals:
- `@__gmpz_init`, `@__gmpz_clear` — lifecycle
- `@__gmpz_set_ui`, `@__gmpz_set_si` — assignment
- `@__gmpz_mul`, `@__gmpz_mul_ui` — multiplication
- `@__gmpz_add`, `@__gmpz_sub` — addition/subtraction
- `@__gmpz_neg` — negation
- `@__gmpz_tdiv_q` — integer division
- `@__gmpf_init2`, `@__gmpf_clear` — float lifecycle with precision
- `@__gmpf_set_z` — convert integer to float
- `@__gmpf_div`, `@__gmpf_mul` — float arithmetic
- `@__gmpf_sqrt` — square root (for sqrt(10005))
- `@__gmp_sprintf` — convert to decimal string

### Core Functions (in .ll)
- `@split(i64 %a, i64 %b, ptr %P, ptr %Q, ptr %T)` — recursive binary splitting
- `@compute_pi(i64 %digits, ptr %buf, i64 %buf_len) -> i64` — public entry point
- `@base_case(i64 %k, ptr %P, ptr %Q, ptr %T)` — single Chudnovsky term

### Optimization Targets
- **Register pressure**: Minimize stack spills in the split function. Keep P, Q, T pointers in registers across the recursion. Use `noalias` and `nocapture` on pointer params.
- **Cache locality**: Process binary split in depth-first order (natural recursion). The working set at each level is 3 mpz_t values (~limb arrays). For large digit counts, the leaf-level mpz_t values fit in L1; merged values grow into L2/L3.
- **Apple Silicon M4 specifics**:
  - Target `aarch64-apple-darwin` with NEON.
  - Use `llc -mcpu=apple-m4` when available, fallback to `llc -mcpu=apple-m1`.
  - GMP auto-detects Apple Silicon and uses optimized assembly routines for mpn_mul.
  - Align mpz_t allocations to 128 bytes (cache line on M4) via GMP's custom allocator if profiling shows benefit.
- **Fallback**: Generic `aarch64` or `x86_64` targets. The `.ll` uses no target-specific intrinsics directly — all platform optimization comes from `llc`'s code generation and GMP's runtime detection.

## 6. Test Strategy

- Embed or include known pi digits (at least 100K) as a reference.
- `./test` compares output at 100, 1000, 10000, 100000 digits.
- Optional: benchmark mode (`./bm`) timing 100K, 1M, 10M digits using CPU time.

## 7. codescan Evaluation

The `.ll` file will be substantial enough to exercise codescan:
- `codescan symbols chudnovsky.ll` — see function/global hierarchy
- `codescan search "gmpz_mul"` — find all multiplication call sites
- `codescan references --file chudnovsky.ll` — trace SSA value flow
- `codescan find-symbol split` — locate the binary splitting function

## 8. Open Questions

- Use GMP's `mpz_t` struct directly in IR or pass opaque pointers? (GMP's struct layout is stable but platform-specific.)
- Worth implementing a non-recursive iterative binary split to avoid deep call stacks at very high digit counts?
- Add MPFR as optional backend for the final sqrt computation?
