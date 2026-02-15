# llvm-pi Design Document

Date: 2026-02-15
Status: Approved

## Summary

Compute millions of digits of pi using hand-written LLVM IR for the core binary splitting Chudnovsky algorithm, with C for CLI/I/O, linked against GMP for arbitrary-precision arithmetic. Secondary goal: evaluate `codescan` for navigating real LLVM IR.

## Approach: Binary Splitting Chudnovsky with GMP mpz_t

All core arithmetic in hand-written `.ll` calling GMP's `mpz_*` (integer) functions. One final `mpf_*` (float) division + sqrt to produce the decimal result. C `main()` handles argument parsing and output.

### Why this approach
- `mpz_t` integer accumulation avoids precision loss from repeated float operations
- Binary splitting is O(N log(N)^2 M(N)) vs O(N^2) for naive summation
- GMP's mpn layer auto-detects Apple Silicon NEON routines at runtime
- Hand-written `.ll` produces substantial IR for codescan evaluation

### Key design decisions
1. **Recursive binary split** (not iterative) — natural depth-first traversal gives good cache behavior; stack depth is O(log N) which is ~20 for 1M digits, ~24 for 10M.
2. **GMP as opaque pointers** — pass `ptr` to mpz_t structs, don't model GMP internals in IR. Safer and portable.
3. **Platform optimization via llc flags** — the `.ll` itself is target-independent. `llc -mcpu=apple-m4` (or fallback) handles instruction selection. `noalias`/`nocapture` hints help register allocation.
4. **Cache alignment** — GMP's default allocator is already well-aligned. Custom allocator only if profiling shows benefit.

## File layout

```
llvm-pi/
  chudnovsky.ll     Core algorithm in LLVM IR
  main.c            CLI entry, I/O
  build             Build script (bash)
  test              Test script (bash)
  flake.nix         Nix dependencies
  flake.lock        Nix lock file
  pi_reference.txt  Known pi digits for testing
  SPEC.md           Project specification
  AGENTS.md         Agent instructions (symlink)
  docs/plans/       Design documents
```

## Build pipeline

```bash
# Assemble .ll to .o
nix develop -c llc -O2 -mcpu=apple-m4 -filetype=obj chudnovsky.ll -o chudnovsky.o

# Compile C
nix develop -c clang -O2 -c main.c -o main.o

# Link
nix develop -c clang -O2 main.o chudnovsky.o -lgmp -o pi
```

Fallback for non-M4:
```bash
nix develop -c llc -O2 -filetype=obj chudnovsky.ll -o chudnovsky.o
```

## Optimization notes

### Register / instruction level
- `noalias` and `nocapture` on all mpz_t pointer params in split()
- `nounwind` on all GMP extern declarations (they don't throw)
- `willreturn` where applicable
- Minimize temporaries — reuse mpz_t variables across iterations where safe

### Cache level
- Depth-first recursion naturally processes small (L1-sized) leaf terms first
- At merge levels, operands grow but are touched sequentially (good prefetch behavior)
- GMP internally manages limb array growth; no explicit cache management needed in IR

### Apple Silicon M4
- Target triple: `aarch64-apple-darwin`
- llc flag: `-mcpu=apple-m4` (falls back gracefully if not recognized)
- GMP compiles with `--enable-assembly` in Nix, picks up NEON mpn routines
- 128-byte cache lines on M4 — GMP's 16-byte aligned mallocs are sufficient for limb arrays

### Fallback
- No target-specific intrinsics in the `.ll` — it's portable IR
- `llc` without `-mcpu` generates generic aarch64 or x86_64
- GMP's runtime CPU detection handles the rest

## Testing

Reference digits from a trusted source (bundled file or hardcoded). Tests at:
- 100 digits (sanity)
- 1,000 digits (basic correctness)
- 10,000 digits (intermediate)
- 100,000 digits (stress)

Optional benchmark at 1M and 10M digits.
