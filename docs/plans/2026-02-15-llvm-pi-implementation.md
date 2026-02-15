# llvm-pi Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Compute millions of digits of pi using hand-written LLVM IR (binary splitting Chudnovsky), C for I/O, GMP for bignum.

**Architecture:** `main.c` (CLI) calls `compute_pi()` defined in `chudnovsky.ll`, which calls GMP `mpz_*`/`mpf_*` externals. Build via `llc` + `clang` + link.

**Tech Stack:** LLVM IR (.ll), C, GMP 6.3, Nix flake, Bash scripts.

---

### Task 1: Create flake.nix

**Files:**
- Create: `flake.nix`

**Step 1: Write flake.nix**

```nix
{
  description = "llvm-pi: compute pi in hand-written LLVM IR";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      supportedSystems = [ "aarch64-darwin" "x86_64-darwin" "aarch64-linux" "x86_64-linux" ];
      forAllSystems = nixpkgs.lib.genAttrs supportedSystems;
    in {
      devShells = forAllSystems (system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
        in {
          default = pkgs.mkShell {
            buildInputs = [
              pkgs.gmp
              pkgs.llvmPackages_19.llvm
              pkgs.clang_19
            ];
            shellHook = ''
              export GMP_INCLUDE="${pkgs.gmp.dev}/include"
              export GMP_LIB="${pkgs.gmp}/lib"
            '';
          };
        });
    };
}
```

**Step 2: Verify flake works**

Run: `nix develop -c bash -c 'llc --version | head -2 && clang --version | head -1 && ls $GMP_LIB/libgmp.*'`
Expected: LLVM 19, clang 19, libgmp files listed.

**Step 3: Commit**

```bash
git add flake.nix
git commit -m "Add flake.nix with GMP, LLVM 19, clang 19"
```

---

### Task 2: Create build script

**Files:**
- Create: `build`

**Step 1: Write build script**

```bash
#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

# Detect CPU for llc
detect_mcpu() {
	if [[ "$(uname -m)" == "arm64" || "$(uname -m)" == "aarch64" ]]; then
		# Try apple-m4, fall back to apple-m1, fall back to generic
		for cpu in apple-m4 apple-m1 generic; do
			if nix develop -c llc -mcpu="$cpu" -mtriple=aarch64-apple-darwin /dev/null -o /dev/null 2>/dev/null; then
				echo "$cpu"
				return
			fi
		done
		echo "generic"
	else
		echo "generic"
	fi
}

OPT="${1:--O2}"
MCPU="$(detect_mcpu)"
TRIPLE=""
if [[ "$(uname -m)" == "arm64" || "$(uname -m)" == "aarch64" ]]; then
	TRIPLE="-mtriple=aarch64-apple-darwin"
fi

echo "Building with -mcpu=$MCPU $OPT"

# Assemble LLVM IR to object file
nix develop -c llc "$OPT" -mcpu="$MCPU" $TRIPLE -filetype=obj chudnovsky.ll -o chudnovsky.o

# Compile C
nix develop -c clang "$OPT" -I"$GMP_INCLUDE" -c main.c -o main.o

# Link
nix develop -c clang "$OPT" main.o chudnovsky.o -L"$GMP_LIB" -lgmp -o pi

echo "Built: ./pi"
```

**Step 2: Make executable**

Run: `chmod +x build`

**Step 3: Commit**

```bash
git add build
git commit -m "Add build script with Apple Silicon CPU detection"
```

---

### Task 3: Create pi_reference.txt

**Files:**
- Create: `pi_reference.txt`

**Step 1: Generate reference digits**

Use `bc -l` to generate 1000+ digits of pi, then augment to 100K+ from a known source. For now, start with 1000 digits (the test script will use prefix matching):

```
3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989
```

Save this as `pi_reference.txt` (the "3." prefix followed by digits, no newlines mid-number).

**Step 2: Commit**

```bash
git add pi_reference.txt
git commit -m "Add pi reference digits (1000+) for testing"
```

---

### Task 4: Write chudnovsky.ll — GMP external declarations

**Files:**
- Create: `chudnovsky.ll`

**Step 1: Write the external declarations and struct types**

This is the foundation — just the declarations, no algorithm yet. We'll add a trivial `compute_pi` stub that sets pi to "3.14" to validate the build pipeline.

```llvm
; chudnovsky.ll — Binary splitting Chudnovsky pi computation
; Core algorithm in LLVM IR, calls GMP for arbitrary-precision arithmetic.

; === GMP struct types ===
; __mpz_struct = { int _mp_alloc, int _mp_size, mp_limb_t* _mp_d }
%struct.mpz = type { i32, i32, ptr }

; __mpf_struct = { int _mp_prec, int _mp_size, mp_exp_t _mp_exp, mp_limb_t* _mp_d }
%struct.mpf = type { i32, i32, i64, ptr }

; === GMP integer function declarations ===
declare void @__gmpz_init(ptr) nounwind
declare void @__gmpz_clear(ptr) nounwind
declare void @__gmpz_set_ui(ptr, i64) nounwind
declare void @__gmpz_set_si(ptr, i64) nounwind
declare void @__gmpz_mul(ptr, ptr, ptr) nounwind
declare void @__gmpz_mul_ui(ptr, ptr, i64) nounwind
declare void @__gmpz_mul_si(ptr, ptr, i64) nounwind
declare void @__gmpz_add(ptr, ptr, ptr) nounwind
declare void @__gmpz_sub(ptr, ptr, ptr) nounwind
declare void @__gmpz_addmul(ptr, ptr, ptr) nounwind
declare void @__gmpz_submul(ptr, ptr, ptr) nounwind
declare void @__gmpz_neg(ptr, ptr) nounwind
declare void @__gmpz_tdiv_q(ptr, ptr, ptr) nounwind

; === GMP float function declarations ===
declare void @__gmpf_init2(ptr, i64) nounwind
declare void @__gmpf_clear(ptr) nounwind
declare void @__gmpf_set_z(ptr, ptr) nounwind
declare void @__gmpf_div(ptr, ptr, ptr) nounwind
declare void @__gmpf_mul(ptr, ptr, ptr) nounwind
declare void @__gmpf_mul_ui(ptr, ptr, i64) nounwind
declare void @__gmpf_sqrt(ptr, ptr) nounwind

; === GMP output ===
declare i32 @__gmp_sprintf(ptr, ptr, ...) nounwind

; === C stdlib ===
declare ptr @malloc(i64) nounwind
declare void @free(ptr) nounwind

; === String constants ===
@.fmt_pi = private unnamed_addr constant [6 x i8] c"%Ff\0A\00", align 1
@.fmt_fixed = private unnamed_addr constant [8 x i8] c"%.*Ff\0A\00", align 1

; === Public entry point ===
; Computes pi to the requested number of decimal digits.
; Writes the decimal string to buf (must be at least digits+10 bytes).
; Returns the number of characters written.
define i64 @compute_pi(i64 %digits, ptr noalias nocapture %buf, i64 %buf_len) nounwind {
entry:
  ; --- STUB: just write "3.14" for build pipeline validation ---
  ; This will be replaced with the real algorithm in Task 6.
  store i8 51, ptr %buf               ; '3'
  %p1 = getelementptr i8, ptr %buf, i64 1
  store i8 46, ptr %p1                ; '.'
  %p2 = getelementptr i8, ptr %buf, i64 2
  store i8 49, ptr %p2                ; '1'
  %p3 = getelementptr i8, ptr %buf, i64 3
  store i8 52, ptr %p3                ; '4'
  %p4 = getelementptr i8, ptr %buf, i64 4
  store i8 0, ptr %p4                 ; null terminator
  ret i64 4
}
```

**Step 2: Verify it assembles**

Run: `nix develop -c llc -O2 -filetype=obj chudnovsky.ll -o chudnovsky.o && echo "OK"`
Expected: `OK`

**Step 3: Commit**

```bash
git add chudnovsky.ll
git commit -m "Add chudnovsky.ll stub with GMP declarations"
```

---

### Task 5: Write main.c and validate end-to-end build

**Files:**
- Create: `main.c`

**Step 1: Write main.c**

```c
/* main.c — CLI entry for llvm-pi */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Defined in chudnovsky.ll */
extern long long compute_pi(long long digits, char *buf, long long buf_len);

int main(int argc, char *argv[]) {
	long long digits = 100; /* default */

	if (argc > 1) {
		digits = atoll(argv[1]);
		if (digits < 1) {
			fprintf(stderr, "Usage: %s [digits]\n", argv[0]);
			return 1;
		}
	}

	/* Buffer: "3." + digits + null + some margin */
	long long buf_len = digits + 64;
	char *buf = malloc((size_t)buf_len);
	if (!buf) {
		fprintf(stderr, "Failed to allocate %lld bytes\n", buf_len);
		return 1;
	}

	long long written = compute_pi(digits, buf, buf_len);
	if (written > 0) {
		puts(buf);
	} else {
		fprintf(stderr, "compute_pi failed\n");
		free(buf);
		return 1;
	}

	free(buf);
	return 0;
}
```

**Step 2: Build end-to-end**

Run: `nix develop -c bash -c './build'`
Expected: `Built: ./pi`

**Step 3: Run the stub**

Run: `nix develop -c ./pi 100`
Expected: `3.14` (stub output)

**Step 4: Commit**

```bash
git add main.c
git commit -m "Add main.c CLI entry point"
```

---

### Task 6: Write the test script (against stub first)

**Files:**
- Create: `test`

**Step 1: Write test script**

```bash
#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

PASS=0
FAIL=0
PI_REF="$(cat pi_reference.txt)"

run_test() {
	local digits="$1"
	local label="$2"
	# Expected: "3." + first $digits decimal digits from reference
	local expected
	expected="$(echo "$PI_REF" | head -c $((digits + 2)))"  # "3." = 2 chars

	local actual
	actual="$(nix develop -c ./pi "$digits" 2>&1)" || {
		echo "FAIL [$label]: pi $digits crashed"
		FAIL=$((FAIL + 1))
		return
	}

	# Compare: actual must start with expected
	if [[ "$actual" == "$expected"* ]] || [[ "$expected" == "$actual"* && ${#actual} -ge ${#expected} ]]; then
		echo "PASS [$label]: $digits digits match"
		PASS=$((PASS + 1))
	else
		echo "FAIL [$label]: $digits digits mismatch"
		echo "  expected prefix: ${expected:0:40}..."
		echo "  actual prefix:   ${actual:0:40}..."
		FAIL=$((FAIL + 1))
	fi
}

# Build first
echo "=== Building ==="
nix develop -c bash -c './build' || { echo "Build failed"; exit 1; }
echo

echo "=== Running tests ==="
run_test 50 "sanity-50"
run_test 100 "basic-100"
run_test 1000 "medium-1000"

echo
echo "=== Results: $PASS passed, $FAIL failed ==="
exit "$FAIL"
```

**Step 2: Make executable**

Run: `chmod +x test`

**Step 3: Run tests (expect failures — stub only returns "3.14")**

Run: `nix develop -c bash -c './test'`
Expected: All tests FAIL (stub output is "3.14", not full digits). This confirms the test harness works and catches wrong output.

**Step 4: Commit**

```bash
git add test
git commit -m "Add test script with digit-comparison harness"
```

---

### Task 7: Implement base_case in chudnovsky.ll

**Files:**
- Modify: `chudnovsky.ll`

**Step 1: Add the base_case function**

This computes a single Chudnovsky term for index k:
- P(k) = (6k-5)(2k-1)(6k-1)
- Q(k) = k^3 * C3_OVER_24  where C3_OVER_24 = 640320^3 / 24 = 10939058860032000
- T(k) = P(k) * (13591409 + 545140134*k), negated if k is odd

Special case: k=0 → P=1, Q=1, T=13591409

Add this function after the declarations, before `compute_pi`:

```llvm
; Chudnovsky constants
; C3_OVER_24 = 640320^3 / 24 = 10939058860032000
; A = 13591409
; B = 545140134

; base_case: compute P, Q, T for a single term k
; For k=0: P=1, Q=1, T=A=13591409
; For k>0: P=(6k-5)(2k-1)(6k-1), Q=k^3*C3_OVER_24, T=P*(A+B*k)*(-1)^k
define void @base_case(i64 %k, ptr noalias nocapture %P, ptr noalias nocapture %Q, ptr noalias nocapture %T) nounwind {
entry:
  %is_zero = icmp eq i64 %k, 0
  br i1 %is_zero, label %zero_case, label %nonzero_case

zero_case:
  ; P = 1, Q = 1, T = 13591409
  call void @__gmpz_set_ui(ptr %P, i64 1)
  call void @__gmpz_set_ui(ptr %Q, i64 1)
  call void @__gmpz_set_si(ptr %T, i64 13591409)
  ret void

nonzero_case:
  ; P = (6k-5) * (2k-1) * (6k-1)
  %k6 = mul i64 %k, 6
  %a1 = sub i64 %k6, 5          ; 6k - 5
  %k2 = mul i64 %k, 2
  %a2 = sub i64 %k2, 1          ; 2k - 1
  %a3 = sub i64 %k6, 1          ; 6k - 1

  call void @__gmpz_set_si(ptr %P, i64 %a1)
  ; tmp = P (reuse T temporarily for multiplication)
  call void @__gmpz_mul_si(ptr %P, ptr %P, i64 %a2)
  call void @__gmpz_mul_si(ptr %P, ptr %P, i64 %a3)

  ; Q = k^3 * C3_OVER_24
  ; C3_OVER_24 = 10939058860032000
  call void @__gmpz_set_si(ptr %Q, i64 %k)
  call void @__gmpz_mul(ptr %Q, ptr %Q, ptr %Q)          ; k^2
  call void @__gmpz_mul_si(ptr %Q, ptr %Q, i64 %k)       ; k^3
  call void @__gmpz_mul_ui(ptr %Q, ptr %Q, i64 10939058860032000) ; k^3 * C3_OVER_24

  ; T = P * (A + B*k)
  ; A + B*k = 13591409 + 545140134*k
  %bk = mul i64 545140134, %k
  %abk = add i64 13591409, %bk
  call void @__gmpz_set_si(ptr %T, i64 %abk)
  call void @__gmpz_mul(ptr %T, ptr %T, ptr %P)

  ; Negate if k is odd
  %k_odd = and i64 %k, 1
  %is_odd = icmp ne i64 %k_odd, 0
  br i1 %is_odd, label %negate, label %done

negate:
  call void @__gmpz_neg(ptr %T, ptr %T)
  br label %done

done:
  ret void
}
```

**Step 2: Verify it assembles**

Run: `nix develop -c llc -O2 -filetype=obj chudnovsky.ll -o chudnovsky.o && echo "OK"`
Expected: `OK`

**Step 3: Commit**

```bash
git add chudnovsky.ll
git commit -m "Add base_case: single Chudnovsky term in LLVM IR"
```

---

### Task 8: Implement split (binary splitting) in chudnovsky.ll

**Files:**
- Modify: `chudnovsky.ll`

**Step 1: Add the split function**

Recursive binary splitting over range [a, b). When b-a == 1, calls base_case. Otherwise splits at midpoint and merges.

```llvm
; split: recursive binary splitting over range [a, b)
; Computes P(a,b), Q(a,b), T(a,b) into caller-provided mpz_t pointers.
; Uses 6 temporary mpz_t values for the right-half results and merge.
define void @split(i64 %a, i64 %b, ptr noalias nocapture %P, ptr noalias nocapture %Q, ptr noalias nocapture %T) nounwind {
entry:
  %diff = sub i64 %b, %a
  %is_base = icmp eq i64 %diff, 1
  br i1 %is_base, label %do_base, label %do_split

do_base:
  call void @base_case(i64 %a, ptr %P, ptr %Q, ptr %T)
  ret void

do_split:
  ; midpoint m = (a + b) / 2
  %sum = add i64 %a, %b
  %m = sdiv i64 %sum, 2

  ; Recurse left: split(a, m, P, Q, T)  — results go into P, Q, T directly
  call void @split(i64 %a, i64 %m, ptr %P, ptr %Q, ptr %T)

  ; Allocate temporaries for right half: Pr, Qr, Tr
  %Pr = alloca %struct.mpz, align 16
  %Qr = alloca %struct.mpz, align 16
  %Tr = alloca %struct.mpz, align 16
  call void @__gmpz_init(ptr %Pr)
  call void @__gmpz_init(ptr %Qr)
  call void @__gmpz_init(ptr %Tr)

  ; Recurse right: split(m, b, Pr, Qr, Tr)
  call void @split(i64 %m, i64 %b, ptr %Pr, ptr %Qr, ptr %Tr)

  ; Merge:
  ; T(a,b) = T(a,m) * Qr + P(a,m) * Tr
  ;   Step 1: T = T * Qr
  call void @__gmpz_mul(ptr %T, ptr %T, ptr %Qr)
  ;   Step 2: temp = P * Tr  (reuse Tr as temp since we're done with it after)
  %temp = alloca %struct.mpz, align 16
  call void @__gmpz_init(ptr %temp)
  call void @__gmpz_mul(ptr %temp, ptr %P, ptr %Tr)
  ;   Step 3: T = T + temp
  call void @__gmpz_add(ptr %T, ptr %T, ptr %temp)
  call void @__gmpz_clear(ptr %temp)

  ; P(a,b) = P(a,m) * Pr
  call void @__gmpz_mul(ptr %P, ptr %P, ptr %Pr)

  ; Q(a,b) = Q(a,m) * Qr
  call void @__gmpz_mul(ptr %Q, ptr %Q, ptr %Qr)

  ; Cleanup right-half temporaries
  call void @__gmpz_clear(ptr %Pr)
  call void @__gmpz_clear(ptr %Qr)
  call void @__gmpz_clear(ptr %Tr)

  ret void
}
```

**Step 2: Verify it assembles**

Run: `nix develop -c llc -O2 -filetype=obj chudnovsky.ll -o chudnovsky.o && echo "OK"`
Expected: `OK`

**Step 3: Commit**

```bash
git add chudnovsky.ll
git commit -m "Add split: recursive binary splitting in LLVM IR"
```

---

### Task 9: Implement compute_pi (replace stub)

**Files:**
- Modify: `chudnovsky.ll`

**Step 1: Replace the compute_pi stub with the real implementation**

This:
1. Calculates N (number of terms needed for requested digits)
2. Calls split(0, N, P, Q, T)
3. Computes pi = Q * 426880 * sqrt(10005) / T using mpf_t
4. Formats result to decimal string

```llvm
; compute_pi: compute pi to `digits` decimal places.
; Writes decimal string to buf, returns chars written.
define i64 @compute_pi(i64 %digits, ptr noalias nocapture %buf, i64 %buf_len) nounwind {
entry:
  ; Number of Chudnovsky terms needed: ~1 term per 14.18 digits
  ; N = digits / 14 + 2  (with safety margin)
  %d14 = sdiv i64 %digits, 14
  %N = add i64 %d14, 2

  ; Precision in bits: digits * log2(10) + 64 safety bits
  ; log2(10) ~ 3.3219, approximate as digits * 4 (slight overestimate is fine)
  %prec_approx = mul i64 %digits, 4
  %prec = add i64 %prec_approx, 128

  ; Initialize P, Q, T
  %P = alloca %struct.mpz, align 16
  %Q = alloca %struct.mpz, align 16
  %T = alloca %struct.mpz, align 16
  call void @__gmpz_init(ptr %P)
  call void @__gmpz_init(ptr %Q)
  call void @__gmpz_init(ptr %T)

  ; Binary splitting: split(0, N, P, Q, T)
  call void @split(i64 0, i64 %N, ptr %P, ptr %Q, ptr %T)

  ; Now compute pi = Q * 426880 * sqrt(10005) / T
  ; All in mpf_t at required precision.

  ; Allocate mpf_t variables
  %pi_f = alloca %struct.mpf, align 16
  %q_f = alloca %struct.mpf, align 16
  %t_f = alloca %struct.mpf, align 16
  %sqrt_f = alloca %struct.mpf, align 16
  call void @__gmpf_init2(ptr %pi_f, i64 %prec)
  call void @__gmpf_init2(ptr %q_f, i64 %prec)
  call void @__gmpf_init2(ptr %t_f, i64 %prec)
  call void @__gmpf_init2(ptr %sqrt_f, i64 %prec)

  ; sqrt_f = sqrt(10005.0)
  ; First set an mpz to 10005 then convert to mpf, then sqrt
  %ten005 = alloca %struct.mpz, align 16
  call void @__gmpz_init(ptr %ten005)
  call void @__gmpz_set_ui(ptr %ten005, i64 10005)
  call void @__gmpf_set_z(ptr %sqrt_f, ptr %ten005)
  call void @__gmpf_sqrt(ptr %sqrt_f, ptr %sqrt_f)
  call void @__gmpz_clear(ptr %ten005)

  ; q_f = Q (convert integer to float)
  call void @__gmpf_set_z(ptr %q_f, ptr %Q)

  ; t_f = T (convert integer to float)
  call void @__gmpf_set_z(ptr %t_f, ptr %T)

  ; q_f = q_f * 426880
  call void @__gmpf_mul_ui(ptr %q_f, ptr %q_f, i64 426880)

  ; q_f = q_f * sqrt_f   (Q * 426880 * sqrt(10005))
  call void @__gmpf_mul(ptr %q_f, ptr %q_f, ptr %sqrt_f)

  ; pi_f = q_f / t_f
  call void @__gmpf_div(ptr %pi_f, ptr %q_f, ptr %t_f)

  ; Format to string: "%.<digits>Ff"
  ; Build format string on stack: "%." + digits_str + "Ff\0"
  ; Simpler approach: use __gmp_sprintf with "%.*Ff" and pass digits as int arg
  ; But gmp_sprintf %.*Ff takes an int, not i64. We'll use a fixed format.
  ; Actually, use gmp_sprintf(buf, "%Ff", pi_f) — prints all significant digits.
  ; Then truncate in C. Simpler and avoids format string construction in IR.
  %written = call i32 (ptr, ptr, ...) @__gmp_sprintf(ptr %buf, ptr @.fmt_pi, ptr %pi_f)

  ; Cleanup
  call void @__gmpf_clear(ptr %pi_f)
  call void @__gmpf_clear(ptr %q_f)
  call void @__gmpf_clear(ptr %t_f)
  call void @__gmpf_clear(ptr %sqrt_f)
  call void @__gmpz_clear(ptr %P)
  call void @__gmpz_clear(ptr %Q)
  call void @__gmpz_clear(ptr %T)

  %written64 = sext i32 %written to i64
  ret i64 %written64
}
```

**Step 2: Build and do a smoke test**

Run: `nix develop -c bash -c './build && ./pi 50'`
Expected: Pi digits starting with 3.14159... (may have precision issues — we'll fix formatting next)

**Step 3: Commit**

```bash
git add chudnovsky.ll
git commit -m "Implement compute_pi: full Chudnovsky binary splitting"
```

---

### Task 10: Fix output formatting in main.c

**Files:**
- Modify: `main.c`

**Step 1: Update main.c to truncate output to requested digit count**

The `%Ff` format may output more or fewer digits than requested. main.c should truncate to exactly the requested number of decimal digits after "3.":

```c
/* main.c — CLI entry for llvm-pi */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Defined in chudnovsky.ll */
extern long long compute_pi(long long digits, char *buf, long long buf_len);

int main(int argc, char *argv[]) {
	long long digits = 100; /* default */

	if (argc > 1) {
		digits = atoll(argv[1]);
		if (digits < 1) {
			fprintf(stderr, "Usage: %s [digits]\n", argv[0]);
			return 1;
		}
	}

	/* Buffer needs room for all GMP output digits + overhead */
	long long buf_len = digits + 256;
	char *buf = malloc((size_t)buf_len);
	if (!buf) {
		fprintf(stderr, "Failed to allocate %lld bytes\n", buf_len);
		return 1;
	}

	long long written = compute_pi(digits, buf, buf_len);
	if (written <= 0) {
		fprintf(stderr, "compute_pi failed\n");
		free(buf);
		return 1;
	}

	/* Find the decimal point and truncate to requested digits */
	char *dot = strchr(buf, '.');
	if (dot && (long long)(strlen(dot + 1)) > digits) {
		/* Truncate: keep "3." + exactly `digits` decimal chars */
		dot[1 + digits] = '\n';
		dot[2 + digits] = '\0';
	}

	/* Remove trailing newline from gmp_sprintf if present before our newline */
	fputs(buf, stdout);
	if (buf[strlen(buf) - 1] != '\n')
		putchar('\n');

	free(buf);
	return 0;
}
```

**Step 2: Build and test**

Run: `nix develop -c bash -c './build && ./pi 50'`
Expected: Exactly 50 decimal digits after "3."

**Step 3: Run the test suite**

Run: `nix develop -c bash -c './test'`
Expected: Tests should start passing (at least the 50-digit sanity test).

**Step 4: Commit**

```bash
git add main.c
git commit -m "Fix output truncation to exact digit count"
```

---

### Task 11: Debug and fix any correctness issues

**Files:**
- Possibly modify: `chudnovsky.ll`, `main.c`

**Step 1: Run tests and diagnose failures**

Run: `nix develop -c bash -c './test'`

Common issues to check:
- **Off-by-one in digit count**: the `%Ff` format may include or exclude the leading `3.`
- **Precision too low**: if digits are wrong partway through, increase the `prec` calculation (try `digits * 5`)
- **Sign issue in T**: the Chudnovsky series alternates sign; verify k=0 gives positive T
- **N too small**: if last digits are wrong, increase the `+2` safety margin in N calculation

**Step 2: Fix whatever broke**

Apply minimal fixes based on test output.

**Step 3: Run tests again**

Run: `nix develop -c bash -c './test'`
Expected: All tests pass.

**Step 4: Commit**

```bash
git add chudnovsky.ll main.c
git commit -m "Fix correctness issues found by test suite"
```

---

### Task 12: Extend pi_reference.txt to 100K digits and add stress test

**Files:**
- Modify: `pi_reference.txt`, `test`

**Step 1: Generate 100K+ reference digits**

Use the now-working `./pi` itself to generate reference digits at high precision, then verify the first 1000 against our known-good reference. Alternatively, use an external source.

Actually, the most trustworthy approach: use our program at *very high* precision (e.g., 110000 digits) and compare the first 1000 against our existing reference. If those match, trust the rest.

Run: `nix develop -c bash -c './build && ./pi 100100 > pi_100k.txt && head -c 1002 pi_100k.txt'`
Verify first 1000 digits match `pi_reference.txt`.

If they match, replace `pi_reference.txt` with the 100K version:
```bash
mv pi_100k.txt pi_reference.txt
```

**Step 2: Add 10K and 100K tests to test script**

Add these lines to the test script:
```bash
run_test 10000 "stress-10k"
run_test 100000 "stress-100k"
```

**Step 3: Run full test suite**

Run: `nix develop -c bash -c './test'`
Expected: All tests pass including 100K.

**Step 4: Commit**

```bash
git add pi_reference.txt test
git commit -m "Extend reference to 100K digits, add stress tests"
```

---

### Task 13: Initialize codescan and verify navigation

**Files:**
- No new files created

**Step 1: Initialize codescan**

Run: `codescan init`

**Step 2: Index the project**

Run: `codescan index`

**Step 3: Verify symbol extraction**

Run: `codescan symbols chudnovsky.ll`
Expected: Shows `compute_pi`, `split`, `base_case` as symbols.

**Step 4: Verify search**

Run: `codescan search "gmpz_mul"` and `codescan find-symbol split`
Expected: Finds relevant lines in chudnovsky.ll.

**Step 5: Commit codescan config if generated**

```bash
git add .codescan/ 2>/dev/null || true
git commit -m "Initialize codescan for IR navigation" 2>/dev/null || true
```

---

### Task 14: Add benchmark script

**Files:**
- Create: `bm`

**Step 1: Write benchmark script**

```bash
#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

echo "=== Building (optimized) ==="
./build -O2

echo
echo "=== Benchmarks ==="
echo "Timing uses CPU time (user+sys) via /usr/bin/time"

for digits in 1000 10000 100000 1000000; do
	printf "%-12s " "${digits} digits:"
	# Use /usr/bin/time for CPU time; discard pi output
	result=$( { /usr/bin/time -p nix develop -c ./pi "$digits" > /dev/null; } 2>&1 )
	user=$(echo "$result" | grep '^user' | awk '{print $2}')
	sys=$(echo "$result" | grep '^sys' | awk '{print $2}')
	total=$(echo "$user + $sys" | bc)
	echo "${total}s cpu (${user}s user, ${sys}s sys)"
done
```

**Step 2: Make executable and run**

Run: `chmod +x bm && nix develop -c bash -c './bm'`
Expected: Timing output for each digit count.

**Step 3: Commit**

```bash
git add bm
git commit -m "Add benchmark script"
```

---

## Summary of Tasks

| # | Task | Key Output |
|---|------|-----------|
| 1 | flake.nix | GMP + LLVM 19 + clang 19 dev environment |
| 2 | build script | Assembles .ll, compiles .c, links with GMP |
| 3 | pi_reference.txt | Known pi digits for test comparison |
| 4 | chudnovsky.ll stub | GMP declarations + build-pipeline-validating stub |
| 5 | main.c | CLI entry point, end-to-end build validation |
| 6 | test script | Digit-comparison test harness |
| 7 | base_case | Single Chudnovsky term in LLVM IR |
| 8 | split | Recursive binary splitting in LLVM IR |
| 9 | compute_pi | Full algorithm: split + mpf sqrt/div + string output |
| 10 | Output formatting | Truncate to exact digit count |
| 11 | Debug & fix | Make all tests green |
| 12 | 100K reference | Extended reference digits + stress tests |
| 13 | codescan | Initialize and verify IR navigation |
| 14 | Benchmark | Timing script for 1K → 1M digits |
