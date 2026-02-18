; chudnovsky.ll — Hand-tuned binary splitting for Chudnovsky pi computation

target datalayout = "e-m:o-p270:32:32-p271:32:32-p272:64:64-i64:64-i128:128-n32:64-S128-Fn32"
target triple = "arm64-apple-macosx14.0.0"

;
; HAND OPTIMIZATIONS (things the C compiler fundamentally cannot do):
;
; 1. SEMANTIC STRENGTH REDUCTION
;    b^3 as native i64, C^3/24 pre-combined into single constant.
;
; 2. CROSS-PRODUCT PRECOMPUTATION
;    (2b-1)*(6b-5) and B*b+A computed as native i64.
;
; 3. BRANCHLESS NEGATION VIA GMP INTERNALS
;    Direct _mp_size field flip via branchless select.
;
; 4. CTZ-BASED BIT STRIPPING
;    @llvm.cttz.i64 replaces branch loop.
;
; 5. SINGLE-LOAD GLOBALS
;    Load @bs_stack/@top once per call, reuse everywhere.
;
; 6. MERGE-PHASE ADDMUL FUSION
;    mpz_addmul replaces separate mul+add (2 calls instead of 3).
;
; 7. DIRECT LIMB MANIPULATION
;    Q, G, T computed via native i128 arithmetic with direct writes to
;    mpz struct limb arrays. 0 GMP calls per base case.
;
; 8. DIFF==2 SPECIAL CASE
;    Two adjacent base cases computed inline with immediate merge,
;    eliminating 2 recursive calls and stack push/pop overhead.
;    ~352K occurrences at 10M digits.
;
; Result: 0 GMP calls per base case (was 14 in original C).
; Merge phase: 2 GMP calls per merge (was 3).
; diff==2 eliminates ~352K recursive call pairs.

; ============================================================================
; Type definitions (must match sieve.c struct layouts exactly)
; ============================================================================

%struct.__mpz_struct = type { i32, i32, ptr }
%struct.fac_s = type { i64, i64, ptr, ptr }
%struct.bs_level_t = type {
  [1 x %struct.__mpz_struct],
  [1 x %struct.__mpz_struct],
  [1 x %struct.__mpz_struct],
  %struct.fac_s,
  %struct.fac_s
}

; ============================================================================
; External globals (thread-local for parallel binary splitting)
; ============================================================================

@bs_stack = external dso_local thread_local global ptr, align 8
@top = external dso_local thread_local global i64, align 8

; ============================================================================
; GMP function declarations
; ============================================================================

declare void @__gmpz_mul(ptr, ptr, ptr) nounwind
declare void @__gmpz_addmul(ptr, ptr, ptr) nounwind

; ============================================================================
; C helper function declarations (defined in sieve.c)
; ============================================================================

declare void @fac_set_bp(ptr, i64, i64) nounwind
declare void @fac_mul_bp(ptr, i64, i64) nounwind
declare void @fac_mul_10005(ptr, i64) nounwind
declare void @fac_mul(ptr, ptr) nounwind
declare void @fac_remove_gcd(ptr, ptr, ptr, ptr) nounwind

; ============================================================================
; C entry points (defined in sieve.c)
; ============================================================================

declare i64 @compute_pi_c(i64, ptr, i64) nounwind

; ============================================================================
; LLVM intrinsics
; ============================================================================

declare i64 @llvm.cttz.i64(i64, i1 immarg) nounwind readnone
declare i64 @llvm.umax.i64(i64, i64) nounwind readnone

; ============================================================================
; bs_base_at: Compute one Chudnovsky base case term
;
; Computes Q, G, T, fp, fg for 1-indexed term %b, storing results at
; stack[%idx]. All arithmetic done via native i128 with direct limb writes.
;
; Marked alwaysinline: called from base_case and diff_two paths.
; ============================================================================

define internal void @bs_base_at(i64 %b, ptr %stk, i64 %idx) alwaysinline nounwind {
entry:
  %q1 = getelementptr inbounds %struct.bs_level_t, ptr %stk, i64 %idx
  %t1 = getelementptr inbounds %struct.bs_level_t, ptr %stk, i64 %idx, i32 1
  %g1 = getelementptr inbounds %struct.bs_level_t, ptr %stk, i64 %idx, i32 2
  %fp1 = getelementptr inbounds %struct.bs_level_t, ptr %stk, i64 %idx, i32 3
  %fg1 = getelementptr inbounds %struct.bs_level_t, ptr %stk, i64 %idx, i32 4

  ; ── Native i64 sub-expressions ──
  %b_sq = mul i64 %b, %b
  %b_cubed = mul i64 %b_sq, %b
  %b_x2 = shl i64 %b, 1
  %b2m1 = sub i64 %b_x2, 1
  %b_x6 = mul i64 %b, 6
  %b6m1 = sub i64 %b_x6, 1
  %b6m5 = sub i64 %b_x6, 5
  %g_partial = mul i64 %b2m1, %b6m5
  %Bb = mul i64 %b, 545140134
  %BbA = add i64 %Bb, 13591409

  ; ── Q = b^3 * C^3/24 via i128 ──
  %b_cubed_128 = zext i64 %b_cubed to i128
  %q_product = mul i128 %b_cubed_128, 10939058860032000
  %q_lo = trunc i128 %q_product to i64
  %q_hi128 = lshr i128 %q_product, 64
  %q_hi = trunc i128 %q_hi128 to i64
  %q_d_field = getelementptr inbounds %struct.__mpz_struct, ptr %q1, i32 0, i32 2
  %q_d = load ptr, ptr %q_d_field, align 8
  store i64 %q_lo, ptr %q_d, align 8
  %q_d_1 = getelementptr inbounds i64, ptr %q_d, i32 1
  store i64 %q_hi, ptr %q_d_1, align 8
  %q_has_hi = icmp ne i64 %q_hi, 0
  %q_size = select i1 %q_has_hi, i32 2, i32 1
  %q_size_ptr = getelementptr inbounds %struct.__mpz_struct, ptr %q1, i32 0, i32 1
  store i32 %q_size, ptr %q_size_ptr, align 4

  ; ── G = (2b-1)*(6b-5)*(6b-1) via i128 ──
  %g_partial_128 = zext i64 %g_partial to i128
  %b6m1_128 = zext i64 %b6m1 to i128
  %g_product = mul i128 %g_partial_128, %b6m1_128
  %g_lo = trunc i128 %g_product to i64
  %g_hi128 = lshr i128 %g_product, 64
  %g_hi = trunc i128 %g_hi128 to i64
  %g_d_field = getelementptr inbounds %struct.__mpz_struct, ptr %g1, i32 0, i32 2
  %g_d = load ptr, ptr %g_d_field, align 8
  store i64 %g_lo, ptr %g_d, align 8
  %g_d_1 = getelementptr inbounds i64, ptr %g_d, i32 1
  store i64 %g_hi, ptr %g_d_1, align 8
  %g_has_hi = icmp ne i64 %g_hi, 0
  %g_size = select i1 %g_has_hi, i32 2, i32 1
  %g_size_ptr = getelementptr inbounds %struct.__mpz_struct, ptr %g1, i32 0, i32 1
  store i32 %g_size, ptr %g_size_ptr, align 4

  ; ── T = (B*b+A) * G * (-1)^b via i128 ──
  %BbA_128 = zext i64 %BbA to i128
  %t_product = mul i128 %BbA_128, %g_product
  %t_lo = trunc i128 %t_product to i64
  %t_hi128 = lshr i128 %t_product, 64
  %t_hi = trunc i128 %t_hi128 to i64
  %t_d_field = getelementptr inbounds %struct.__mpz_struct, ptr %t1, i32 0, i32 2
  %t_d = load ptr, ptr %t_d_field, align 8
  store i64 %t_lo, ptr %t_d, align 8
  %t_d_1 = getelementptr inbounds i64, ptr %t_d, i32 1
  store i64 %t_hi, ptr %t_d_1, align 8
  %t_has_hi = icmp ne i64 %t_hi, 0
  %t_abs_size = select i1 %t_has_hi, i32 2, i32 1
  %t_neg_size = sub i32 0, %t_abs_size
  %b_low = and i64 %b, 1
  %is_odd = icmp ne i64 %b_low, 0
  %t_final_size = select i1 %is_odd, i32 %t_neg_size, i32 %t_abs_size
  %t_size_ptr = getelementptr inbounds %struct.__mpz_struct, ptr %t1, i32 0, i32 1
  store i32 %t_final_size, ptr %t_size_ptr, align 4

  ; ── CTZ strip + factorization ──
  %tz = call i64 @llvm.cttz.i64(i64 %b, i1 true)
  %i = lshr i64 %b, %tz
  call void @fac_set_bp(ptr %fp1, i64 %i, i64 3)
  call void @fac_mul_10005(ptr %fp1, i64 3)
  %pow_field = getelementptr inbounds %struct.fac_s, ptr %fp1, i32 0, i32 3
  %pow_ptr = load ptr, ptr %pow_field, align 8
  %pow0 = load i64, ptr %pow_ptr, align 8
  %pow0_dec = sub i64 %pow0, 1
  store i64 %pow0_dec, ptr %pow_ptr, align 8
  call void @fac_set_bp(ptr %fg1, i64 %b2m1, i64 1)
  call void @fac_mul_bp(ptr %fg1, i64 %b6m1, i64 1)
  call void @fac_mul_bp(ptr %fg1, i64 %b6m5, i64 1)
  ret void
}

; ============================================================================
; bs: Hand-tuned binary splitting over Chudnovsky terms (a, b]
;
; Computes Q, T, G for terms in the half-open interval (a, b].
; Results stored in bs_stack[top].{q, t, g, fp, fg}.
;
; Parameters:
;   %a     - left bound (exclusive)
;   %b     - right bound (inclusive)
;   %gflag - 1 = maintain G product (needed by caller's merge), 0 = skip
;   %level - recursion depth (GCD removal enabled at depth >= 4)
; ============================================================================

define void @bs(i64 %a, i64 %b, i32 %gflag, i64 %level) nounwind {
entry:
  %diff = sub i64 %b, %a
  %is_base = icmp eq i64 %diff, 1
  br i1 %is_base, label %base_case, label %check_two

; ──────────────────────────────────────────────────────────────────────────────
; BASE CASE (diff==1): Single term
; ──────────────────────────────────────────────────────────────────────────────

base_case:
  %bc_stack = load ptr, ptr @bs_stack, align 8
  %bc_top = load i64, ptr @top, align 8
  call void @bs_base_at(i64 %b, ptr %bc_stack, i64 %bc_top)
  ret void

; ──────────────────────────────────────────────────────────────────────────────
; DIFF==2 SPECIAL CASE: Two adjacent terms with inline merge
;
; Eliminates 2 recursive calls + stack push/pop.
; ~352K occurrences at 10M digits.
; ──────────────────────────────────────────────────────────────────────────────

check_two:
  %is_two = icmp eq i64 %diff, 2
  br i1 %is_two, label %diff_two, label %recursive

diff_two:
  %dt_stack = load ptr, ptr @bs_stack, align 8
  %dt_top = load i64, ptr @top, align 8
  %dt_top1 = add i64 %dt_top, 1

  ; Left base case: term (a+1) at stack[top]
  %dt_b_left = add i64 %a, 1
  call void @bs_base_at(i64 %dt_b_left, ptr %dt_stack, i64 %dt_top)

  ; Right base case: term b (= a+2) at stack[top+1]
  call void @bs_base_at(i64 %b, ptr %dt_stack, i64 %dt_top1)

  ; ── Merge pointers ──
  %dt_q1 = getelementptr inbounds %struct.bs_level_t, ptr %dt_stack, i64 %dt_top
  %dt_t1 = getelementptr inbounds %struct.bs_level_t, ptr %dt_stack, i64 %dt_top, i32 1
  %dt_g1 = getelementptr inbounds %struct.bs_level_t, ptr %dt_stack, i64 %dt_top, i32 2
  %dt_fp1 = getelementptr inbounds %struct.bs_level_t, ptr %dt_stack, i64 %dt_top, i32 3
  %dt_fg1 = getelementptr inbounds %struct.bs_level_t, ptr %dt_stack, i64 %dt_top, i32 4

  %dt_q2 = getelementptr inbounds %struct.bs_level_t, ptr %dt_stack, i64 %dt_top1
  %dt_t2 = getelementptr inbounds %struct.bs_level_t, ptr %dt_stack, i64 %dt_top1, i32 1
  %dt_g2 = getelementptr inbounds %struct.bs_level_t, ptr %dt_stack, i64 %dt_top1, i32 2
  %dt_fp2 = getelementptr inbounds %struct.bs_level_t, ptr %dt_stack, i64 %dt_top1, i32 3
  %dt_fg2 = getelementptr inbounds %struct.bs_level_t, ptr %dt_stack, i64 %dt_top1, i32 4

  ; GCD removal at depth >= 4
  %dt_deep = icmp sge i64 %level, 4
  br i1 %dt_deep, label %dt_gcd, label %dt_merge

dt_gcd:
  call void @fac_remove_gcd(ptr %dt_q2, ptr %dt_fp2, ptr %dt_g1, ptr %dt_fg1)
  br label %dt_merge

dt_merge:
  ; T = T_left * Q_right + G_left * T_right
  call void @__gmpz_mul(ptr %dt_t1, ptr %dt_t1, ptr %dt_q2)
  call void @__gmpz_addmul(ptr %dt_t1, ptr %dt_t2, ptr %dt_g1)
  ; Q = Q_left * Q_right
  call void @__gmpz_mul(ptr %dt_q1, ptr %dt_q1, ptr %dt_q2)
  ; fp merge
  call void @fac_mul(ptr %dt_fp1, ptr %dt_fp2)
  ; Conditional G merge
  %dt_want_g = icmp ne i32 %gflag, 0
  br i1 %dt_want_g, label %dt_merge_g, label %dt_done

dt_merge_g:
  call void @__gmpz_mul(ptr %dt_g1, ptr %dt_g1, ptr %dt_g2)
  call void @fac_mul(ptr %dt_fg1, ptr %dt_fg2)
  br label %dt_done

dt_done:
  ret void

; ──────────────────────────────────────────────────────────────────────────────
; RECURSIVE CASE (diff >= 3): Divide and conquer with asymmetric split
; ──────────────────────────────────────────────────────────────────────────────

recursive:
  ; ── Compute asymmetric midpoint (0.5224 ratio) ──
  %diff_f = uitofp i64 %diff to double
  %mid_f = fmul double %diff_f, 5.224000e-01
  %mid_offset = fptoui double %mid_f to i64
  %mid_at_least_1 = call i64 @llvm.umax.i64(i64 %mid_offset, i64 1)
  %mid_unclamped = add i64 %mid_at_least_1, %a
  %bm1 = sub i64 %b, 1
  %in_range = icmp ult i64 %mid_unclamped, %b
  %mid = select i1 %in_range, i64 %mid_unclamped, i64 %bm1

  %level1 = add i64 %level, 1

  ; ── Left child: bs(a, mid, 1, level+1) ──
  call void @bs(i64 %a, i64 %mid, i32 1, i64 %level1)

  ; ── top++ ──
  %old_top = load i64, ptr @top, align 8
  %inc_top = add i64 %old_top, 1
  store i64 %inc_top, ptr @top, align 8

  ; ── Right child: bs(mid, b, gflag, level+1) ──
  call void @bs(i64 %mid, i64 %b, i32 %gflag, i64 %level1)

  ; ── top-- ──
  %post_top = load i64, ptr @top, align 8
  %dec_top = sub i64 %post_top, 1
  store i64 %dec_top, ptr @top, align 8

  ; ════════════════════════════════════════════════════════════════════════════
  ; MERGE PHASE
  ; ════════════════════════════════════════════════════════════════════════════

  %m_stack = load ptr, ptr @bs_stack, align 8

  ; ── Left child pointers (level = dec_top) ──
  %m_q1 = getelementptr inbounds %struct.bs_level_t, ptr %m_stack, i64 %dec_top
  %m_t1 = getelementptr inbounds %struct.bs_level_t, ptr %m_stack, i64 %dec_top, i32 1
  %m_g1 = getelementptr inbounds %struct.bs_level_t, ptr %m_stack, i64 %dec_top, i32 2
  %m_fp1 = getelementptr inbounds %struct.bs_level_t, ptr %m_stack, i64 %dec_top, i32 3
  %m_fg1 = getelementptr inbounds %struct.bs_level_t, ptr %m_stack, i64 %dec_top, i32 4

  ; ── Right child pointers (level = post_top = dec_top + 1) ──
  %m_q2 = getelementptr inbounds %struct.bs_level_t, ptr %m_stack, i64 %post_top
  %m_t2 = getelementptr inbounds %struct.bs_level_t, ptr %m_stack, i64 %post_top, i32 1
  %m_g2 = getelementptr inbounds %struct.bs_level_t, ptr %m_stack, i64 %post_top, i32 2
  %m_fp2 = getelementptr inbounds %struct.bs_level_t, ptr %m_stack, i64 %post_top, i32 3
  %m_fg2 = getelementptr inbounds %struct.bs_level_t, ptr %m_stack, i64 %post_top, i32 4

  ; ── GCD removal at depth >= 4 ──
  %deep_enough = icmp sge i64 %level, 4
  br i1 %deep_enough, label %do_gcd, label %merge

do_gcd:
  call void @fac_remove_gcd(ptr %m_q2, ptr %m_fp2, ptr %m_g1, ptr %m_fg1)
  br label %merge

merge:
  ; T = T_left * Q_right + G_left * T_right (addmul fusion)
  call void @__gmpz_mul(ptr %m_t1, ptr %m_t1, ptr %m_q2)
  call void @__gmpz_addmul(ptr %m_t1, ptr %m_t2, ptr %m_g1)
  ; Q = Q_left * Q_right
  call void @__gmpz_mul(ptr %m_q1, ptr %m_q1, ptr %m_q2)
  ; fp merge
  call void @fac_mul(ptr %m_fp1, ptr %m_fp2)
  ; Conditional G merge
  %want_g = icmp ne i32 %gflag, 0
  br i1 %want_g, label %merge_g, label %done

merge_g:
  call void @__gmpz_mul(ptr %m_g1, ptr %m_g1, ptr %m_g2)
  call void @fac_mul(ptr %m_fg1, ptr %m_fg2)
  br label %done

done:
  ret void
}

; ============================================================================
; binary_split: Entry point for the binary splitting phase
; ============================================================================

define void @binary_split(i64 %N) nounwind {
entry:
  store i64 0, ptr @top, align 8
  call void @bs(i64 0, i64 %N, i32 0, i64 0)
  ret void
}

; ============================================================================
; compute_pi: Public entry point (called from main.c)
; ============================================================================

define i64 @compute_pi(i64 %digits, ptr noalias nocapture %buf, i64 %buf_len) nounwind {
entry:
  %result = call i64 @compute_pi_c(i64 %digits, ptr %buf, i64 %buf_len)
  ret i64 %result
}
