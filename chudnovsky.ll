; chudnovsky.ll — Binary splitting Chudnovsky pi computation with GCD optimization
; Core algorithm in LLVM IR, calls GMP for arbitrary-precision arithmetic.
; Prime sieve helpers in sieve.c for factorized-number operations.

; === GMP struct types ===
; __mpz_struct = { int _mp_alloc, int _mp_size, mp_limb_t* _mp_d }
%struct.mpz = type { i32, i32, ptr }

; __mpf_struct = { int _mp_prec, int _mp_size, mp_exp_t _mp_exp, mp_limb_t* _mp_d }
%struct.mpf = type { i32, i32, i64, ptr }

; fac_s = { unsigned long max_facs, num_facs, *fac, *pow }
%struct.fac = type { i64, i64, ptr, ptr }

; === GMP integer function declarations ===
declare void @__gmpz_init(ptr) nounwind
declare void @__gmpz_clear(ptr) nounwind
declare void @__gmpz_set_ui(ptr, i64) nounwind
declare void @__gmpz_set_si(ptr, i64) nounwind
declare void @__gmpz_set(ptr, ptr) nounwind
declare void @__gmpz_mul(ptr, ptr, ptr) nounwind
declare void @__gmpz_mul_ui(ptr, ptr, i64) nounwind
declare void @__gmpz_mul_si(ptr, ptr, i64) nounwind
declare void @__gmpz_add(ptr, ptr, ptr) nounwind
declare void @__gmpz_sub(ptr, ptr, ptr) nounwind
declare void @__gmpz_addmul_ui(ptr, ptr, i64) nounwind
declare void @__gmpz_neg(ptr, ptr) nounwind
declare void @__gmpz_tdiv_q(ptr, ptr, ptr) nounwind

; === GMP float function declarations ===
declare void @__gmpf_init2(ptr, i64) nounwind
declare void @__gmpf_clear(ptr) nounwind
declare void @__gmpf_set_z(ptr, ptr) nounwind
declare void @__gmpf_set_ui(ptr, i64) nounwind
declare void @__gmpf_div(ptr, ptr, ptr) nounwind
declare void @__gmpf_mul(ptr, ptr, ptr) nounwind
declare void @__gmpf_mul_ui(ptr, ptr, i64) nounwind
declare void @__gmpf_sqrt(ptr, ptr) nounwind

; === GMP output ===
declare i32 @__gmp_sprintf(ptr, ptr, ...) nounwind

; === C stdlib ===
declare ptr @malloc(i64) nounwind
declare void @free(ptr) nounwind

; === Sieve helpers (sieve.c) ===
declare void @sieve_build(i64) nounwind
declare void @sieve_free() nounwind
declare void @sieve_helpers_init() nounwind
declare void @sieve_helpers_free() nounwind
declare void @fac_init(ptr) nounwind
declare void @fac_clear(ptr) nounwind
declare void @fac_reset(ptr) nounwind
declare void @fac_set_bp(ptr, i64, i64) nounwind
declare void @fac_mul_bp(ptr, i64, i64) nounwind
declare void @fac_mul(ptr, ptr) nounwind
declare void @fac_remove_gcd(ptr, ptr, ptr, ptr) nounwind

; === String constants ===
@.fmt_Ff = private unnamed_addr constant [6 x i8] c"%.*Ff\00", align 1

; === Constants ===
; C3_OVER_24 = 640320^3 / 24 = 10939058860032000
; A = 13591409, B = 545140134, C = 640320
; C_FACS = 3*5*23*29 = 10005 (product of odd prime factors of C)

; ============================================================================
; base_case: compute single Chudnovsky term for index k
;   P, Q, T, G are pre-initialized mpz_t
;   fpQ, fgP are pre-initialized fac_s
;
;   k=0: P=1, Q=1, T=A=13591409, G=1, fpQ={}, fgP={}
;   k>0: P=(6k-5)(2k-1)(6k-1), Q=k^3*C3_OVER_24
;        T=P*(A+B*k)*(-1)^k, G=P
;        fpQ=factorized(Q), fgP=factorized(P)
; ============================================================================
define void @base_case(i64 %k, ptr %P, ptr %Q, ptr %T, ptr %G, ptr %fpQ, ptr %fgP) nounwind {
entry:
  %tmp_lin = alloca %struct.mpz, align 8

  %is_zero = icmp eq i64 %k, 0
  br i1 %is_zero, label %case_zero, label %case_nonzero

case_zero:
  call void @__gmpz_set_ui(ptr %P, i64 1)
  call void @__gmpz_set_ui(ptr %Q, i64 1)
  call void @__gmpz_set_si(ptr %T, i64 13591409)
  call void @__gmpz_set_ui(ptr %G, i64 1)
  call void @fac_reset(ptr %fpQ)
  call void @fac_reset(ptr %fgP)
  br label %done

case_nonzero:
  ; P = (6k-5) * (2k-1) * (6k-1)
  %k6 = mul i64 %k, 6
  %k6m5 = sub i64 %k6, 5
  %k2 = mul i64 %k, 2
  %k2m1 = sub i64 %k2, 1
  %k6m1 = sub i64 %k6, 1

  call void @__gmpz_set_si(ptr %P, i64 %k6m5)
  call void @__gmpz_mul_si(ptr %P, ptr %P, i64 %k2m1)
  call void @__gmpz_mul_si(ptr %P, ptr %P, i64 %k6m1)

  ; Q = k^3 * C3_OVER_24
  ; Compute Q using individual multiplications to keep within i64 range
  call void @__gmpz_set_ui(ptr %Q, i64 %k)
  call void @__gmpz_mul_ui(ptr %Q, ptr %Q, i64 %k)
  call void @__gmpz_mul_ui(ptr %Q, ptr %Q, i64 %k)
  ; C^3/24 = 640320^3/24 = (C/24)*(C/24)*(C*24) but we use factored form:
  ; 640320 = 2^7 * 3 * 5 * 23 * 29 = 2^7 * 10005
  ; C^3/24 = 2^21 * 10005^3 / 24 = 2^21 * 10005^3 / (2^3 * 3)
  ;        = 2^18 * 10005^3 / 3 = 2^18 * 10005^2 * 3335
  ; Actually just use the constant directly:
  call void @__gmpz_mul_ui(ptr %Q, ptr %Q, i64 10939058860032000)

  ; G = P (copy)
  call void @__gmpz_set(ptr %G, ptr %P)

  ; T = P * (A + B*k)
  %bk = mul i64 545140134, %k
  %a_plus_bk = add i64 13591409, %bk
  call void @__gmpz_init(ptr %tmp_lin)
  call void @__gmpz_set_si(ptr %tmp_lin, i64 %a_plus_bk)
  call void @__gmpz_mul(ptr %T, ptr %P, ptr %tmp_lin)
  call void @__gmpz_clear(ptr %tmp_lin)

  ; If k is odd, negate T
  %k_and_1 = and i64 %k, 1
  %is_odd = icmp ne i64 %k_and_1, 0
  br i1 %is_odd, label %negate, label %set_facs

negate:
  call void @__gmpz_neg(ptr %T, ptr %T)
  br label %set_facs

set_facs:
  ; fpQ = factorized form of Q = k^3 * C3_OVER_24
  ; Strip factors of 2 from k to get k_odd
  br label %strip_loop

strip_loop:
  %kv = phi i64 [ %k, %set_facs ], [ %kv_shr, %strip_continue ]
  %kv_and_1 = and i64 %kv, 1
  %kv_even = icmp eq i64 %kv_and_1, 0
  br i1 %kv_even, label %strip_continue, label %strip_done

strip_continue:
  %kv_shr = lshr i64 %kv, 1
  br label %strip_loop

strip_done:
  ; kv is now the odd part of k
  ; fpQ = (k_odd)^3 * (3*5*23*29)^3, then decrement first power by 1
  ; The constant 3*5*23*29 = 10005
  call void @fac_set_bp(ptr %fpQ, i64 %kv, i64 3)
  call void @fac_mul_bp(ptr %fpQ, i64 10005, i64 3)
  ; Adjust: C^3/24 = (2^7 * 10005)^3 / 24 = 2^21 * 10005^3 / 24
  ; = 2^21 * 10005^3 / (8*3) = 2^18 * 10005^3 / 3
  ; We tracked 10005^3 but the actual constant includes /3.
  ; In gmp-chudnovsky.c they use 3*5*23*29 and decrement pow[0] by 1.
  ; pow[0] is for the smallest prime factor of 10005 = 3*5*23*29,
  ; which is 3. So pow[0]-- means dividing by 3^1, giving us the /3.
  ; Access fpQ->pow[0] and decrement
  %fpQ_pow_ptr = getelementptr %struct.fac, ptr %fpQ, i32 0, i32 3
  %fpQ_pow = load ptr, ptr %fpQ_pow_ptr
  %fpQ_pow0 = load i64, ptr %fpQ_pow
  %fpQ_pow0_dec = sub i64 %fpQ_pow0, 1
  store i64 %fpQ_pow0_dec, ptr %fpQ_pow

  ; fgP = factorized form of P = (2k-1)(6k-1)(6k-5)
  call void @fac_set_bp(ptr %fgP, i64 %k2m1, i64 1)
  call void @fac_mul_bp(ptr %fgP, i64 %k6m1, i64 1)
  call void @fac_mul_bp(ptr %fgP, i64 %k6m5, i64 1)
  br label %done

done:
  ret void
}

; ============================================================================
; split: recursive binary splitting over range [a, b) with GCD optimization
;   P, Q, T, G are pre-initialized mpz_t
;   fpQ, fgP are pre-initialized fac_s
;   level: recursion depth from top (0 = initial call)
;   gflag: 1 = maintain G value, 0 = skip (top-level optimization)
;
;   Merge: remove GCD(Q_R, P_L) first, then
;     T = T_L*Q_R + G_L*T_R, Q = Q_L*Q_R, P = P_L*P_R
;     G = G_L*G_R (if gflag), fpQ/fgP updated
; ============================================================================
define void @split(i64 %a, i64 %b, ptr %P, ptr %Q, ptr %T, ptr %G, ptr %fpQ, ptr %fgP, i64 %level, i64 %gflag) nounwind {
entry:
  ; All allocas in entry block
  %Pr = alloca %struct.mpz, align 8
  %Qr = alloca %struct.mpz, align 8
  %Tr = alloca %struct.mpz, align 8
  %Gr = alloca %struct.mpz, align 8
  %fpQr = alloca %struct.fac, align 8
  %fgPr = alloca %struct.fac, align 8
  %tmp_merge = alloca %struct.mpz, align 8

  %diff = sub i64 %b, %a
  %is_base = icmp eq i64 %diff, 1
  br i1 %is_base, label %do_base, label %do_split

do_base:
  call void @base_case(i64 %a, ptr %P, ptr %Q, ptr %T, ptr %G, ptr %fpQ, ptr %fgP)
  br label %done

do_split:
  %sum = add i64 %a, %b
  %m = sdiv i64 %sum, 2

  call void @__gmpz_init(ptr %Pr)
  call void @__gmpz_init(ptr %Qr)
  call void @__gmpz_init(ptr %Tr)
  call void @__gmpz_init(ptr %Gr)
  call void @fac_init(ptr %fpQr)
  call void @fac_init(ptr %fgPr)

  ; Left half: always maintain G (gflag=1) since merge needs G_L
  %next_level = add i64 %level, 1
  call void @split(i64 %a, i64 %m, ptr %P, ptr %Q, ptr %T, ptr %G, ptr %fpQ, ptr %fgP, i64 %next_level, i64 1)

  ; Right half: pass through gflag
  call void @split(i64 %m, i64 %b, ptr %Pr, ptr %Qr, ptr %Tr, ptr %Gr, ptr %fpQr, ptr %fgPr, i64 %next_level, i64 %gflag)

  ; GCD removal: at level >= 4, remove GCD(Q_R, G_L) before multiplying
  ; This makes subsequent multiplications faster on smaller numbers
  %do_gcd = icmp sge i64 %level, 4
  br i1 %do_gcd, label %gcd_remove, label %merge

gcd_remove:
  ; fac_remove_gcd(Q_R, fpQ_R, G_L, fgP_L)
  ; This divides Q_R and G_L by their common prime factors
  call void @fac_remove_gcd(ptr %Qr, ptr %fpQr, ptr %G, ptr %fgP)
  br label %merge

merge:
  ; T = T_L * Q_R + G_L * T_R
  call void @__gmpz_init(ptr %tmp_merge)
  call void @__gmpz_mul(ptr %tmp_merge, ptr %G, ptr %Tr)     ; tmp = G_L * T_R
  call void @__gmpz_mul(ptr %T, ptr %T, ptr %Qr)             ; T = T_L * Q_R
  call void @__gmpz_add(ptr %T, ptr %T, ptr %tmp_merge)      ; T = T_L*Q_R + G_L*T_R

  ; Q = Q_L * Q_R  (Q_R may be reduced by GCD)
  call void @__gmpz_mul(ptr %Q, ptr %Q, ptr %Qr)

  ; P = P_L * P_R  (P is NOT affected by GCD removal)
  call void @__gmpz_mul(ptr %P, ptr %P, ptr %Pr)

  ; Merge factorized forms
  call void @fac_mul(ptr %fpQ, ptr %fpQr)

  ; Conditionally maintain G and fgP
  %need_g = icmp ne i64 %gflag, 0
  br i1 %need_g, label %merge_g, label %cleanup

merge_g:
  call void @__gmpz_mul(ptr %G, ptr %G, ptr %Gr)
  call void @fac_mul(ptr %fgP, ptr %fgPr)
  br label %cleanup

cleanup:
  call void @__gmpz_clear(ptr %tmp_merge)
  call void @__gmpz_clear(ptr %Pr)
  call void @__gmpz_clear(ptr %Qr)
  call void @__gmpz_clear(ptr %Tr)
  call void @__gmpz_clear(ptr %Gr)
  call void @fac_clear(ptr %fpQr)
  call void @fac_clear(ptr %fgPr)
  br label %done

done:
  ret void
}

; ============================================================================
; compute_pi: public entry point
;   Computes pi to the requested number of decimal digits.
;   Writes the decimal string "3.xxxxx" to buf.
;   Returns the number of characters written.
;
;   pi = (Q * C/D * sqrt(C)) / (T + Q*A)
;   where C=640320, D=12, A=13591409
;   N = digits/14 + 2, precision = digits*4 + 128 bits
; ============================================================================
define i64 @compute_pi(i64 %digits, ptr noalias nocapture %buf, i64 %buf_len) nounwind {
entry:
  ; N = digits / 14 + 2  (each Chudnovsky term gives ~14.18 digits)
  %d14 = sdiv i64 %digits, 14
  %N = add i64 %d14, 2

  ; Precision in bits
  %prec_base = mul i64 %digits, 4
  %prec = add i64 %prec_base, 128

  ; Build prime sieve: size = max(10005+1, N*6)
  %sieve_from_n = mul i64 %N, 6
  %sieve_min = icmp sgt i64 %sieve_from_n, 10006
  %sieve_size = select i1 %sieve_min, i64 %sieve_from_n, i64 10006
  call void @sieve_build(i64 %sieve_size)
  call void @sieve_helpers_init()

  ; Allocate P, Q, T, G and their factorized forms
  %P = alloca %struct.mpz, align 8
  %Q = alloca %struct.mpz, align 8
  %T = alloca %struct.mpz, align 8
  %G = alloca %struct.mpz, align 8
  %fpQ = alloca %struct.fac, align 8
  %fgP = alloca %struct.fac, align 8

  call void @__gmpz_init(ptr %P)
  call void @__gmpz_init(ptr %Q)
  call void @__gmpz_init(ptr %T)
  call void @__gmpz_init(ptr %G)
  call void @fac_init(ptr %fpQ)
  call void @fac_init(ptr %fgP)

  ; Binary splitting: gflag=0 at top level (don't need G after this)
  call void @split(i64 0, i64 %N, ptr %P, ptr %Q, ptr %T, ptr %G, ptr %fpQ, ptr %fgP, i64 0, i64 0)

  ; Free sieve resources (no longer needed)
  call void @fac_clear(ptr %fpQ)
  call void @fac_clear(ptr %fgP)
  call void @__gmpz_clear(ptr %G)
  call void @sieve_helpers_free()
  call void @sieve_free()

  ; Final formula: pi = Q * 426880 * sqrt(10005) / T
  ; Our binary splitting includes k=0, so T already contains the A term.
  ; (gmp-chudnovsky.c starts from k=1 and adds A back separately.)
  ; With GCD optimization, Q and T are reduced by the same cumulative factor,
  ; preserving their ratio.

  ; Float conversion
  %pi_f = alloca %struct.mpf, align 8
  %q_f = alloca %struct.mpf, align 8
  %t_f = alloca %struct.mpf, align 8
  %sqrt_f = alloca %struct.mpf, align 8
  %c10005_f = alloca %struct.mpf, align 8

  call void @__gmpf_init2(ptr %pi_f, i64 %prec)
  call void @__gmpf_init2(ptr %q_f, i64 %prec)
  call void @__gmpf_init2(ptr %t_f, i64 %prec)
  call void @__gmpf_init2(ptr %sqrt_f, i64 %prec)
  call void @__gmpf_init2(ptr %c10005_f, i64 %prec)

  call void @__gmpf_set_z(ptr %q_f, ptr %Q)
  call void @__gmpf_set_z(ptr %t_f, ptr %T)

  ; sqrt(10005)
  call void @__gmpf_set_ui(ptr %c10005_f, i64 10005)
  call void @__gmpf_sqrt(ptr %sqrt_f, ptr %c10005_f)

  ; pi = Q * 426880 * sqrt(10005) / T
  call void @__gmpf_mul_ui(ptr %pi_f, ptr %q_f, i64 426880)
  call void @__gmpf_mul(ptr %pi_f, ptr %pi_f, ptr %sqrt_f)
  call void @__gmpf_div(ptr %pi_f, ptr %pi_f, ptr %t_f)

  ; Format to decimal string
  %tmp_buf_size = add i64 %digits, 256
  %tmp_buf = call ptr @malloc(i64 %tmp_buf_size)

  %digits_extra = add i64 %digits, 20
  %digits_extra_i32 = trunc i64 %digits_extra to i32
  %written_raw = call i32 (ptr, ptr, ...) @__gmp_sprintf(ptr %tmp_buf, ptr @.fmt_Ff, i32 %digits_extra_i32, ptr %pi_f)

  ; Copy exactly digits+2 chars ("3." + digits) to output buf
  %copy_len = add i64 %digits, 2
  br label %copy_loop

copy_loop:
  %i = phi i64 [ 0, %entry ], [ %i_next, %copy_continue ]
  %cmp_len = icmp uge i64 %i, %copy_len
  br i1 %cmp_len, label %copy_done, label %copy_check_src

copy_check_src:
  %written_i64 = sext i32 %written_raw to i64
  %cmp_src = icmp uge i64 %i, %written_i64
  br i1 %cmp_src, label %copy_done, label %copy_continue

copy_continue:
  %src_ptr = getelementptr i8, ptr %tmp_buf, i64 %i
  %dst_ptr = getelementptr i8, ptr %buf, i64 %i
  %ch = load i8, ptr %src_ptr
  store i8 %ch, ptr %dst_ptr
  %i_next = add i64 %i, 1
  br label %copy_loop

copy_done:
  %term_ptr = getelementptr i8, ptr %buf, i64 %i
  store i8 0, ptr %term_ptr

  call void @free(ptr %tmp_buf)

  ; Cleanup
  call void @__gmpf_clear(ptr %pi_f)
  call void @__gmpf_clear(ptr %q_f)
  call void @__gmpf_clear(ptr %t_f)
  call void @__gmpf_clear(ptr %sqrt_f)
  call void @__gmpf_clear(ptr %c10005_f)

  call void @__gmpz_clear(ptr %P)
  call void @__gmpz_clear(ptr %Q)
  call void @__gmpz_clear(ptr %T)

  ret i64 %i
}
