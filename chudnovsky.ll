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
declare void @__gmpz_set(ptr, ptr) nounwind
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

; === String constants ===
@.fmt_Ff = private unnamed_addr constant [6 x i8] c"%.*Ff\00", align 1

; === Constants ===
; C3_OVER_24 = 640320^3 / 24 = 10939058860032000
; A = 13591409
; B = 545140134

; ============================================================================
; base_case: compute single Chudnovsky term for index k
;   P, Q, T are pre-initialized mpz_t (caller owns init/clear)
;
;   k=0: P=1, Q=1, T=13591409
;   k>0: P=(6k-5)(2k-1)(6k-1), Q=k^3*C3_OVER_24, T=P*(A+B*k)*(-1)^k
; ============================================================================
define void @base_case(i64 %k, ptr %P, ptr %Q, ptr %T) nounwind {
entry:
  ; Alloca in entry block (required for correct codegen)
  %tmp_lin = alloca %struct.mpz, align 8

  %is_zero = icmp eq i64 %k, 0
  br i1 %is_zero, label %case_zero, label %case_nonzero

case_zero:
  call void @__gmpz_set_ui(ptr %P, i64 1)
  call void @__gmpz_set_ui(ptr %Q, i64 1)
  call void @__gmpz_set_si(ptr %T, i64 13591409)
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
  %k2_val = mul i64 %k, %k
  %k3_val = mul i64 %k2_val, %k
  call void @__gmpz_set_si(ptr %Q, i64 %k3_val)
  call void @__gmpz_mul_ui(ptr %Q, ptr %Q, i64 10939058860032000)

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
  br i1 %is_odd, label %negate, label %done

negate:
  call void @__gmpz_neg(ptr %T, ptr %T)
  br label %done

done:
  ret void
}

; ============================================================================
; split: recursive binary splitting over range [a, b)
;   P, Q, T are pre-initialized mpz_t (caller owns init/clear)
;
;   Merge: T = T_L*Q_R + P_L*T_R, P = P_L*P_R, Q = Q_L*Q_R
; ============================================================================
define void @split(i64 %a, i64 %b, ptr %P, ptr %Q, ptr %T) nounwind {
entry:
  ; All allocas in entry block (required for correct codegen with optimization)
  %Pr = alloca %struct.mpz, align 8
  %Qr = alloca %struct.mpz, align 8
  %Tr = alloca %struct.mpz, align 8
  %tmp_merge = alloca %struct.mpz, align 8

  %diff = sub i64 %b, %a
  %is_base = icmp eq i64 %diff, 1
  br i1 %is_base, label %do_base, label %do_split

do_base:
  call void @base_case(i64 %a, ptr %P, ptr %Q, ptr %T)
  br label %done

do_split:
  %sum = add i64 %a, %b
  %m = sdiv i64 %sum, 2

  call void @__gmpz_init(ptr %Pr)
  call void @__gmpz_init(ptr %Qr)
  call void @__gmpz_init(ptr %Tr)

  ; Left half: results go into P, Q, T
  call void @split(i64 %a, i64 %m, ptr %P, ptr %Q, ptr %T)

  ; Right half: results go into Pr, Qr, Tr
  call void @split(i64 %m, i64 %b, ptr %Pr, ptr %Qr, ptr %Tr)

  ; Merge
  call void @__gmpz_init(ptr %tmp_merge)
  call void @__gmpz_mul(ptr %tmp_merge, ptr %P, ptr %Tr)    ; tmp = P_L * T_R
  call void @__gmpz_mul(ptr %T, ptr %T, ptr %Qr)            ; T = T_L * Q_R
  call void @__gmpz_add(ptr %T, ptr %T, ptr %tmp_merge)     ; T = T_L*Q_R + P_L*T_R
  call void @__gmpz_mul(ptr %P, ptr %P, ptr %Pr)            ; P = P_L * P_R
  call void @__gmpz_mul(ptr %Q, ptr %Q, ptr %Qr)            ; Q = Q_L * Q_R

  call void @__gmpz_clear(ptr %tmp_merge)
  call void @__gmpz_clear(ptr %Pr)
  call void @__gmpz_clear(ptr %Qr)
  call void @__gmpz_clear(ptr %Tr)

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
;   pi = Q(0,N) * 426880 * sqrt(10005) / T(0,N)
;   N = digits/14 + 2, precision = digits*4 + 128 bits
; ============================================================================
define i64 @compute_pi(i64 %digits, ptr noalias nocapture %buf, i64 %buf_len) nounwind {
entry:
  ; N = digits / 14 + 2  (each Chudnovsky term gives ~14.18 digits)
  %d14 = sdiv i64 %digits, 14
  %N = add i64 %d14, 2

  ; Precision in bits (generous overestimate)
  %prec_base = mul i64 %digits, 4
  %prec = add i64 %prec_base, 128

  ; Integer binary splitting: split(0, N, P, Q, T)
  %P = alloca %struct.mpz, align 8
  %Q = alloca %struct.mpz, align 8
  %T = alloca %struct.mpz, align 8
  call void @__gmpz_init(ptr %P)
  call void @__gmpz_init(ptr %Q)
  call void @__gmpz_init(ptr %T)

  call void @split(i64 0, i64 %N, ptr %P, ptr %Q, ptr %T)

  ; Float conversion: pi = Q * 426880 * sqrt(10005) / T
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

  call void @__gmpf_set_ui(ptr %c10005_f, i64 10005)
  call void @__gmpf_sqrt(ptr %sqrt_f, ptr %c10005_f)

  call void @__gmpf_mul_ui(ptr %pi_f, ptr %q_f, i64 426880)
  call void @__gmpf_mul(ptr %pi_f, ptr %pi_f, ptr %sqrt_f)
  call void @__gmpf_div(ptr %pi_f, ptr %pi_f, ptr %t_f)

  ; Format to decimal string using gmp_sprintf("%.*Ff", digits+20, pi_f)
  ; Request extra digits to avoid rounding artifacts, then truncate below
  %tmp_buf_size = add i64 %digits, 256
  %tmp_buf = call ptr @malloc(i64 %tmp_buf_size)

  %digits_extra = add i64 %digits, 20
  %digits_extra_i32 = trunc i64 %digits_extra to i32
  %written_raw = call i32 (ptr, ptr, ...) @__gmp_sprintf(ptr %tmp_buf, ptr @.fmt_Ff, i32 %digits_extra_i32, ptr %pi_f)

  ; Copy exactly min(digits+2, written) chars to output buf
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
