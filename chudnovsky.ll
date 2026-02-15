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
@.fmt_pi = private unnamed_addr constant [5 x i8] c"%Ff\0A\00", align 1

; === Public entry point (STUB) ===
; Computes pi to the requested number of decimal digits.
; Writes the decimal string to buf (must be at least digits+10 bytes).
; Returns the number of characters written.
define i64 @compute_pi(i64 %digits, ptr noalias nocapture %buf, i64 %buf_len) nounwind {
entry:
  ; --- STUB: just write "3.14" for build pipeline validation ---
  ; This will be replaced with the real algorithm in Task 9.
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
