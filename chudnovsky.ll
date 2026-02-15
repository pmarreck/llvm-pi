; chudnovsky.ll — Binary splitting Chudnovsky pi computation with GCD optimization
; Core algorithm orchestrated from LLVM IR, heavy lifting in sieve.c:
; - Prime sieve + factorized-number GCD optimization
; - Iterative binary splitting with pre-allocated stacks
; - Newton's method sqrt for final computation

declare i64 @compute_pi_c(i64, ptr, i64) nounwind

; ============================================================================
; compute_pi: public entry point (called from main.c)
;   Delegates to compute_pi_c in sieve.c which handles:
;   - Sieve construction, binary splitting, Newton sqrt, formatting
; ============================================================================
define i64 @compute_pi(i64 %digits, ptr noalias nocapture %buf, i64 %buf_len) nounwind {
entry:
  %result = call i64 @compute_pi_c(i64 %digits, ptr %buf, i64 %buf_len)
  ret i64 %result
}
