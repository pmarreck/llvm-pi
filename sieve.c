/* sieve.c — Prime sieve and factorized-number helpers for Chudnovsky GCD optimization.
 *
 * Adapted from gmp-chudnovsky.c by Hanhong Xue.
 * Called from chudnovsky.ll via external function declarations.
 *
 * The key optimization: at each merge step in binary splitting, we remove
 * the GCD of Q_right and P_left before multiplying. This makes the numbers
 * smaller, speeding up subsequent GMP multiplications.
 *
 * To do this efficiently, we maintain the prime factorization of P and Q
 * alongside their mpz_t values, enabling fast GCD computation via sorted
 * prime-factor list intersection.
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

/* ── Sieve data structure ─────────────────────────────────────────────── */

typedef struct {
	long fac;   /* smallest prime factor */
	long pow;   /* its power */
	long nxt;   /* index of next factor in chain (after dividing out this one) */
} sieve_t;

static sieve_t *sieve;
static long sieve_size;

/* Build prime factorization sieve for odd numbers up to n.
 * After building, sieve[k/2] gives the smallest prime factor of k,
 * its power, and a link to the remaining factorization. */
static void sieve_build(long n)
{
	long m, i, j, k;

	sieve_size = n;
	sieve = (sieve_t *)malloc(sizeof(sieve_t) * (n / 2 + 1));
	memset(sieve, 0, sizeof(sieve_t) * (n / 2 + 1));
	m = (long)sqrt((double)n);

	sieve[1 / 2].fac = 1;
	sieve[1 / 2].pow = 1;

	for (i = 3; i <= n; i += 2) {
		if (sieve[i / 2].fac == 0) {
			sieve[i / 2].fac = i;
			sieve[i / 2].pow = 1;
			if (i <= m) {
				for (j = i * i, k = i / 2; j <= n; j += i + i, k++) {
					if (sieve[j / 2].fac == 0) {
						sieve[j / 2].fac = i;
						if (sieve[k].fac == i) {
							sieve[j / 2].pow = sieve[k].pow + 1;
							sieve[j / 2].nxt = sieve[k].nxt;
						} else {
							sieve[j / 2].pow = 1;
							sieve[j / 2].nxt = k;
						}
					}
				}
			}
		}
	}
}

static void sieve_free(void)
{
	free(sieve);
	sieve = NULL;
}

/* ── Factorized number representation ─────────────────────────────────── */

/* A factorized number: product of fac[i]^pow[i] for i in [0, num_facs).
 * The fac[] array is always sorted in ascending order. */
typedef struct {
	unsigned long max_facs;
	unsigned long num_facs;
	unsigned long *fac;
	unsigned long *pow;
} fac_s;

#define INIT_FACS 32

static fac_s ftmp;  /* scratch for fac_mul_bp */
static fac_s fmul;  /* scratch for fac_mul and fac_remove_gcd */

static void fac_init(fac_s *f)
{
	f->fac = (unsigned long *)malloc(INIT_FACS * sizeof(unsigned long) * 2);
	f->pow = f->fac + INIT_FACS;
	f->max_facs = INIT_FACS;
	f->num_facs = 0;
}

static void fac_clear(fac_s *f)
{
	free(f->fac);
	f->fac = NULL;
	f->pow = NULL;
}

static void fac_reset(fac_s *f)
{
	f->num_facs = 0;
}

static void fac_resize(fac_s *f, long s)
{
	if (f->max_facs < (unsigned long)s) {
		fac_clear(f);
		f->fac = (unsigned long *)malloc(s * sizeof(unsigned long) * 2);
		f->pow = f->fac + s;
		f->max_facs = s;
		f->num_facs = 0;
	}
}

/* f = base^pow, factorized using the sieve */
static void fac_set_bp(fac_s *f, unsigned long base, long pow)
{
	long i;
	assert(base < (unsigned long)sieve_size);
	for (i = 0, base /= 2; base > 0; i++, base = sieve[base].nxt) {
		f->fac[i] = sieve[base].fac;
		f->pow[i] = sieve[base].pow * pow;
	}
	f->num_facs = i;
	assert((unsigned long)i <= f->max_facs);
}

/* r = f * g (merge sorted prime factor lists) */
static void fac_mul2(fac_s *r, fac_s *f, fac_s *g)
{
	unsigned long i, j, k;
	for (i = j = k = 0; i < f->num_facs && j < g->num_facs; k++) {
		if (f->fac[i] == g->fac[j]) {
			r->fac[k] = f->fac[i];
			r->pow[k] = f->pow[i] + g->pow[j];
			i++;
			j++;
		} else if (f->fac[i] < g->fac[j]) {
			r->fac[k] = f->fac[i];
			r->pow[k] = f->pow[i];
			i++;
		} else {
			r->fac[k] = g->fac[j];
			r->pow[k] = g->pow[j];
			j++;
		}
	}
	for (; i < f->num_facs; i++, k++) {
		r->fac[k] = f->fac[i];
		r->pow[k] = f->pow[i];
	}
	for (; j < g->num_facs; j++, k++) {
		r->fac[k] = g->fac[j];
		r->pow[k] = g->pow[j];
	}
	r->num_facs = k;
	assert(k <= r->max_facs);
}

/* f *= g */
static void fac_mul(fac_s *f, fac_s *g)
{
	fac_s tmp;
	fac_resize(&fmul, f->num_facs + g->num_facs);
	fac_mul2(&fmul, f, g);
	tmp = *f;
	*f = fmul;
	fmul = tmp;
}

/* f *= base^pow */
static void fac_mul_bp(fac_s *f, unsigned long base, unsigned long pow)
{
	fac_set_bp(&ftmp, base, pow);
	fac_mul(f, &ftmp);
}

/* Remove factors of power 0 */
static void fac_compact(fac_s *f)
{
	unsigned long i, j;
	for (i = 0, j = 0; i < f->num_facs; i++) {
		if (f->pow[i] > 0) {
			if (j < i) {
				f->fac[j] = f->fac[i];
				f->pow[j] = f->pow[i];
			}
			j++;
		}
	}
	f->num_facs = j;
}

/* Convert factorized form (from fmul scratch) to mpz_t via binary splitting */
static void bs_mul(mpz_t r, long a, long b)
{
	long i, j;
	if (b - a <= 32) {
		mpz_set_ui(r, 1);
		for (i = a; i < b; i++)
			for (j = 0; j < (long)fmul.pow[i]; j++)
				mpz_mul_ui(r, r, fmul.fac[i]);
	} else {
		mpz_t r2;
		mpz_init(r2);
		bs_mul(r2, a, (a + b) / 2);
		bs_mul(r, (a + b) / 2, b);
		mpz_mul(r, r, r2);
		mpz_clear(r2);
	}
}

static mpz_t gcd_val;
static int gcd_initialized = 0;

/* Remove GCD(fp, fg) from both p and g.
 * fp and fg are updated (common powers subtracted).
 * p and g are divided by the computed GCD. */
static void fac_remove_gcd(mpz_ptr p, fac_s *fp, mpz_ptr g, fac_s *fg)
{
	unsigned long i, j, k, c;

	if (!gcd_initialized) {
		mpz_init(gcd_val);
		gcd_initialized = 1;
	}

	fac_resize(&fmul, fp->num_facs < fg->num_facs ? fp->num_facs : fg->num_facs);
	for (i = j = k = 0; i < fp->num_facs && j < fg->num_facs;) {
		if (fp->fac[i] == fg->fac[j]) {
			c = fp->pow[i] < fg->pow[j] ? fp->pow[i] : fg->pow[j];
			fp->pow[i] -= c;
			fg->pow[j] -= c;
			fmul.fac[k] = fp->fac[i];
			fmul.pow[k] = c;
			i++;
			j++;
			k++;
		} else if (fp->fac[i] < fg->fac[j]) {
			i++;
		} else {
			j++;
		}
	}
	fmul.num_facs = k;
	assert(k <= fmul.max_facs);

	if (fmul.num_facs) {
		bs_mul(gcd_val, 0, fmul.num_facs);
		mpz_divexact(p, p, gcd_val);
		mpz_divexact(g, g, gcd_val);
		fac_compact(fp);
		fac_compact(fg);
	}
}

/* Initialize scratch variables. Call once before using fac_mul or fac_remove_gcd. */
static void sieve_helpers_init(void)
{
	fac_init(&ftmp);
	fac_init(&fmul);
}

/* Free scratch variables. */
static void sieve_helpers_free(void)
{
	fac_clear(&ftmp);
	fac_clear(&fmul);
	if (gcd_initialized) {
		mpz_clear(gcd_val);
		gcd_initialized = 0;
	}
}

/* ── Pre-allocated recursion stacks ─────────────────────────────────────── */
/* Array-of-structs layout for cache locality: all data for a level is contiguous. */

typedef struct {
	mpz_t q, t, g;
	fac_s fp, fg;
} bs_level_t;

static bs_level_t *bs_stack;
static long stack_depth;

/* Compute recursion depth for asymmetric split (ratio 0.5224).
 * Max depth = ceil(log(terms) / log(1/0.5224)) + safety margin. */
static long compute_depth(long terms)
{
	if (terms <= 1) return 2;
	long depth = (long)(log((double)terms) / log(1.0/0.5224)) + 4;
	return depth;
}

/* Allocate stacks for binary splitting. Call after sieve_build. */
static void stacks_init(long terms)
{
	long i;
	stack_depth = compute_depth(terms);
	bs_stack = (bs_level_t *)malloc(sizeof(bs_level_t) * stack_depth);
	for (i = 0; i < stack_depth; i++) {
		mpz_init(bs_stack[i].q);
		mpz_init(bs_stack[i].t);
		mpz_init(bs_stack[i].g);
		fac_init(&bs_stack[i].fp);
		fac_init(&bs_stack[i].fg);
	}
}

/* Free stacks. */
static void stacks_free(void)
{
	long i;
	for (i = 0; i < stack_depth; i++) {
		mpz_clear(bs_stack[i].q);
		mpz_clear(bs_stack[i].t);
		mpz_clear(bs_stack[i].g);
		fac_clear(&bs_stack[i].fp);
		fac_clear(&bs_stack[i].fg);
	}
	free(bs_stack);
}

/* ── Binary splitting ──────────────────────────────────────────────────── */

/* Current stack position. Left child reuses top; right child uses top+1. */
static long top = 0;

/* Convenience macros for current and next stack level */
#define q1 (bs_stack[top].q)
#define t1 (bs_stack[top].t)
#define g1 (bs_stack[top].g)
#define fp1 (bs_stack[top].fp)
#define fg1 (bs_stack[top].fg)

#define q2 (bs_stack[top+1].q)
#define t2 (bs_stack[top+1].t)
#define g2 (bs_stack[top+1].g)
#define fp2 (bs_stack[top+1].fp)
#define fg2 (bs_stack[top+1].fg)

/*
 * bs: binary splitting over terms (a, b] (1-indexed).
 * Recursive with explicit top counter for stack indexing.
 * Left child reuses top, right child uses top+1.
 *
 * gflag: 1 = maintain G (needed by caller's merge), 0 = skip
 * level: recursion depth (for GCD threshold)
 *
 * After return, results are in qstack[top]/tstack[top]/gstack[top].
 */
__attribute__((flatten))
static void bs(unsigned long a, unsigned long b, int gflag, long level)
{
	unsigned long mid;

	if (b - a == 1) {
		/*
		 * Base case: single Chudnovsky term for 1-indexed term b.
		 *   Q = b^3 * C^3/24
		 *   G = (6b-5)(2b-1)(6b-1)
		 *   T = G * (A + B*b) * (-1)^b
		 */
		unsigned long i;

		mpz_set_ui(q1, b);
		mpz_mul_ui(q1, q1, b);
		mpz_mul_ui(q1, q1, b);
		mpz_mul_ui(q1, q1, (640320UL/24)*(640320UL/24));
		mpz_mul_ui(q1, q1, 640320UL*24);

		mpz_set_ui(g1, 2*b-1);
		mpz_mul_ui(g1, g1, 6*b-1);
		mpz_mul_ui(g1, g1, 6*b-5);

		mpz_set_ui(t1, b);
		mpz_mul_ui(t1, t1, 545140134UL);
		mpz_add_ui(t1, t1, 13591409UL);
		mpz_mul(t1, t1, g1);
		if (b % 2)
			mpz_neg(t1, t1);

		i = b;
		while ((i & 1) == 0) i >>= 1;
		fac_set_bp(&fp1, i, 3);
		fac_mul_bp(&fp1, 3*5*23*29, 3);
		fp1.pow[0]--;

		fac_set_bp(&fg1, 2*b-1, 1);
		fac_mul_bp(&fg1, 6*b-1, 1);
		fac_mul_bp(&fg1, 6*b-5, 1);
		return;
	}

	/* Asymmetric split (tuning parameter from gmp-chudnovsky) */
	mid = a + (unsigned long)((b - a) * 0.5224);
	if (mid == a) mid = a + 1;
	if (mid >= b) mid = b - 1;

	/* Left half: always maintain G (needed for merge) */
	bs(a, mid, 1, level + 1);

	/* Right half: uses top+1 */
	top++;
	bs(mid, b, gflag, level + 1);
	top--;

	/* GCD removal at depth >= 4 */
	if (level >= 4) {
		fac_remove_gcd(q2, &fp2, g1, &fg1);
	}

	/* Merge:
	 *   T = T_L * Q_R + G_L * T_R
	 *   Q = Q_L * Q_R
	 *   G = G_L * G_R (if gflag)
	 */
	mpz_mul(t1, t1, q2);
	mpz_mul(t2, t2, g1);
	mpz_add(t1, t1, t2);
	mpz_mul(q1, q1, q2);

	fac_mul(&fp1, &fp2);

	if (gflag) {
		mpz_mul(g1, g1, g2);
		fac_mul(&fg1, &fg2);
	}
}

/*
 * binary_split: entry point.
 * Splits terms 1..N. Results in bs_stack[0].q (Q) and bs_stack[0].t (T).
 * The k=0 term (constant A) is added in pi_final.
 */
static void binary_split(long N)
{
	top = 0;
	bs(0, (unsigned long)N, 0, 0);
}

/* ── Newton's method sqrt (precision doubling) ─────────────────────────── */
/* Adapted from gmp-chudnovsky.c by Hanhong Xue.
 * Computes r = sqrt(x) using Newton iteration on 1/sqrt(x),
 * doubling precision at each step. Much faster than mpf_sqrt
 * for large precisions. */

#define DOUBLE_PREC 53

static mpf_t nt1, nt2;
static int newton_initialized = 0;

static void newton_init(unsigned long prec)
{
	mpf_init2(nt1, prec);
	mpf_init2(nt2, prec);
	newton_initialized = 1;
}

static void newton_free(void)
{
	if (newton_initialized) {
		mpf_clear(nt1);
		mpf_clear(nt2);
		newton_initialized = 0;
	}
}

static void my_sqrt_ui(mpf_t r, unsigned long x)
{
	unsigned long prec, bits, prec0;

	prec0 = mpf_get_prec(r);

	if (prec0 <= DOUBLE_PREC) {
		mpf_set_d(r, sqrt(x));
		return;
	}

	bits = 0;
	for (prec = prec0; prec > DOUBLE_PREC;) {
		int bit = prec & 1;
		prec = (prec + bit) / 2;
		bits = bits * 2 + bit;
	}

	mpf_set_prec_raw(nt1, DOUBLE_PREC);
	mpf_set_d(nt1, 1.0 / sqrt((double)x));

	while (prec < prec0) {
		prec *= 2;
		if (prec < prec0) {
			/* nt1 = nt1 + nt1*(1 - x*nt1*nt1)/2 */
			mpf_set_prec_raw(nt2, prec);
			mpf_mul(nt2, nt1, nt1);
			mpf_mul_ui(nt2, nt2, x);
			mpf_ui_sub(nt2, 1, nt2);
			mpf_set_prec_raw(nt2, prec / 2);
			mpf_div_2exp(nt2, nt2, 1);
			mpf_mul(nt2, nt2, nt1);
			mpf_set_prec_raw(nt1, prec);
			mpf_add(nt1, nt1, nt2);
		} else {
			break;
		}
		prec -= (bits & 1);
		bits /= 2;
	}
	/* nt2 = x*nt1, r = nt2 + nt1*(x - nt2*nt2)/2 */
	mpf_set_prec_raw(nt2, prec0 / 2);
	mpf_mul_ui(nt2, nt1, x);
	mpf_mul(r, nt2, nt2);
	mpf_ui_sub(r, x, r);
	mpf_mul(nt1, nt1, r);
	mpf_div_2exp(nt1, nt1, 1);
	mpf_add(r, nt1, nt2);
}

/* ── Final pi computation ──────────────────────────────────────────────── */
/* Combines: Q * (C/D) * sqrt(C) / T
 * where C=640320, D=12, C/D=53360.
 * Does Q *= C/D in integer space to save a float multiply.
 * Uses Newton sqrt for speed. Called from LLVM IR after binary_split. */

#define C 640320
#define D 12

long pi_final(long digits, char *buf, long buf_len)
{
	/* log2(10) ≈ 3.32192809489; exact precision needed for decimal digits */
	unsigned long prec = (unsigned long)(digits * 3.32192809489) + 256;
	mpf_t pi_f, q_f, sqrt_f;

	newton_init(prec);

	/* Add k=0 term: T += A*Q (integer addmul, since k=0 not in binary split) */
	mpz_addmul_ui(bs_stack[0].t, bs_stack[0].q, 13591409UL);

	/* Q *= C/D in integer space (cheaper than float mul_ui) */
	mpz_mul_ui(bs_stack[0].q, bs_stack[0].q, C/D);

	mpf_init2(pi_f, prec);
	mpf_init2(q_f, prec);
	mpf_init2(sqrt_f, prec);

	mpf_set_z(q_f, bs_stack[0].q);

	/* Division: pi_f = Q*(C/D) / T */
	{
		mpf_t t_f;
		mpf_init2(t_f, prec);
		mpf_set_z(t_f, bs_stack[0].t);
		mpf_div(pi_f, q_f, t_f);
		mpf_clear(t_f);
	}

	/* sqrt(C) via Newton's method */
	my_sqrt_ui(sqrt_f, C);

	/* pi = Q*(C/D)/T * sqrt(C) — single float multiply */
	mpf_mul(pi_f, pi_f, sqrt_f);

	/* Format to string using mpf_get_str */
	{
		mp_exp_t exp;
		char *str = mpf_get_str(NULL, &exp, 10, digits + 2, pi_f);
		buf[0] = str[0];
		buf[1] = '.';
		long copy_digits = digits;
		long str_len = (long)strlen(str);
		if (copy_digits > str_len - 1) copy_digits = str_len - 1;
		memcpy(buf + 2, str + 1, copy_digits);
		buf[copy_digits + 2] = '\0';
		free(str);
	}
	long copy_len = digits + 2;
	mpf_clear(pi_f);
	mpf_clear(q_f);
	mpf_clear(sqrt_f);
	newton_free();

	return copy_len;
}

/* ── Full pi computation ───────────────────────────────────────────────── */
/* Called from LLVM IR compute_pi. Handles everything: sieve, split, final. */

#define DIGITS_PER_ITER 14.1816474627254776555

long compute_pi_c(long digits, char *buf, long buf_len)
{
	long N = (long)(digits / DIGITS_PER_ITER) + 1;
	long sieve_sz = N * 6;
	if (sieve_sz < 10006) sieve_sz = 10006;

	sieve_build(sieve_sz);
	sieve_helpers_init();
	stacks_init(N);

	binary_split(N);

	sieve_helpers_free();
	sieve_free();

	long result = pi_final(digits, buf, buf_len);

	stacks_free();
	return result;
}
