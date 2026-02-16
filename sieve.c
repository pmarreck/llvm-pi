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
#include <time.h>
#include <gmp.h>
#include <pthread.h>

#ifdef PHASE_TIMING
#include <stdio.h>
static double phase_now(void) {
	struct timespec ts;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
	return ts.tv_sec + ts.tv_nsec * 1e-9;
}
#endif

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

static __thread fac_s ftmp;  /* scratch for fac_mul_bp (thread-local for parallel bs) */
static __thread fac_s fmul;  /* scratch for fac_mul and fac_remove_gcd (thread-local) */

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
void fac_set_bp(fac_s *f, unsigned long base, long pow)
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
void fac_mul(fac_s *f, fac_s *g)
{
	fac_s tmp;
	fac_resize(&fmul, f->num_facs + g->num_facs);
	fac_mul2(&fmul, f, g);
	tmp = *f;
	*f = fmul;
	fmul = tmp;
}

/* f *= base^pow */
void fac_mul_bp(fac_s *f, unsigned long base, unsigned long pow)
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

static __thread mpz_t gcd_val;
static __thread int gcd_initialized = 0;

/* Remove GCD(fp, fg) from both p and g.
 * fp and fg are updated (common powers subtracted).
 * p and g are divided by the computed GCD. */
void fac_remove_gcd(mpz_ptr p, fac_s *fp, mpz_ptr g, fac_s *fg)
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

__thread bs_level_t *bs_stack;
static __thread long stack_depth;

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
		/* Pre-allocate 128 bits (2 limbs) so hand-written LLVM IR can
		   write directly to the limb data without triggering GMP realloc.
		   Base case Q/G/T all fit in 2 limbs; GMP will grow as needed
		   for larger merge results at higher recursion levels. */
		mpz_init2(bs_stack[i].q, 128);
		mpz_init2(bs_stack[i].t, 128);
		mpz_init2(bs_stack[i].g, 128);
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

/* Current stack position (non-static: accessed from hand-written LLVM IR).
 * Thread-local for parallel binary splitting. */
__thread long top = 0;

/* bs() and binary_split() are implemented in hand-tuned LLVM IR (chudnovsky.ll).
 * The IR loads @bs_stack and @top once per call, eliminating ~30 redundant
 * global reloads per base case that the C compiler conservatively inserts
 * after every GMP function call. */
extern void binary_split(long N);
extern void bs(long a, long b, int gflag, long level);

/* ── Parallel binary splitting ─────────────────────────────────────────── */
/* Split work into chunks, each computed by a separate thread.
 * Controlled by PI_THREADS environment variable (default: 1).
 * Each thread gets its own TLS bs_stack/top/fac_s scratch. */

#include <stdio.h>

typedef struct {
	long a, b;        /* half-open interval of terms */
	long terms;       /* b - a, for stack depth */
	mpz_t q, t, g;   /* results (pre-initialized by main thread) */
} chunk_t;

static void *bs_worker(void *arg)
{
	chunk_t *chunk = (chunk_t *)arg;
	long depth = compute_depth(chunk->terms);
	long i;

	/* Initialize thread-local state */
	bs_stack = (bs_level_t *)malloc(sizeof(bs_level_t) * depth);
	stack_depth = depth;
	for (i = 0; i < depth; i++) {
		mpz_init2(bs_stack[i].q, 128);
		mpz_init2(bs_stack[i].t, 128);
		mpz_init2(bs_stack[i].g, 128);
		fac_init(&bs_stack[i].fp);
		fac_init(&bs_stack[i].fg);
	}
	fac_init(&ftmp);
	fac_init(&fmul);
	gcd_initialized = 0;
	top = 0;

	/* Compute this chunk (gflag=1: merge needs G) */
	bs(chunk->a, chunk->b, 1, 0);

	/* Transfer results via swap (O(1), no limb copying) */
	mpz_swap(chunk->q, bs_stack[0].q);
	mpz_swap(chunk->t, bs_stack[0].t);
	mpz_swap(chunk->g, bs_stack[0].g);

	/* Cleanup thread-local state */
	for (i = 0; i < depth; i++) {
		mpz_clear(bs_stack[i].q);
		mpz_clear(bs_stack[i].t);
		mpz_clear(bs_stack[i].g);
		fac_clear(&bs_stack[i].fp);
		fac_clear(&bs_stack[i].fg);
	}
	free(bs_stack);
	bs_stack = NULL;
	fac_clear(&ftmp);
	fac_clear(&fmul);
	if (gcd_initialized) {
		mpz_clear(gcd_val);
		gcd_initialized = 0;
	}

	return NULL;
}

static void parallel_binary_split(long N, int n_threads)
{
	chunk_t *chunks;
	pthread_t *threads;
	int i;

	chunks = (chunk_t *)malloc(n_threads * sizeof(chunk_t));
	threads = (pthread_t *)malloc((n_threads - 1) * sizeof(pthread_t));

	/* Divide terms into chunks */
	for (i = 0; i < n_threads; i++) {
		chunks[i].a = (long)i * N / n_threads;
		chunks[i].b = (long)(i + 1) * N / n_threads;
		chunks[i].terms = chunks[i].b - chunks[i].a;
		mpz_init(chunks[i].q);
		mpz_init(chunks[i].t);
		mpz_init(chunks[i].g);
	}

	/* Spawn worker threads for chunks 1..n_threads-1 */
	for (i = 1; i < n_threads; i++) {
		pthread_create(&threads[i - 1], NULL, bs_worker, &chunks[i]);
	}

	/* Main thread does chunk 0 using its own TLS state */
	top = 0;
	bs(chunks[0].a, chunks[0].b, 1, 0);
	mpz_swap(chunks[0].q, bs_stack[0].q);
	mpz_swap(chunks[0].t, bs_stack[0].t);
	mpz_swap(chunks[0].g, bs_stack[0].g);

	/* Wait for all workers */
	for (i = 1; i < n_threads; i++) {
		pthread_join(threads[i - 1], NULL);
	}

	/* Merge results sequentially: left-to-right accumulation
	 * T = T_left * Q_right + G_left * T_right
	 * Q = Q_left * Q_right
	 * G = G_left * G_right (only if needed for next merge) */
	for (i = 1; i < n_threads; i++) {
		mpz_mul(chunks[0].t, chunks[0].t, chunks[i].q);
		mpz_addmul(chunks[0].t, chunks[i].t, chunks[0].g);
		mpz_mul(chunks[0].q, chunks[0].q, chunks[i].q);
		if (i < n_threads - 1) {
			mpz_mul(chunks[0].g, chunks[0].g, chunks[i].g);
		}
		mpz_clear(chunks[i].q);
		mpz_clear(chunks[i].t);
		mpz_clear(chunks[i].g);
	}

	/* Put merged result back in main thread's bs_stack[0] */
	mpz_swap(bs_stack[0].q, chunks[0].q);
	mpz_swap(bs_stack[0].t, chunks[0].t);
	mpz_clear(chunks[0].q);
	mpz_clear(chunks[0].t);
	mpz_clear(chunks[0].g);

	free(threads);
	free(chunks);
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

/* ── Binary-splitting base conversion ────────────────────────────────────── */
/* Replaces GMP's mpf_get_str with a parallelizable divide-and-conquer
 * algorithm. GMP's internal conversion is already subquadratic, but it
 * cannot be parallelized. Our implementation enables multi-threaded
 * formatting via the same PI_THREADS env var that controls binary splitting.
 *
 * Algorithm: precompute 10^(2^k) table, then recursively divide-and-conquer:
 *   to_decimal(Z, n) = to_decimal(Z / 10^half, n-half) || to_decimal(Z % 10^half, half)
 * Parallelism: after each split, spawn a thread for the larger half. */

static mpz_t *fmt_pow10;
static int fmt_depth;

static void fmt_pow10_build(long digits)
{
	int depth = 0;
	long d = 1;
	while (d < digits) { d <<= 1; depth++; }
	fmt_depth = depth;
	fmt_pow10 = (mpz_t *)malloc((depth + 1) * sizeof(mpz_t));
	mpz_init_set_ui(fmt_pow10[0], 10);
	for (int k = 1; k <= depth; k++) {
		mpz_init(fmt_pow10[k]);
		mpz_mul(fmt_pow10[k], fmt_pow10[k - 1], fmt_pow10[k - 1]);
	}
}

static void fmt_pow10_free(void)
{
	for (int k = 0; k <= fmt_depth; k++)
		mpz_clear(fmt_pow10[k]);
	free(fmt_pow10);
	fmt_pow10 = NULL;
}

/* Compute 10^n from precomputed table by binary decomposition of n */
static void fmt_pow10_exp(mpz_t result, long n)
{
	mpz_set_ui(result, 1);
	int k = 0;
	while (n > 0) {
		if (n & 1) mpz_mul(result, result, fmt_pow10[k]);
		n >>= 1;
		k++;
	}
}

#define FMT_BASE_CASE 1024

/* Sequential binary-splitting conversion.
 * Writes exactly ndigits decimal characters to buf, left-padded with '0'. */
static void fmt_convert(char *buf, long ndigits, mpz_t z, int level)
{
	if (ndigits <= FMT_BASE_CASE) {
		char *s = mpz_get_str(NULL, 10, z);
		long len = (long)strlen(s);
		long pad = ndigits - len;
		if (pad > 0) memset(buf, '0', pad);
		memcpy(buf + (pad > 0 ? pad : 0), s, len);
		free(s);
		return;
	}

	/* Find level where 2^(level-1) < ndigits */
	while (level > 0 && (1L << (level - 1)) >= ndigits)
		level--;
	if (level <= 0) {
		char *s = mpz_get_str(NULL, 10, z);
		long len = (long)strlen(s);
		long pad = ndigits - len;
		if (pad > 0) memset(buf, '0', pad);
		memcpy(buf + (pad > 0 ? pad : 0), s, len);
		free(s);
		return;
	}

	long half = 1L << (level - 1);
	mpz_t q, r;
	mpz_init(q);
	mpz_init(r);
	mpz_tdiv_qr(q, r, z, fmt_pow10[level - 1]);

	fmt_convert(buf, ndigits - half, q, level - 1);
	fmt_convert(buf + (ndigits - half), half, r, level - 1);

	mpz_clear(q);
	mpz_clear(r);
}

/* Tree-parallel conversion: after each split, spawn a thread for the right
 * (larger) half and recurse on the left with remaining threads. */
typedef struct {
	char *buf;
	long ndigits;
	mpz_t z;
	int level;
	int n_threads;
} fmt_par_t;

static void fmt_par_convert(char *buf, long ndigits, mpz_t z,
                            int level, int n_threads);

static void *fmt_par_worker(void *arg)
{
	fmt_par_t *p = (fmt_par_t *)arg;
	fmt_par_convert(p->buf, p->ndigits, p->z, p->level, p->n_threads);
	return NULL;
}

static void fmt_par_convert(char *buf, long ndigits, mpz_t z,
                            int level, int n_threads)
{
	if (n_threads <= 1 || ndigits <= FMT_BASE_CASE * 2) {
		fmt_convert(buf, ndigits, z, level);
		return;
	}

	/* Find split level */
	while (level > 0 && (1L << (level - 1)) >= ndigits)
		level--;
	if (level <= 0) {
		fmt_convert(buf, ndigits, z, level);
		return;
	}

	long half = 1L << (level - 1);
	mpz_t q, r;
	mpz_init(q);
	mpz_init(r);
	mpz_tdiv_qr(q, r, z, fmt_pow10[level - 1]);

	/* Allocate threads proportional to digit count */
	long left_d = ndigits - half;
	int right_t = (int)(0.5 + (double)half / ndigits * n_threads);
	if (right_t < 1) right_t = 1;
	if (right_t >= n_threads) right_t = n_threads - 1;
	int left_t = n_threads - right_t;

	/* Spawn thread for right (larger) half */
	fmt_par_t right_arg;
	right_arg.buf = buf + left_d;
	right_arg.ndigits = half;
	mpz_init_set(right_arg.z, r);
	right_arg.level = level - 1;
	right_arg.n_threads = right_t;

	pthread_t thread;
	pthread_create(&thread, NULL, fmt_par_worker, &right_arg);

	/* Left half on current thread */
	fmt_par_convert(buf, left_d, q, level - 1, left_t);

	pthread_join(thread, NULL);
	mpz_clear(right_arg.z);
	mpz_clear(q);
	mpz_clear(r);
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

#ifdef PHASE_TIMING
	double pt0, pt1;
	pt0 = phase_now();
#endif

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

#ifdef PHASE_TIMING
	pt1 = phase_now();
	fprintf(stderr, "  div   time = %.3f\n", pt1 - pt0);
	pt0 = pt1;
#endif

	/* sqrt(C) via Newton's method */
	my_sqrt_ui(sqrt_f, C);

#ifdef PHASE_TIMING
	pt1 = phase_now();
	fprintf(stderr, "  sqrt  time = %.3f\n", pt1 - pt0);
	pt0 = pt1;
#endif

	/* pi = Q*(C/D)/T * sqrt(C) — single float multiply */
	mpf_mul(pi_f, pi_f, sqrt_f);

#ifdef PHASE_TIMING
	pt1 = phase_now();
	fprintf(stderr, "  mul   time = %.3f\n", pt1 - pt0);
	pt0 = pt1;
#endif

	/* Format pi to decimal string */
	{
		int fmt_threads = 1;
		{
			const char *env = getenv("PI_THREADS");
			if (env) fmt_threads = atoi(env);
			if (fmt_threads < 1) fmt_threads = 1;
		}

		if (fmt_threads > 1) {
			/* Parallel binary-splitting base conversion:
			 * Scale pi to integer, then divide-and-conquer with threads. */
			long guard = 5;
			long total_d = digits + guard;

			fmt_pow10_build(total_d + 1);

#ifdef PHASE_TIMING
			pt1 = phase_now();
			fprintf(stderr, "  p10   time = %.3f\n", pt1 - pt0);
			pt0 = pt1;
#endif

			/* Scale pi to integer: pi_z = floor(pi_f * 10^total_d) */
			{
				mpz_t ten_d;
				mpz_init(ten_d);
				fmt_pow10_exp(ten_d, total_d);

				mpf_t sf;
				mpf_init2(sf, prec + 64);
				mpf_set_z(sf, ten_d);
				mpf_mul(pi_f, pi_f, sf);
				mpf_clear(sf);
				mpz_clear(ten_d);
			}

			mpz_t pi_z;
			mpz_init(pi_z);
			mpz_set_f(pi_z, pi_f);

#ifdef PHASE_TIMING
			pt1 = phase_now();
			fprintf(stderr, "  scale time = %.3f\n", pt1 - pt0);
			pt0 = pt1;
#endif

			char *raw = (char *)malloc(total_d + 10);
			fmt_par_convert(raw, total_d + 1, pi_z, fmt_depth, fmt_threads);
			mpz_clear(pi_z);

			buf[0] = raw[0];
			buf[1] = '.';
			memcpy(buf + 2, raw + 1, digits);
			buf[digits + 2] = '\0';
			free(raw);

			fmt_pow10_free();
		} else {
			/* Single-threaded: use GMP's optimized mpf_get_str directly */
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
	}

#ifdef PHASE_TIMING
	pt1 = phase_now();
	fprintf(stderr, "  fmt   time = %.3f\n", pt1 - pt0);
#endif
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

#ifdef PHASE_TIMING
	double t0, t1;
	t0 = phase_now();
#endif

	sieve_build(sieve_sz);
	sieve_helpers_init();
	stacks_init(N);

#ifdef PHASE_TIMING
	t1 = phase_now();
	fprintf(stderr, "sieve   time = %.3f\n", t1 - t0);
	t0 = t1;
#endif

	{
		int n_threads = 1;
		const char *env = getenv("PI_THREADS");
		if (env) n_threads = atoi(env);
		if (n_threads < 1) n_threads = 1;

		if (n_threads > 1) {
			parallel_binary_split(N, n_threads);
		} else {
			binary_split(N);
		}
	}

#ifdef PHASE_TIMING
	t1 = phase_now();
	fprintf(stderr, "bs      time = %.3f\n", t1 - t0);
	t0 = t1;
#endif

	sieve_helpers_free();
	sieve_free();

	long result = pi_final(digits, buf, buf_len);

#ifdef PHASE_TIMING
	t1 = phase_now();
	fprintf(stderr, "final   time = %.3f\n", t1 - t0);
#endif

	stacks_free();
	return result;
}
