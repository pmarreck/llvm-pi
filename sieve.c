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
void sieve_build(long n)
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

void sieve_free(void)
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

void fac_init(fac_s *f)
{
	f->fac = (unsigned long *)malloc(INIT_FACS * sizeof(unsigned long) * 2);
	f->pow = f->fac + INIT_FACS;
	f->max_facs = INIT_FACS;
	f->num_facs = 0;
}

void fac_clear(fac_s *f)
{
	free(f->fac);
	f->fac = NULL;
	f->pow = NULL;
}

void fac_reset(fac_s *f)
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

static mpz_t gcd_val;
static int gcd_initialized = 0;

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
void sieve_helpers_init(void)
{
	fac_init(&ftmp);
	fac_init(&fmul);
}

/* Free scratch variables. */
void sieve_helpers_free(void)
{
	fac_clear(&ftmp);
	fac_clear(&fmul);
	if (gcd_initialized) {
		mpz_clear(gcd_val);
		gcd_initialized = 0;
	}
}
