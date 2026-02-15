/* Pi computation using MPFR's mpfr_const_pi (Gauss-Legendre / AGM algorithm).
 *
 * This is a fundamentally different algorithm from Chudnovsky — it uses
 * the arithmetic-geometric mean (AGM) which converges quadratically
 * (doubling correct digits each iteration).
 *
 * Usage: mpfr-pi <digits> [output]
 *   digits: number of decimal digits to compute
 *   output: 1 to print digits (default: 0, timing only)
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>

int main(int argc, char *argv[])
{
	long digits = 1000;
	int output = 0;

	if (argc > 1)
		digits = atol(argv[1]);
	if (argc > 2)
		output = atoi(argv[2]);

	/* bits needed: digits * log2(10) + margin */
	mpfr_prec_t bits = (mpfr_prec_t)(digits * 3.32192809488736235) + 64;

	mpfr_t pi;
	mpfr_init2(pi, bits);
	mpfr_const_pi(pi, MPFR_RNDN);

	if (output) {
		/* Print "3." followed by exactly 'digits' decimal digits */
		mpfr_printf("%.*Rf\n", digits, pi);
	}

	mpfr_clear(pi);
	mpfr_free_cache();
	return 0;
}
