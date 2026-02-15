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
