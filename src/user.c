#include "../headers/ecm.h"

/* User function : find a non-trivial big factor of N */
void factor(mpz_t d, mpz_t N, unsigned long B, unsigned long *primes, unsigned long *differences)
{	
	if (mpz_cmp_ui(N, 0) == 0)
	{
		printf("Input is 0.\n");
		return;
	}
	gmp_printf("Factoring %Zd.\n", N);

	trial_division(N, primes);

	if ((mpz_cmp_ui(N,1) == 0)|| mpz_probab_prime_p(N,25) > 0)
	{	
		printf("Factorized by trial division only.\n");
		return;
	}

	gmp_printf("Calling ECM on : %Zd.\n", N);

	clock_t st = clock();

	ECM_factor(d, N, B, primes, differences);

	gmp_printf("\nFactor found : %Zd.\n", d);
	printf("\nElapsed time for ECM: %f\n", (double)(clock() - st) / CLOCKS_PER_SEC);
	return;
}