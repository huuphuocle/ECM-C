/* C file containing auxiliary functions : 
- trial division
- Miller-Rabin primality test
*/

#include "../headers/ecm.h"

/* Trial division for N by a list of precomputed primes */

void trial_division(mpz_t N, unsigned long *primes)
{	
	unsigned long l = primes[0];
	int e;
	mpz_t tmp;
	mpz_init(tmp);
	for (int i = 1; i <= l; i++)
	{
		mpz_mod_ui(tmp, N, primes[i]);
		e = 0;
		while (mpz_cmp_ui(tmp, 0) == 0)
		{
			e++;
			mpz_divexact_ui(N, N, primes[i]);
			mpz_mod_ui(tmp, N, primes[i]);
		}
	}
	mpz_clear(tmp);
	return;
}

/* Miller-Rabin test : return 1 if N is composite, 0 if N is probably prime */

int is_composite(mpz_t N, int c)
{

	mpz_t q, a, b, t, e, tmp;
	mpz_inits(q, a, b, t, e, tmp, NULL);

	mpz_sub_ui(q, N, 1); /* q = N-1 */
	mpz_set_ui(t, 0);

	while (mpz_divisible_2exp_p(q, 1))
	{
		mpz_fdiv_q_2exp(q, q, 1);
		mpz_add_ui(t, t, 1);
	}
	mpz_sub_ui(t, t, 1);

	while (c > 0)
	{
		mpz_sub_ui(tmp, N, 3);
		mpz_set_ui(a, rand());
		mpz_mod(a, a, tmp);

		mpz_add_ui(a, a, 2); // a = a+2
		mpz_sub_ui(tmp, N, 1);

		mpz_set_ui(e, 0);
		mpz_powm(b, a, q, N);
		if (mpz_cmp_ui(b, 1) != 0)
		{
			while ((mpz_cmp_ui(b, 1) != 0) && (mpz_cmp(b, tmp) != 0) && (mpz_cmp(e, t) < 0))
			{
				mpz_powm_ui(b, b, 2, N);
				mpz_add_ui(e, e, 1);
			}
			if (mpz_cmp(b, tmp) != 0)
			{
				mpz_clears(q, a, b, t, e, tmp, NULL);
				return 1;
			}
		}
		c = c - 1;
	}
	mpz_clears(q, a, b, t, e, tmp, NULL);
	return 0;
}

/* Precompute 
- an array T of primes p up to B2
- an array D of the difference / 2 of primes in [B1,B2]
- T[0] is the number of primes up to B2
- D[0] is the number of primes in [B1,B2] */

void precompute(unsigned long *T, unsigned long B1, unsigned long *D, unsigned long B2)
{
	unsigned long n = 2, c = 0, j = 0, i;
	clock_t st = clock();
	int flag = 1;
	while (n < B1)
	{
		flag = 1;
		for (unsigned long i = 1; i <= c; i++)
		{
			if (T[i] > sqrt(n))
				break;
			if ((n % T[i]) == 0)
			{
				flag = 0;
				break;
			}
		}
		if (flag == 1)
		{
			c++;
			T[c] = n;
		}
		n++;
	}
	T[0] = c;
	unsigned long max_d = 0;
	while (n < B2)
	{
		flag = 1;
		for (i = 1; i <= c; i++)
		{
			if (T[i] > sqrt(n))
				break;
			if ((n % T[i]) == 0)
			{
				flag = 0;
				break;
			}
		}
		if (flag == 1)
		{
			c++;
			j++;
			T[c] = n;
			D[j] = (T[c] - T[c - 1]) >> 1;
			if (D[j] > max_d) max_d = D[j];
		}
		n++;
	}
	D[0] = j;
	D[j+1] = max_d;
	
	printf("Precompute : %lu primes B1 ; %lu primes B2 ; max_difference = %lu.\nElapsed time: %f (s) \n", T[0], D[0], D[D[0]+1], (double)(clock() - st) / CLOCKS_PER_SEC);
	return;
}
