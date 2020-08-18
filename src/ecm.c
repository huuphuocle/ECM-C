#include "../headers/ecm.h"

void Suyama_param(mpz_t a24, mpz_t A, mpz_t X, mpz_t Z, unsigned int m){
	mpz_t u, v, tmp;
	mpz_inits(u, v, tmp, NULL);

	mpz_set_ui(u, m*m-5); // u = m^2-5
	mpz_mul_2exp(v, m << 2); // v = 4m
	mpz_powm_ui(X, u, 3, N);
	mpz_powm_ui(Z, v, 3, N);
	
	mpz_mul(tmp,X,v);
	mpz_mul_2exp(tmp, tmp, 4);
	mpz_invert(tmp, tmp, N);
	mpz_sub(a24, v, u);
	mpz_powm_ui(a24, a24, 3, N);
	mpz_mul(a24, a24, tmp);
	mpz_mul_ui(tmp, u, 3);
	mpz_add(tmp, tmp, v);
	mpz_mul(a24, a24, tmp);
	mpz_mod(a24, a24, N);
	
	mpz_mul_2exp(A, a24, 2);
	mpz_add_ui(A, A, 2);
	mpz_mod(A, A, N);
	mpz_clears(tmp, NULL);
	return;
}

void normal_param(mpz_t a24, mpz_t A, mpz_t X, mpz_t Z, unsigned int m, mpz_t mpz_4){
	mpz_t m2, k2, k;
	mpz_inits(m2, k2, k, NULL);
	
	mpz_set_ui(m2, m); 
	mpz_powm_ui(m2, m2, 2, N); // m2 = m^2

	mpz_sub_ui(k2, m2, 1); 
	mpz_mul_2exp(k2, k2, 1);
	mpz_invert(k2, k2, N); // k2 = 1/2*(m^2-1);

	
	mpz_sub_ui(k ,m2,4);
	mpz_neg(k, k);
	mpz_mod(k, k, N);
	mpz_mul(k, k, k2); // k = (4 - m^2)/2*(m^2 - 1)
	
	mpz_invert(A, k, N);
	mpz_add(A, A, k);
	mpz_mod(A, A, N); // A = k + 1/k

	mpz_add_ui(a24, A, 2);
	mpz_mul(a24, a24, mpz_4);
	mpz_mod(a24, a24, N);

	mpz_set_ui(X, 2);
	mpz_set_ui(Z, 1);

	mpz_clears(m2, k2, k, NULL);
	return;
}

/* Split N, using the bound B, e = lg(B), T is the table of primes up to B */
void ECM_factor(mpz_t d, mpz_t N, unsigned long B1, unsigned long * primes, unsigned long * diff)
{	
	int flag;
	unsigned int i, j, m = 2;
	unsigned long q, q1, q2;

	mpz_t X, Y, Z, x, y, z, x1, y1, z1, m2, k2, k, A, B, a24, mpz_4, Zt, X3, Z2, Z3, tmp, a;
	mpz_inits(X, Y, Z, x, y, z, x1, y1, z1, m2, k2, k, A, B, a24, mpz_4, Zt, X3, Z2, Z3, tmp, a, NULL);
	

	unsigned int l1 = primes[0], l2 = diff[0]; // number of primes for stage 1

	// largest difference in the array diff
	unsigned long maxd = diff[l2+1];
	
    /* allocate memory for the table of powers of Q in Stage 2 */
	mpz_t * powX = malloc(maxd * sizeof(mpz_t));
	mpz_t * powY = malloc(maxd * sizeof(mpz_t));
	mpz_t * powZ = malloc(maxd * sizeof(mpz_t));

	for(i=0; i < maxd; i++){
		mpz_init(powX[i]);
		mpz_init(powY[i]);
		mpz_init(powZ[i]);
	}
	mpz_set_ui(mpz_4, 4);
	mpz_invert(mpz_4, mpz_4, N);
	/*---------------------------------------------*/
    
	/* choose the curve and base point using Suyama's parametrization */
loop_curve:
	
	m++;
	printf("[%u]",m);
	
	normal_param(a24, A, X, Z, m, mpz_4);

	flag = 0;

	goto step1;

loop_point:
	if(flag == 1){
		mpz_set_ui(X, 2);
		// mpz_add_ui(Y, m2, 2);
		// mpz_invert(Y, Y, N);
		// mpz_mul_ui(Y, Y, 6*m);
		// mpz_mod(Y, Y, N);
		mpz_set_ui(Z, 1);
	}else{
		goto loop_curve;
	}
	
	/* Montgomery's curve */

	/* Stage 1 */

step1:
	/* for loop computes M = prod(p_i^l_i) */
	gmp_printf("%Zd %Zd \n",X,Z);
	for (i = 1; i <= l1; i++){
		printf("i = %d \n", i);
		q = primes[i];
		q2 = q;
		while ( q2 <= B1 ){
			q1 = q2;
			q2 = q2 * q;
		}
		
		ladder_ui(X, Z, q1, X, Z, N, a24);
		gmp_printf("%Zd %Zd \n",X,Z);
		if(mpz_cmp_ui(Z,0) == 0){
			flag = flag ^ 1;
			goto loop_point;
		}else{
			mpz_gcd(d, Z, N);
			if (mpz_cmp_ui(d, 1)){
				gmp_printf("d = %Zd \n", d);
				goto final_step;
			}
			goto loop_curve;
		}
	}
	
	// (X, Z) = Q = [M](P) 

// step2:	
// 	printf("there\n");
// 	/* Stage 2 - List of primes : p_1, ... , p_s */	

// 	// convertWM(x,y,z,X,Z); // convert (X,Z) to Weierstrass point

// 	mpz_set_ui(Y, 1); // set Y = 1; Q = (X, 1, Z)

// 	// compute B
// 	mpz_invert(Zt, Z, N);
// 	mpz_powm_ui(X3, X, 3, N);
// 	mpz_powm_ui(Z2, Z, 2, N);
// 	mpz_mul(Z3, Z2, Z);
// 	mpz_mul(B, X, Z2);
// 	mpz_mul(B, B, A);
// 	mpz_add(B, B, Z3);
// 	mpz_add(B, B, X3);
// 	mpz_mul(B, B, Zt);
// 	mpz_mod(B, B, N);
	
// 	// compute a = (3-A^2)/(3B^2)
// 	mpz_powm_ui(tmp, B, 2, N);
// 	mpz_mul_ui(tmp, tmp, 3);
// 	mpz_invert(tmp, tmp, N);
// 	mpz_powm_ui(a, A, 2, N);
// 	mpz_sub_ui(a, a, 3);
// 	mpz_neg(a, a);
// 	mpz_mul(a, a, tmp);


// 	mpz_mul(X, X, Zt);
// 	mpz_mul(Y, Y, Zt);
	

// 	// (x,y,z) = Q
// 	mpz_mul_ui(x, x, 3);
// 	mpz_add(x, x, A);
// 	mpz_mul_ui(y, Y, 3);
// 	mpz_mul_ui(z, B, 3);

// 	wDBL(d, x1, y1, z1, x, y, z, N, a); // (x1,y1,z1) = [2](Q)
	
// 	if (mpz_cmp_ui(d, 1)){
// 		goto final_step;
// 	}

// 	mpz_set(powX[0], x1);
// 	mpz_set(powY[0], y1);
// 	mpz_set(powZ[0], z1);

// 	// pow[2..maxd-1] ; pow[i] = [2*(i+1)](Q)
// 	for (i = 1; i < maxd - 1; i++){
// 		j = i + 1;
// 		wADD(d,powX[j],powY[j],powZ[j],powX[i],powY[i],powZ[i],x1,y1,z1,N,a);
// 		if (mpz_cmp_ui(d, 1)){
// 			goto final_step;
// 		}
// 	}
	
// 	// (x,y,z) = [p_1](Q)
// 	wLADDER_ui(d, x, y, z, primes[l1]+(diff[1]<<1), x, y, z, N, a); // compute [p_1] (Q)
	
// 	if (mpz_cmp_ui(d, 1)){
// 		goto final_step;
// 	}

// 	for (i = 2; i <= l2; i++){	
//         wADD(d, x, y, z, x, y, z, powX[diff[i]-1], powY[diff[i]-1], powZ[diff[i]-1], N, a);
    
// 		if (mpz_cmp_ui(d, 1)){
//             goto final_step;
// 		}
// 	}


final_step:
	
	/* free the memory */

	mpz_clears(X, Y, Z, x, y, z, x1, y1, z1, m2, k2, k, A, B, a24, mpz_4, Zt, X3, Z2, Z3, tmp, a, NULL);
	for(i=0; i< maxd; i++){
		mpz_clear(powX[i]);
		mpz_clear(powY[i]);
		mpz_clear(powZ[i]);
	}
	free(powX);
	free(powY);
	free(powZ);

	return;
}

/* User function : find a non-trivial big factor of N */
void factor(mpz_t d, mpz_t N, unsigned long B, unsigned long *primes, unsigned long *differences)
{
	printf("====================================================\n");
	if (mpz_cmp_ui(N, 0) == 0)
	{
		printf("Input is 0. Exit! \n \n");
		return;
	}
	if ((mpz_cmp_ui(N,1) == 0)|| mpz_probab_prime_p(N,25) > 0)
	{	
		printf("Factorized by trial division. Exit! \n \n");
		return;
	}
	gmp_printf("Factoring %Zd \n\n", N);
	trial_division(N, primes);
	gmp_printf("ECM : %Zd\n\n", N);
	clock_t st = clock();
	ECM_factor(d, N, B, primes, differences);
	printf("Elapsed time: %f \n\n", (double)(clock() - st) / CLOCKS_PER_SEC);
	return;
}