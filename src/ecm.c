#include "../headers/ecm.h"

/* Split N, using the bound B, e = lg(B), T is the table of primes up to B */
void ECM_factor(mpz_t d, mpz_t N, unsigned long B1, unsigned long * primes, unsigned long * diff)
{	
	// clock_t st;
	int flag;
	unsigned int i, m = 5;
	unsigned long q;

	mpz_t X, Z, u, v, u1, v1, A, B, a24, a, mpz_4, tmp;
	mpz_inits(X, Z, u, v, u1, v1, A, B, a24, a, mpz_4, tmp, NULL);
	
	mpz_set_ui(mpz_4, 4);
	mpz_invert(mpz_4, mpz_4, N); // mpz_4 = 1/4 mod N

	unsigned int l1 = primes[0], l2 = diff[0]; // number of primes for stage 1

	// largest difference in the array diff
	unsigned long maxd = diff[l2+1];
	
    /* allocate memory for the table of powers of Q in Stage 2 */
	mpz_t * powU = malloc(maxd * sizeof(mpz_t));
	mpz_t * powV = malloc(maxd * sizeof(mpz_t));

	for(i=0; i < maxd; i++){
		mpz_init(powU[i]);
		mpz_init(powV[i]);
	}
	
	/*---------------------------------------------*/
    
	/* choose the curve and base point using Suyama's parametrization */
loop_curve:
	
	m++;
	printf("[%u]",m);
	
	// Suyama_param(a24, A, X, Z, m, mpz_4, N);
	normal_param(a24, A, X, Z, m, mpz_4, N);

	flag = 0;

	goto stage1;

loop_point:
	if(flag == 1){
		mpz_set_ui(X, 2);
		mpz_set_ui(Z, 1);
	}else{
		goto loop_curve;
	}

	/* Stage 1 using Montgomery's ladder */

stage1:
	/* for loop computes M = prod(p_i^l_i) */
	
	for (i = 1; i <= l1; i++){
		q = primes[i];
		while ( q <= B1 ){
			q = q * primes[i];
		}

		ladder_ui(X, Z, q, X, Z, N, a24); // compute [q1](P)

		if(mpz_cmp_ui(Z,0) == 0){ // if Z = 0 then need to choose another point
			printf("here\n");
			flag = flag ^ 1;
			goto loop_point;
		}else{
			invert(X,d,X,Z,N);
			if (mpz_cmp_ui(d, 1)){
				goto final_step;
			}else{
				mpz_set_ui(Z,1);
			}
		}
	}
	printf("Stage 2\n");
	// (X : Z) = Q = [M](P) 

	/* move to Weierstrass curve: v^2 = u^3 + a * u + b
	 we set Y = 1, then Q = (X, 1, Z) */

	mpz_powm_ui(B, X, 3, N);
	mpz_mul(tmp, X, A);
	mpz_add(B, B, tmp);
	mpz_add(B, B, X);
	mpz_mod(B, B, N); // B = X * (X^2 + A * X + 1)
	
	mpz_powm_ui(tmp, B, 2, N);
	mpz_mul_ui(tmp, tmp, 3);
	mpz_invert(tmp, tmp, N);
	mpz_powm_ui(a, A, 2, N);
	mpz_sub_ui(a, a, 3);
	mpz_neg(a, a);
	mpz_mul(a, a, tmp);
	mpz_mod(a, a, N); // compute a = (3-A^2)/(3B^2)

	
	// t = X/B + A/3*B
	mpz_mul_ui(u,X,3);
	mpz_add(u,u,A);
	mpz_mul_ui(tmp,B,3);
	mpz_invert(tmp, tmp, N); 
	mpz_mul(u,u,tmp);
	mpz_mod(u,u,N);
	// v = Y/B = 1/B
	mpz_mul_ui(v,tmp,3);
	mpz_mod(v,v,N);

	wDBL(d, u1, v1, u, v, N, a); // (u1,v1,z1) = [2](Q)
	if (mpz_cmp_ui(d, 1)){
		goto final_step;
	}

	// pow[i] = [i+1] ([2](Q))
	mpz_set(powU[0], u1);
	mpz_set(powV[0], v1);
	for (i = 0; i < maxd - 1; i++){
		wADD(d,powU[i+1],powV[i+1],powU[i],powV[i],u1,v1,N,a);
		if (mpz_cmp_ui(d, 1)){
			goto final_step;
		}
	}
	
	// (u,v) = [p_1](Q)
	wLADDER_ui(d, u, v, primes[l1]+(diff[1]<<1), u, v, N, a); // compute [p_1] (Q)
	
	if (mpz_cmp_ui(d, 1)){
		goto final_step;
	}

	for (i = 2; i <= l2; i++){	
        wADD(d, u, v, powU[diff[i]-1], powV[diff[i]-1], u, v, N, a);
		if (mpz_cmp_ui(d, 1)){
            goto final_step;
		}
	}

	goto loop_curve;

final_step:
	if(mpz_cmp_ui(d,0) == 0){
		goto loop_curve;
	}
	
	/* free the memory */

	mpz_clears(X, Z, u, v, u1, v1, A, B, a24, a, mpz_4, tmp, NULL);
	for(i=0; i< maxd; i++){
		mpz_clear(powU[i]);
		mpz_clear(powV[i]);
	}
	free(powU);
	free(powV);
	return;
}