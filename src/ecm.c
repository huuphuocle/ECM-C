#include "../headers/ecm.h"

/* Split N, using the bound B, e = lg(B), T is the table of primes up to B */
void ECM_factor(mpz_t d, mpz_t N, unsigned long B1, unsigned long * primes, unsigned long * diff){
	printf("Use only stage 1 with normal parametrization\n");
	unsigned int i, m = 2;
	unsigned long q;

	mpz_t X, Z, u, v, u1, v1, A, B, a24, a, b, mpz_4, tmp;
	mpz_inits(X, Z, u, v, u1, v1, A, B, a24, a, b, mpz_4, tmp, NULL);
	
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
	printf("\r[%u]",m);
	// Suyama_param(a24, A, X, Z, m, mpz_4, N);
	normal_param(a24, A, X, Z, m, mpz_4, N);
	
	for (i = 1; i <= l1; i++){
		q = primes[i];
		while ( q <= B1 ){
			q = q * primes[i];
		}
		ladder_ui(X, Z, q, X, Z, N, a24); // compute [q1](P)

		if(mpz_cmp_ui(Z,0) == 0){
			goto loop_curve;
		}else{
			invert(X,d,X,Z,N);
			if (mpz_cmp_ui(d, 1)){
				goto final_step;
			}else{
				mpz_set_ui(Z,1);
			}
		}
	}
	// if(mpz_cmp_ui(Z,0) == 0){
	// 	goto loop_curve;
	// }else{
	// 	invert(X,d,X,Z,N);
	// 	if (mpz_cmp_ui(d, 1)){
	// 		goto final_step;
	// 	}else{
	// 		mpz_set_ui(Z,1);
	// 	}
	// }

	goto loop_curve;

final_step:
	if(mpz_cmp_ui(d,0) == 0) goto loop_curve;
	
	/* free the memory */

	mpz_clears(X, Z, u, v, u1, v1, A, B, a24, a, b, mpz_4, tmp, NULL);
	for(i=0; i< maxd; i++){
		mpz_clear(powU[i]);
		mpz_clear(powV[i]);
	}
	free(powU);
	free(powV);
	return;
}

void ECM_factor2(mpz_t d, mpz_t N, unsigned long B1, unsigned long * primes, unsigned long * diff){	
	printf("Use stage 1 and 2 with normal parametrization\n");
	int flag, switch_curve = 0;
	unsigned int i, m = 2;
	unsigned long q;

	mpz_t X, Z, u, v, u1, v1, A, B, a24, a, b, mpz_4, tmp;
	mpz_inits(X, Z, u, v, u1, v1, A, B, a24, a, b, mpz_4, tmp, NULL);
	
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
	// if(switch_curve == 0){
	// 	m++;
	// 	printf("[%u]",m);
	// 	normal_param(a24, A, X, Z, m, mpz_4, N);
	// }else Suyama_param(a24, A, X, Z, m, mpz_4, N);

	m++;
	printf("\r[%u]",m);
	// Suyama_param(a24, A, X, Z, m, mpz_4, N);
	normal_param(a24, A, X, Z, m, mpz_4, N);

/* stage 1 using Montgomery's ladder */
	
	for (i = 1; i <= l1; i++){
		q = primes[i];
		while ( q <= B1 ){
			q = q * primes[i];
		}
		ladder_ui(X, Z, q, X, Z, N, a24); // compute [q1](P)

		if(mpz_cmp_ui(Z,0) == 0){ // if Z = 0 then need to choose another point
			switch_curve = switch_curve ^ 1;
			goto loop_curve;
		}else{
			invert(X,d,X,Z,N);
			if (mpz_cmp_ui(d, 1)){
				goto final_step;
			}else{
				mpz_set_ui(Z,1);
			}
		}
	}

	/* move to Weierstrass curve: v^2 = u^3 + a * u + b
	  we set Y = 1, then Q = (X, 1, Z) */

	mpz_powm_ui(B, X, 2, N);
	mpz_mul(tmp, X, A);
	mpz_add(B, B, tmp);
	mpz_add_ui(B, B, 1);
	mpz_mul(B, B, X);
	mpz_mod(B, B, N); // B = X * (X^2 + A * X + 1)

	mpz_mul_ui(v,B,3);
	mpz_invert(v, v, N); // v = 1/(3 * B)

	mpz_powm_ui(tmp, v, 2, N);
	mpz_powm_ui(a, A, 2, N);
	mpz_sub_ui(a, a, 3);
	mpz_neg(a, a);
	mpz_mul(a, a, tmp);
	mpz_mul_ui(a, a, 3);
	mpz_mod(a, a, N); // a = (3-A^2)/(3B^2)

	// mpz_powm_ui(b,A,3,N);
	// mpz_mul_ui(b,b,2);
	// mpz_mul_ui(tmp, A, 9);
	// mpz_sub(b, b, tmp);
	// mpz_powm_ui(tmp, v, 3, N);
	// mpz_mul(b, b, tmp);
	// mpz_mod(b, b, N); // b = (2*A^3 - 9*A)/(27*B^3)

	// t = X/B + A/3*B
	mpz_mul_ui(u, X, 3);
	mpz_add(u, u, A);
	mpz_mul(u, u, v);
	mpz_mod(u, u, N);
	
	// v = Y/B = 1/B
	mpz_mul_ui(v, v, 3);
	mpz_mod(v, v, N);

	// Q = (u,v) = [M] P is a point on the Weierstrass' curve v^2 = u^3 + a*u + b
	flag = wDBL(d, u1, v1, u, v, N, a); // (u1,v1) = [2] Q
	if (flag == 1) goto final_step;

	// pow[i] = [i+1] Q' = [2(i+1)] Q
	mpz_set(powU[0], u1);
	mpz_set(powV[0], v1);
	for (i = 0; i < maxd - 1; i++){
		flag = wADD(d,powU[i+1],powV[i+1],powU[i],powV[i],u1,v1,N,a);
		if (flag == 1) goto final_step;
	}

	// (u,v) = [p_1](Q)
	flag = wLADDER_ui(d, u, v, primes[l1]+(diff[1]<<1), u, v, N, a); // compute [p_1] Q
	if (flag == 1) goto final_step;

	for (i = 2; i <= l2; i++){	
        flag = wADD(d, u, v, powU[diff[i]-1], powV[diff[i]-1], u, v, N, a);
		if (flag == 1) goto final_step;
	}

	switch_curve = switch_curve ^ 1;
	goto loop_curve;

final_step:
	if(mpz_cmp_ui(d,0) == 0){
		switch_curve = switch_curve ^ 1;
		goto loop_curve;
	}
	
	/* free the memory */

	mpz_clears(X, Z, u, v, u1, v1, A, B, a24, a, b, mpz_4, tmp, NULL);
	for(i=0; i< maxd; i++){
		mpz_clear(powU[i]);
		mpz_clear(powV[i]);
	}
	free(powU);
	free(powV);
	return;
}

#include "../headers/ecm.h"

/* Split N, using the bound B, e = lg(B), T is the table of primes up to B */
void ECM_factor_Suyama(mpz_t d, mpz_t N, unsigned long B1, unsigned long * primes, unsigned long * diff){
	printf("Use only stage 1 with Suyama\n");
	unsigned int i, m = 2;
	unsigned long q;

	mpz_t X, Z, u, v, u1, v1, A, B, a24, a, b, mpz_4, tmp;
	mpz_inits(X, Z, u, v, u1, v1, A, B, a24, a, b, mpz_4, tmp, NULL);
	
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
	printf("\r[%u]",m);
	Suyama_param(a24, A, X, Z, m, mpz_4, N);
	// normal_param(a24, A, X, Z, m, mpz_4, N);
	
	for (i = 1; i <= l1; i++){
		q = primes[i];
		while ( q <= B1 ){
			q = q * primes[i];
		}
		ladder_ui(X, Z, q, X, Z, N, a24); // compute [q1](P)

		if(mpz_cmp_ui(Z,0) == 0){
			goto loop_curve;
		}else{
			invert(X,d,X,Z,N);
			if (mpz_cmp_ui(d, 1)){
				goto final_step;
			}else{
				mpz_set_ui(Z,1);
			}
		}
	}
	// if(mpz_cmp_ui(Z,0) == 0){
	// 	goto loop_curve;
	// }else{
	// 	invert(X,d,X,Z,N);
	// 	if (mpz_cmp_ui(d, 1)){
	// 		goto final_step;
	// 	}else{
	// 		mpz_set_ui(Z,1);
	// 	}
	// }

	goto loop_curve;

final_step:
	if(mpz_cmp_ui(d,0) == 0) goto loop_curve;
	
	/* free the memory */

	mpz_clears(X, Z, u, v, u1, v1, A, B, a24, a, b, mpz_4, tmp, NULL);
	for(i=0; i< maxd; i++){
		mpz_clear(powU[i]);
		mpz_clear(powV[i]);
	}
	free(powU);
	free(powV);
	return;
}

void ECM_factor2_Suyama(mpz_t d, mpz_t N, unsigned long B1, unsigned long * primes, unsigned long * diff){	
	printf("Use stage 1 and 2 with Suyama\n");
	int flag, switch_curve = 0;
	unsigned int i, m = 2;
	unsigned long q;

	mpz_t X, Z, u, v, u1, v1, A, B, a24, a, b, mpz_4, tmp;
	mpz_inits(X, Z, u, v, u1, v1, A, B, a24, a, b, mpz_4, tmp, NULL);
	
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
	// if(switch_curve == 0){
	// 	m++;
	// 	printf("[%u]",m);
	// 	normal_param(a24, A, X, Z, m, mpz_4, N);
	// }else Suyama_param(a24, A, X, Z, m, mpz_4, N);

	m++;
	printf("\r[%u]",m);
	Suyama_param(a24, A, X, Z, m, mpz_4, N);
	// normal_param(a24, A, X, Z, m, mpz_4, N);

/* stage 1 using Montgomery's ladder */
	
	for (i = 1; i <= l1; i++){
		q = primes[i];
		while ( q <= B1 ){
			q = q * primes[i];
		}
		ladder_ui(X, Z, q, X, Z, N, a24); // compute [q1](P)

		if(mpz_cmp_ui(Z,0) == 0){ // if Z = 0 then need to choose another point
			switch_curve = switch_curve ^ 1;
			goto loop_curve;
		}else{
			invert(X,d,X,Z,N);
			if (mpz_cmp_ui(d, 1)){
				goto final_step;
			}else{
				mpz_set_ui(Z,1);
			}
		}
	}

	/* move to Weierstrass curve: v^2 = u^3 + a * u + b
	  we set Y = 1, then Q = (X, 1, Z) */

	mpz_powm_ui(B, X, 2, N);
	mpz_mul(tmp, X, A);
	mpz_add(B, B, tmp);
	mpz_add_ui(B, B, 1);
	mpz_mul(B, B, X);
	mpz_mod(B, B, N); // B = X * (X^2 + A * X + 1)

	mpz_mul_ui(v,B,3);
	mpz_invert(v, v, N); // v = 1/(3 * B)

	mpz_powm_ui(tmp, v, 2, N);
	mpz_powm_ui(a, A, 2, N);
	mpz_sub_ui(a, a, 3);
	mpz_neg(a, a);
	mpz_mul(a, a, tmp);
	mpz_mul_ui(a, a, 3);
	mpz_mod(a, a, N); // a = (3-A^2)/(3B^2)

	// mpz_powm_ui(b,A,3,N);
	// mpz_mul_ui(b,b,2);
	// mpz_mul_ui(tmp, A, 9);
	// mpz_sub(b, b, tmp);
	// mpz_powm_ui(tmp, v, 3, N);
	// mpz_mul(b, b, tmp);
	// mpz_mod(b, b, N); // b = (2*A^3 - 9*A)/(27*B^3)

	// t = X/B + A/3*B
	mpz_mul_ui(u, X, 3);
	mpz_add(u, u, A);
	mpz_mul(u, u, v);
	mpz_mod(u, u, N);
	
	// v = Y/B = 1/B
	mpz_mul_ui(v, v, 3);
	mpz_mod(v, v, N);

	// Q = (u,v) = [M] P is a point on the Weierstrass' curve v^2 = u^3 + a*u + b
	flag = wDBL(d, u1, v1, u, v, N, a); // (u1,v1) = [2] Q
	if (flag == 1) goto final_step;

	// pow[i] = [i+1] Q' = [2(i+1)] Q
	mpz_set(powU[0], u1);
	mpz_set(powV[0], v1);
	for (i = 0; i < maxd - 1; i++){
		flag = wADD(d,powU[i+1],powV[i+1],powU[i],powV[i],u1,v1,N,a);
		if (flag == 1) goto final_step;
	}

	// (u,v) = [p_1](Q)
	flag = wLADDER_ui(d, u, v, primes[l1]+(diff[1]<<1), u, v, N, a); // compute [p_1] Q
	if (flag == 1) goto final_step;

	for (i = 2; i <= l2; i++){	
        flag = wADD(d, u, v, powU[diff[i]-1], powV[diff[i]-1], u, v, N, a);
		if (flag == 1) goto final_step;
	}

	switch_curve = switch_curve ^ 1;
	goto loop_curve;

final_step:
	if(mpz_cmp_ui(d,0) == 0){
		switch_curve = switch_curve ^ 1;
		goto loop_curve;
	}
	
	/* free the memory */

	mpz_clears(X, Z, u, v, u1, v1, A, B, a24, a, b, mpz_4, tmp, NULL);
	for(i=0; i< maxd; i++){
		mpz_clear(powU[i]);
		mpz_clear(powV[i]);
	}
	free(powU);
	free(powV);
	return;
}