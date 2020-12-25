#include "../headers/ecm.h"

int main(int argc, char * argv[]){
	setbuf(stdout, NULL);

	// if (argc == 1){
	// 	return 0;
	// }
	
	// srand((unsigned) time(NULL));

	// char * filename = argv[1];
	// FILE * fp;
	// fp = fopen(filename, "r");
	// if (fp == NULL){
	// 	exit(EXIT_FAILURE);
	// }
	// printf("Reading input from './%s'\n", filename);
	
	// char * line = NULL;
	// size_t len = 0;
 //    ssize_t read;
	// read = getline(&line, &len, fp);
	// int l = atoi(line), i = 0;
	// mpz_t * input = malloc(sizeof(mpz_t)*l);
	// while ((read = getline(&line, &len, fp)) != -1) {
 //        mpz_set_str(input[i],line,10);
	// 	i++;
 //    }
	// fclose(fp);
	// if (argc >= 3) B1 = atoi(argv[2]);
	// if (argc == 4) B2 = atoi(argv[3]);

	unsigned long B1 = 100000, B2 = 2000000;

	
	printf("B1 = %lu , B2 = %lu\n", B1, B2);
	unsigned long *primes,*differences;
	unsigned long size_array = B2/2; /* change the size of array - using PNT ln(B) */
	primes = (unsigned long *)malloc(size_array * sizeof(unsigned long));
	differences = (unsigned long *)malloc(size_array * sizeof(unsigned long));
	
	// precompute primes up to limit and  
	precompute(primes,B1,differences,B2);
	printf("====================================================\n\n\n");
	
	// for(i = 0 ; i < l; i++){
	// 	factor(input[i], B1, primes,differences);
	// }
	mpz_t d, N;
	mpz_inits(d, N, NULL);
	mpz_set_str(N, "8124698401622845318488774792861997627459268027692201251242048839048",10);
	factor(d, N, B1, primes, differences);
	gmp_printf("%Zd \n",d);
	// mpz_t X1, Z1, N, a24, X, Z, X2, Z2, X3, Z3, m;
	// mpz_inits(X1,Z1,N,a24,X,Z,X2, Z2, X3, Z3, m, NULL);
	// mpz_set_str(N,"1067656525553",10);
	// mpz_set_str(a24,"498699635029",10);
	// mpz_set_str(X,"750406921487",10);
	// mpz_set_ui(Z,15);
	// mpz_set_ui(X2,19);
	// mpz_set_ui(Z2,29);
	// mpz_set_ui(X3,9);
	// mpz_set_ui(Z3,5);
	// mpz_set_str(m, "4999",10);
	// ladder( X1, Z1, N, a24, m, X);
	// gmp_printf(" X1 = %Zd Z1 = %Zd \n", X1, Z1);
	// mpz_clears(X1,Z1,N,a24,X,Z, X2, Z2, X3, Z3, m, NULL);
	mpz_clears(d, N, NULL);
	free(primes);
	free(differences);
	
	return 0;
}