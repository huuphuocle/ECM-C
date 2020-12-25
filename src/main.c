#include "../headers/ecm.h"

int main(int argc, char * argv[]){
	setbuf(stdout, NULL);

	if (argc == 1){
		return 0;
	}
	
	srand((unsigned) time(NULL));

	char * filename = argv[1];
	FILE * fp;
	fp = fopen(filename, "r");
	if (fp == NULL){
		exit(EXIT_FAILURE);
	}
	
	char * line = NULL;
	size_t len = 0;
	
	mpz_t d, N;
	mpz_inits(d, N, NULL);

	getline(&line, &len, fp);
    mpz_set_str(N,line,10);
	
	fclose(fp);
	
	unsigned long B1, B2;

	if (argc >= 3){
		B1 = atoi(argv[2]);
		if (argc == 4)
			B2 = atoi(argv[3]);
	}else{
		B1 = 100000;
		B2 = 5000000;
	}
	
	unsigned long size_array = (int) (1.5*B2/log(B2));
	unsigned long *primes,*differences;
	primes = (unsigned long *)malloc(size_array * sizeof(unsigned long));
	differences = (unsigned long *)malloc(size_array * sizeof(unsigned long));
	
	// precompute primes up to limit B1 and B2 
	precompute(primes,B1,differences,B2);
	
	factor(d, N, B1, primes, differences);

	mpz_clears(d, N, NULL);
	free(primes);
	free(differences);
	
	return 0;
}