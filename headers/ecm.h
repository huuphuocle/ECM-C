#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* auxiliary.c */

void trial_division(mpz_t N, unsigned long *primes);
int is_composite(mpz_t N, int c);
void precompute(unsigned long *T, unsigned long B1, unsigned long *D, unsigned long B2);

/* montgomery.c */

void printmCurve(mpz_t A, mpz_t B);
void printmPoint(mpz_t X, mpz_t Z);
void Suyama_param(mpz_t a24, mpz_t A, mpz_t X, mpz_t Z, unsigned int m, mpz_t mpz_4, mpz_t N);
void normal_param(mpz_t a24, mpz_t A, mpz_t X, mpz_t Z, unsigned int m, mpz_t mpz_4, mpz_t N);
void invert(mpz_t X, mpz_t d, mpz_t x, mpz_t z, mpz_t N);
void xADD(mpz_t X, mpz_t Z, mpz_t X1, mpz_t Z1, mpz_t X2, mpz_t Z2, mpz_t X3, mpz_t Z3, mpz_t N);
// void xADD2(mpz_t X, mpz_t Z, mpz_t X1, mpz_t Z1, mpz_t X2, mpz_t Z2, mpz_t X3, mpz_t  N);
void xDBL(mpz_t X, mpz_t Z, mpz_t X1, mpz_t Z1, mpz_t N, mpz_t a24);
void ladder(mpz_t X, mpz_t Z, mpz_t m, mpz_t X1, mpz_t Z1, mpz_t N, mpz_t a24);
void ladder_ui(mpz_t X, mpz_t Z, unsigned long m, mpz_t X1, mpz_t Z1, mpz_t N, mpz_t a24);
// void ladder2_ui(mpz_t X, mpz_t Z, unsigned long m, mpz_t X1, mpz_t N, mpz_t a24);
// void ladder2(mpz_t X, mpz_t Z, mpz_t m, mpz_t X1, mpz_t N, mpz_t a24);

/* weierstrass.c */

void printwCurve(mpz_t a, mpz_t b);
void printwPoint(mpz_t X, mpz_t Y, mpz_t Z);
int wDBL(mpz_t d, mpz_t X, mpz_t Y, mpz_t X1, mpz_t Y1, mpz_t N, mpz_t a);
int wADD(mpz_t d, mpz_t X, mpz_t Y, mpz_t X1, mpz_t Y1, mpz_t X2, mpz_t Y2, mpz_t N, mpz_t a);
int wLADDER(mpz_t d, mpz_t X, mpz_t Y, mpz_t m, mpz_t X1, mpz_t Y1, mpz_t N, mpz_t a);
int wLADDER_ui(mpz_t d, mpz_t X, mpz_t Y, unsigned int m, mpz_t X1, mpz_t Y1, mpz_t N, mpz_t a);
int on_wCurve(mpz_t a, mpz_t b, mpz_t x, mpz_t y, mpz_t N);
// void computeB(mpz_t b, mpz_t a, mpz_t x, mpz_t y, mpz_t N);

/* ecm.c */

void ECM_factor(mpz_t d, mpz_t N, unsigned long B, unsigned long * primes, unsigned long * diff);
void ECM_factor2(mpz_t d, mpz_t N, unsigned long B1, unsigned long * primes, unsigned long * diff);
void ECM_factor_Suyama(mpz_t d, mpz_t N, unsigned long B, unsigned long * primes, unsigned long * diff);
void ECM_factor2_Suyama(mpz_t d, mpz_t N, unsigned long B1, unsigned long * primes, unsigned long * diff);

/* user.c */

void factor(mpz_t d, mpz_t N, unsigned long B, unsigned long *primes, unsigned long *differences, int stage, int param);
