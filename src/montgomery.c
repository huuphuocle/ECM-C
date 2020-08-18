#include "../headers/ecm.h"

void printmCurve(mpz_t A, mpz_t B){
    gmp_printf("Montgomery curve %Zd Y^2 = X^3 + %Zd X^2 + X\n",B, A);
    return;
}

void printmPoint(mpz_t X, mpz_t Z){
    gmp_printf("Point (%Zd, %Zd) \n",X,Z);
    return;
}

void Suyama_param(mpz_t a24, mpz_t A, mpz_t X, mpz_t Z, unsigned int m, mpz_t N){
    mpz_t u, v, tmp;
    mpz_inits(u, v, tmp, NULL);

    mpz_set_ui(u, m*m-5); // u = m^2-5
    mpz_set_ui(v, m << 2); // v = 4m
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

    mpz_clears(u, v, tmp, NULL);
    return;
}

void normal_param(mpz_t a24, mpz_t A, mpz_t X, mpz_t Z, unsigned int m, mpz_t N, mpz_t mpz_4){
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

/* [x:z] -> [X:1] */

void normalisation(mpz_t X, mpz_t x, mpz_t z, mpz_t N){
	mpz_t X0;
	mpz_init(X0);
    mpz_invert(X0, z, N);
    mpz_mul(X0,X0,x);
    mpz_mod(X,X0,N);
    mpz_clear(X0);
    return;
}

/* compute (X,Z) = (X1,Z1)+(X2,Z2), with (X3,Z3) = (X1,Z1)-(X2,Z2) */

void xADD(mpz_t X, mpz_t Z, mpz_t X1, mpz_t Z1, mpz_t X2, mpz_t Z2, mpz_t X3, mpz_t Z3, mpz_t N){
    mpz_t X0, Z0, a, b, tmp;
    mpz_inits(X0, Z0, a,b,tmp,NULL);
    mpz_sub(a,X1,Z1); // a = X1 - Z1
    mpz_add(tmp,X2,Z2); // tmp = X2 + Z2
    mpz_mul(a,a,tmp); // a = (X1 - Z1)*(X2 + Z2)
    mpz_mod(a,a,N); 
    mpz_add(b,X1,Z1);
    mpz_sub(tmp,X2,Z2);
    mpz_mul(b,b,tmp); // b = (X1 + Z1)*(X2 - Z2)
    mpz_mod(b,b,N);
    /*
	a = ((X1-Z1)*(X2+Z2)) % N
	b = ((X1+Z1)*(X2-Z2)) % N
    */
    
    mpz_add(X0,a,b);
    mpz_powm_ui(X0, X0, 2, N);
    mpz_mul(X0, X0, Z3);
    mpz_mod(X0, X0, N); // X0 = Z3 * (a + b)^2

    mpz_sub(Z0, a, b);
    mpz_powm_ui(Z0, Z0, 2, N);
    mpz_mul(Z0, Z0, X3);
    mpz_mod(Z0,Z0,N); // Z0 = X3 * (a-b)^2

    mpz_set(X,X0);
    mpz_set(Z,Z0);

    mpz_clears(X0, Z0, a,b,tmp,NULL);
	return;
}

/* xADD with (X3,Z3) fixed */

void xADD2(mpz_t X, mpz_t Z, mpz_t X1, mpz_t Z1, mpz_t X2, mpz_t Z2, mpz_t X3, mpz_t  N){
    mpz_t X0, Z0, a, b, tmp;
    mpz_inits(X0, Z0, a,b,tmp,NULL);
    mpz_sub(a,X1,Z1);
    mpz_add(tmp,X2,Z2);
    mpz_mul(a,a,tmp);
    mpz_mod(a,a,N);
    mpz_add(b,X1,Z1);
    mpz_sub(tmp,X2,Z2);
    mpz_mul(b,b,tmp);
    mpz_mod(b,b,N);

    mpz_add(X0,a,b);
    mpz_powm_ui(X0, X0, 2, N); // one multiplication saved here as Z3 = 1

    mpz_sub(Z0, a, b);
    mpz_powm_ui(Z0, Z0, 2, N);
    mpz_mul(Z0, Z0, X3);
    mpz_mod(Z0, Z0, N);

    mpz_set(X, X0);
    mpz_set(Z, Z0);

    mpz_clears(X0, Z0, a,b,tmp,NULL);
	return;
}

/* xDBL : compute (X1,Z1) = [2](X,Z) mod N */

void xDBL(mpz_t X, mpz_t Z, mpz_t X1, mpz_t Z1, mpz_t N, mpz_t a24){
    mpz_t X0, Z0, a, b, c, tmp;
    mpz_inits(X0, Z0, a, b, c, tmp, NULL);
    
    mpz_add(a, X1, Z1); 
    mpz_powm_ui(a, a, 2, N); // a = (X1 + Z1)^2
    
    mpz_sub(b, X1, Z1); 
    mpz_powm_ui(b, b, 2, N); // b = (X1 - Z1)^2
    mpz_mul(X0,a,b);
    mpz_mod(X0,X0,N); // X = ((X1 + Z1) * (X1 - Z1))^2

    mpz_sub(c,a,b);
    mpz_mod(c,c,N); // c = 4 * X1 * Z1
    mpz_mul(tmp,c,a24);
    mpz_add(tmp,tmp,b);
    mpz_mul(Z0,c,tmp);
    mpz_mod(Z0,Z0,N); // Z1 = (c * (b + a24*c)) mod N

    mpz_set(X, X0);
    mpz_set(Z, Z0);
    mpz_clears(X0, Z0, a, b, c, tmp, NULL);
    return;
}

/* Montgomery ladder : (X,Z) = [m] (X1,1) : m <> 0 */

void ladder2(mpz_t X, mpz_t Z, mpz_t m, mpz_t X1, mpz_t N, mpz_t a24){
    mpz_t X0, Z0, X2, Z2, m_t;
    mpz_inits(X0, Z0, X2, Z2, m_t, NULL);
    mpz_set_ui(X0, 1);
    mpz_set_ui(Z0, 0);
    mpz_set(X2, X1);
    mpz_set_ui(Z2, 1);
    unsigned int swap = 0;

    unsigned int bcnt = mpz_sizeinbase(m, 2);
    
    for(int t = bcnt - 1; t >= 0; t--){
        mpz_fdiv_q_2exp(m_t, m, t);
        swap = mpz_get_ui(m_t) & 1;
        if (swap){
            mpz_swap(X2, X0);
            mpz_swap(Z2, Z0);
        }
        xADD2(X2, Z2, N, X0, Z0, X2, Z2, X1);
        xDBL(X0, Z0, N, a24, X0, Z0);
        if (swap){
            mpz_swap(X2, X0);
            mpz_swap(Z2, Z0);
        }
    }
    mpz_set(X, X0);
    mpz_set(Z, Z0);
    mpz_clears(X0, Z0, X2, Z2, m_t, NULL);
    return;
}

void ladder(mpz_t X, mpz_t Z, mpz_t m, mpz_t x, mpz_t z, mpz_t N, mpz_t a24){
    mpz_t X0, Z0, X1, Z1, mt;
    mpz_inits(X0, Z0, X1, Z1, mt, NULL);
    mpz_set(X0, x);
    mpz_set(Z0, z);
    xDBL(X1, Z1, x, z, N, a24);
    unsigned int  m_t;

    unsigned int bcnt = mpz_sizeinbase(m, 2);
    
    for(int t = bcnt - 2; t >= 0; t--){
        mpz_fdiv_q_2exp(mt, m, t);
        m_t = mpz_get_ui(mt) & 1;
        if (m_t == 0){
            xADD(X1, Z1, X1, Z1, X0, Z0, x, z, N);
            xDBL(X0, Z0, X0, Z0, N, a24);
        }else{
            xADD(X0, Z0, X0, Z0, X1, Z1, x, z, N);
            xDBL(X1, Z1, X1, Z1, N, a24);
        }
        
    }
    
    mpz_set(X, X0);
    mpz_set(Z, Z0);
    mpz_clears(X0, Z0, X1, Z1, mt, NULL);
    return;
}

void ladder_ui(mpz_t X, mpz_t Z, unsigned long m, mpz_t x, mpz_t z, mpz_t N, mpz_t a24){
    mpz_t X0, Z0, X1, Z1;
    mpz_inits(X0, Z0, X1, Z1, NULL);
    mpz_set(X0, x);
    mpz_set(Z0, z);
    xDBL(X1, Z1, x, z, N, a24);
    unsigned int  m_t = m, bcnt = 0;
    
    while(m_t > 0){
        bcnt++;
        m_t = m_t >> 1;
    }
    
    for(int t = bcnt - 2; t >= 0; t--){
        m_t = (m >> t) & 1;
        if (m_t == 0){
            xADD(X1, Z1, X1, Z1, X0, Z0, x, z, N);
            xDBL(X0, Z0, X0, Z0, N, a24);
        }else{
            xADD(X0, Z0, X0, Z0, X1, Z1, x, z, N);
            xDBL(X1, Z1, X1, Z1, N, a24);
        }
        
    }
    
    mpz_set(X, X0);
    mpz_set(Z, Z0);
    mpz_clears(X0, Z0, X1, Z1, NULL);
    return;
}

void ladder2_ui(mpz_t X, mpz_t Z, unsigned long m, mpz_t X1, mpz_t N, mpz_t a24){
    mpz_t X0, Z0, X2, Z2;
    mpz_inits(X0, Z0, X2, Z2, NULL);
    mpz_set_ui(X0, 1);
    mpz_set_ui(Z0, 0);
    mpz_set(X2, X1);
    mpz_set_ui(Z2, 1);
    unsigned int  m_t = m, bcnt = 0;
    
    while(m_t > 0){
        bcnt++;
        m_t = m_t >> 1;
    }
    
    for(int t = bcnt - 1; t >= 0; t--){
        m_t = (m >> t) & 1;
        if (m_t){
            mpz_swap(X2, X0);
            mpz_swap(Z2, Z0);
        }
        xADD2(X2, Z2, X0, Z0, X2, Z2, X1, N);
        xDBL(X0, Z0, X0, Z0, N, a24);
        if (m_t){
            mpz_swap(X2, X0);
            mpz_swap(Z2, Z0);
        }
    }
    
    mpz_set(X, X0);
    mpz_set(Z, Z0);
    mpz_clears(X0, Z0, X2, Z2, NULL);
    return;
}

// int main(){
//     mpz_t a, b, a24, x1, y1, z1, x2, y2, z2, p, x0, y0, z0;
//     mpz_inits(a,b,a24, x1,y1,z1,x2,y2,z2,p,x0,y0,z0,NULL);
//     mpz_set_ui(a, 1);
//     mpz_set_ui(b, 1);
//     mpz_set_ui(a24, 2);
//     mpz_set_ui(p, 5);
//     mpz_set_ui(x1, 4);
//     mpz_set_ui(y1, 2);
//     mpz_set_ui(z1, 1);
//     mpz_set_ui(x2, 0);
//     mpz_set_ui(y2, 4);
//     mpz_set_ui(z2, 1);
//     xDBL(x0,z0,x1,z1,p,a24);
//     // xADD(x0,z0,x1,z1,p,a24);
//     printmCurve(a,b);
//     printPoint(x0,z0);
//     mpz_clears(a,b,x1,y1,z1,x2,y2,z2,p,x0,y0,z0,NULL);
// }