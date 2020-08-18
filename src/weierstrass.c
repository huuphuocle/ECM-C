/* Weierstrass curve */
#include "../headers/ecm.h"

void printwCurve(mpz_t a, mpz_t b){
    gmp_printf("Weierstrass curve Y^2 = X^3 + %Zd X + %Zd\n",a,b);
    return;
}

void printPoint(mpz_t X, mpz_t Y, mpz_t Z){
    gmp_printf("Point (%Zd, %Zd, %Zd) \n",X,Y,Z);
    return;
}

void computeB(mpz_t b, mpz_t a, mpz_t x, mpz_t y, mpz_t N){
    mpz_t y2, x3;
    mpz_inits(y2, x3, NULL);
    mpz_powm_ui(y2, y, 2, N);
    mpz_powm_ui(x3, x, 3, N);
    mpz_sub(y2, y2, x3);
    mpz_mul(x3, x, a);
    mpz_sub(y2, y2, x3);
    mpz_mod(b, y2, N);
    mpz_clears(y2, x3, NULL);
}

void on_wCurve(mpz_t a, mpz_t b, mpz_t x, mpz_t y, mpz_t N){
    mpz_t y2, x3;
    mpz_inits(y2, x3, NULL);
    mpz_powm_ui(y2, y, 2, N);
    mpz_powm_ui(x3, x, 3, N);
    mpz_sub(y2, y2, x3);
    mpz_mul(x3, x, a);
    mpz_sub(y2, y2, x3);
    mpz_sub(y2, y2, b);
    mpz_mod(y2, y2, N);
    if(mpz_cmp_ui(y2, 0) == 0){
        //printf("true\n");
    }else{
        printf("false\n");
    }
    mpz_clears(y2, x3, NULL);
}

/* double a point in Weierstrass curve */

void wDBL(mpz_t d, mpz_t X, mpz_t Y, mpz_t Z, mpz_t X1, mpz_t Y1, mpz_t Z1, mpz_t N, mpz_t a){
    
    if(mpz_cmp_ui(Z1,0) == 0){
        mpz_set_ui(X, 0);
        mpz_set_ui(Y, 1);
        mpz_set_ui(Z, 0);
        return;
    }
    if(mpz_cmp_ui(Y1,0)==0){
        mpz_set_ui(X, 0);
        mpz_set_ui(Y, 1);
        mpz_set_ui(Z, 0);
        return;
    }
    mpz_t X0, Y0, s, tmp, tmp2;

    mpz_inits(X0, Y0, s, tmp, tmp2, NULL);
    
    mpz_powm_ui(s, X1, 2, N);
    mpz_mul_ui(s, s, 3);
    mpz_add(s, s, a);
    mpz_mul_2exp(tmp, Y1, 1);
    
    mpz_gcdext (d, tmp, tmp2, tmp, N); 
    // tmp * 2Y1 + tmp2 * N = d
    if (mpz_cmp_ui(d, 1) != 0){
        return;
    }

    mpz_mul(s, s, tmp);
    mpz_mod(s, s, N); // s = (3*X1^2+a)/2*Y1

    mpz_mul_2exp(tmp, X1, 1);
    mpz_powm_ui(X0, s, 2, N);
    mpz_sub(X0, X0, tmp);
    mpz_mod(X0, X0, N); // X = s^2 - (X1+X2)

    mpz_sub(Y0, X1, X0);
    mpz_mul(Y0, Y0, s);
    mpz_sub(Y0, Y0, Y1);
    mpz_mod(Y0, Y0, N); // Y = s*(X1-X) - Y1

    mpz_set(X, X0);
    mpz_set(Y, Y0);
    mpz_set_ui(Z, 1);
    mpz_clears(X0, Y0, s, tmp, tmp2, NULL);
    return;
}

/* add two points on a Weierstrass curve */

void wADD(mpz_t d, mpz_t X, mpz_t Y, mpz_t Z, mpz_t X1, mpz_t Y1, mpz_t Z1, mpz_t X2, mpz_t Y2, mpz_t Z2, mpz_t N, mpz_t a){

    if(mpz_cmp_ui(Z1, 0) == 0){
        mpz_set(X, X2);
        mpz_set(Y, Y2);
        mpz_set(Z, Z2);
        return;
    }
    if(mpz_cmp_ui(Z2, 0) == 0){
        mpz_set(X, X1);
        mpz_set(Y, Y1);
        mpz_set(Z, Z1);
        return;
    }
    if(mpz_cmp(X1, X2) == 0){
        if(mpz_cmp(Y1, Y2) != 0){
            mpz_set_ui(X, 0);
            mpz_set_ui(Y, 1);
            mpz_set_ui(Z, 0);
            return;
        }else{
            wDBL(d, X, Y, Z, X1, Y1, Z1, N, a);
            return;
        }
    }
    mpz_t X0, Y0, s, tmp, tmp2;

    mpz_inits(X0, Y0, s, tmp, tmp2, NULL);
    
    mpz_sub(s, Y2, Y1);
    mpz_sub(tmp, X2, X1);
    mpz_gcdext (d, tmp, tmp2, tmp, N);
    if (mpz_cmp_ui(d, 1) != 0){
        return;
    }
        
    mpz_mul(s, s, tmp);
    mpz_mod(s, s, N); // s = (Y2-Y1)/(X2-X1)

    mpz_add(tmp, X1, X2);
    mpz_powm_ui(X0, s, 2, N);
    mpz_sub(X0, X0, tmp);
    mpz_mod(X0, X0, N); // X = s^2 - (X1+X2)
    
    mpz_sub(Y0, X1, X0);
    mpz_mul(Y0, Y0, s);
    mpz_sub(Y0, Y0, Y1);
    mpz_mod(Y0, Y0, N); // Y = s*(X1-X) - Y1

    mpz_set(X, X0);
    mpz_set(Y, Y0);
    mpz_set_ui(Z, 1);
    mpz_clears(X0, Y0, s, tmp, tmp2, NULL);
    return;
}

void wLADDER(mpz_t d, mpz_t X, mpz_t Y, mpz_t Z, mpz_t m, mpz_t X1, mpz_t Y1, mpz_t Z1, mpz_t N, mpz_t a){
    mpz_t X0, Y0, Z0, mt;
    mpz_inits(mt, X0, Y0, Z0, NULL);
    mpz_set(X0, X1);
    mpz_set(Y0, Y1);
    mpz_set(Z0, Z1);

    unsigned int m_t;

    unsigned int bcnt = mpz_sizeinbase(m, 2);
    
    for(int t = bcnt - 2; t >= 0; t--){
        wDBL(d, X0, Y0, Z0, X0, Y0, Z0, N, a);
        if(mpz_cmp_ui(d, 1) != 0){
            return;
        }
        mpz_fdiv_q_2exp(mt, m, t);
        m_t = mpz_get_ui(mt) & 1;
        if(m_t == 0){
            wADD(d, X0, Y0, Z0, X0, Y0, Z0, X1, Y1, Z1, N, a);    
            if(mpz_cmp_ui(d, 1) != 0){
                return;
            }
        }
    }
    
    mpz_set(X, X0);
    mpz_set(Y, Y0);
    mpz_set(Z, Z0);
    mpz_clears(mt, X0, Y0, Z0, NULL);
    return;
}

void wLADDER_ui(mpz_t d, mpz_t X, mpz_t Y, mpz_t Z, unsigned int m, mpz_t X1, mpz_t Y1, mpz_t Z1, mpz_t N, mpz_t a){
    mpz_t X0, Y0, Z0;
    mpz_inits(X0, Y0, Z0, NULL);
    mpz_set(X0, X1);
    mpz_set(Y0, Y1);
    mpz_set(Z0, Z1);

    unsigned int  m_t = m;
    unsigned int bcnt=0;
    
    while(m_t > 0){
        bcnt++;
        m_t = m_t >> 1;
    }
    
    for(int t = bcnt - 2; t >= 0; t--){
        wDBL(d, X0, Y0, Z0, X0, Y0, Z0, N, a);
        if(mpz_cmp_ui(d, 1) != 0){
            return;
        }
        m_t = (m >> t) & 1;
        if(m_t){
            wADD(d, X0, Y0, Z0, X0, Y0, Z0, X1, Y1, Z1, N, a);    
            if(mpz_cmp_ui(d, 1) != 0){
                return;
            }
        }
    }
    
    mpz_set(X, X0);
    mpz_set(Y, Y0);
    mpz_set(Z, Z0);
    mpz_clears(X0, Y0, Z0, NULL);
    return;
}

/* convert a Montgomery curve to a Weierstrass curve */
// void convertMW(mpz_t x, mpz_t y, mpz_t z, mpz_t X, mpz_t Y mpz_t Z, mpz_t A, mpz_t B){
//     mpz_invert(a, B, N);
//     mpz_mul(a, A, a);
//     mpz_mod(a, a, N); // a = A/B
//     mpz_powm_ui(b, B, 2, N);
//     mpz_invert(b, b, N); // b = 1/B^2
//     mpz_powm_ui(tmp, a, 2, N);
//     mpz_invert(three, 3, N);
//     mpz_mul(tmp, xx, tmp);
//     mpz_sub(b, b, tmp);
//     mpz_powm_ui(tmp, a, 3, N);
//     mpz_set(z, B);
// }

// int main(){
//     mpz_t a, b, x1, y1, z1, x2, y2, z2, p, x0, y0, z0;
//     mpz_inits(a,b,x1,y1,z1,x2,y2,z2,p,x0,y0,z0,NULL);
//     mpz_set_ui(a, 2);
//     mpz_set_ui(b, 1);
//     mpz_set_ui(p, 5);
//     mpz_set_ui(x1, 0);
//     mpz_set_ui(y1, 1);
//     mpz_set_ui(z1, 1);
//     mpz_set_ui(x2, 0);
//     mpz_set_ui(y2, 4);
//     mpz_set_ui(z2, 1);
//     printwCurve(a,b);
//     //wDBL(x0, y0, z0, x1, y1, z1, p, a);
//     wLADDER_ui(x0,y0,z0,7,x1,y1,z1,p,a);
//     gmp_printf(" x0 = %Zd ; y0 = %Zd ; z0 = %Zd \n", x0, y0, z0);
//     mpz_clears(a,b,x1,y1,z1,x2,y2,z2,p,x0,y0,z0,NULL);
// }