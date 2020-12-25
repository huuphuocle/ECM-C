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