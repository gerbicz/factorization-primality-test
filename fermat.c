#include <stdio.h>
#include "gmp.h"


int main()  {

    int i,n;
    
    mpz_t maradek;
    mpz_t q;
    mpz_t r;
    mpz_t N;
    mpz_init(maradek);
    mpz_init(q);
    mpz_init(r);
    mpz_init(N);
    
    printf("Pepin teszt Fermat szamokra (2^(2^n)+1-re).\n");
    printf("Kerem az n kitevot: ");
    scanf("%d",&n);
    if(n<0)  {
       printf("n-nek nemnagativ egesznek kell lennie!\n");
       return 0;
    }
    if(n==0)  {
       printf("Pepin teszt nem alkalmazhato, de 2^(2^0)+1=2 prim.\n");
       return 0;
    }
    // Fermat(n)=2^(2^n)+1
    mpz_ui_pow_ui(N,2,1<<n);
    mpz_add_ui(N,N,1);
    
    // Pepin teszt
    mpz_set_ui(maradek,3);

    for(i=1;i<=(1<<n)-1;i++)  {
        mpz_mul(maradek,maradek,maradek);
              
        // Tenyleges osztas elkerulheto, mert a redukcio megoldhato bitshift-bitmaskkal
        mpz_fdiv_q_2exp(q,maradek,1<<n);  // bitshift
        mpz_fdiv_r_2exp(r,maradek,1<<n);  // bitmask
              
        mpz_sub(maradek,r,q);
        if(mpz_sgn(maradek)==-1)  {
           mpz_add(maradek,maradek,N);
        }
    }
    
    mpz_add_ui(maradek,maradek,1);
    if(mpz_cmp(maradek,N)!=0)  printf("Teszt befejezve, 2^(2^%d)+1 Fermat szam osszetett.\n",n);
    else  printf("Teszt befejezve, 2^(2^%d)+1 Fermat prim!!!!!\n",n);
    
    mpz_clear(maradek);
    mpz_clear(q);
    mpz_clear(r);
    mpz_clear(N);
    
    return 0;
}
