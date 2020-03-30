#include <stdio.h>
#include "gmp.h"


int prim(int n)  {
// naiv modszer, de futasi ido nem lenyeges
     if(n<2)  return 0;
     int i;
     for(i=2;i*i<=n;i++)
         if(n%i==0)  return 0;
     return 1;
}

int main()  {

    int i,n,res;
    
    mpz_t maradek;
    mpz_t q;
    mpz_t r;
    mpz_t N;
    mpz_t ketto;
    mpz_init(maradek);
    mpz_init(q);
    mpz_init(r);
    mpz_init(N);
    mpz_init(ketto);
    
    printf("Mersenne szamok tesztelese.\n");
    printf("Kerem az n kitevot: ");
    scanf("%d",&n);
    // Eloszor ellenorzes, hogy n az prim-e
    if(prim(n)==0)  printf("Kitevo nem prim, igy a 2^%d-1 sem.\n",n);
    else {
          // primosztoi q=2*k*p+1 alakuak,ahol q==1,7 mod 8 is teljesul
          mpz_set_ui(ketto,2);
          mpz_set_ui(q,1);
          for(i=1;i<=n;i++)  {
              mpz_add_ui(q,q,2*n);
              res=mpz_mod_ui(r,q,8);
              if((res==1)||(res==7))  {
                  mpz_powm_ui(r,ketto,n,q);
                  if((mpz_cmp_ui(r,1)==0)&&(mpz_sizeinbase(q,2)!=n))  {
                     gmp_printf("2%d-1 osszetett, mert egy valodi osztoja: %Zd\n",n,q);
                     return 0;
                  }
              }
          }
          printf("Primosztot nem talaltam, Lucas-Lehmer teszt elinditva.\n");
          // N=2^n-1
          mpz_ui_pow_ui(N,2,n);
          mpz_sub_ui(N,N,1);
          // Lucas tesztben kezdoelem=4, de ezt már mod N veszem
          if(n!=2)  mpz_set_ui(maradek,4);
          else      mpz_set_ui(maradek,0);

          for(i=1;i<=n-2;i++)  {
              mpz_mul(maradek,maradek,maradek);
              
              // Tenyleges osztas elkerulheto, mert a redukcio megoldhato bitshift-bitmaskkal
              mpz_fdiv_q_2exp(q,maradek,n);  // bitshift
              mpz_fdiv_r_2exp(r,maradek,n);  // bitmask
              
              mpz_add(maradek,q,r);
              if(mpz_cmp(maradek,N)>=0)  {
                 mpz_sub(maradek,maradek,N);
              }
              mpz_sub_ui(maradek,maradek,2);
              if(mpz_sgn(maradek)==-1)  {
                 mpz_add(maradek,maradek,N);
              }
          }
          if(mpz_cmp_ui(maradek,0)!=0)  printf("Teszt befejezve, 2^%d-1 osszetett.\n",n);
          else  printf("Teszt befejezve, 2^%d-1 Mersenne prim!!!!!\n",n);
    }
    mpz_clear(maradek);
    mpz_clear(q);
    mpz_clear(r);
    mpz_clear(N);
    mpz_clear(ketto);
    return 0;
}
