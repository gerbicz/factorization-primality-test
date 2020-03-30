#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "gmp.h"


int main()  {
    
    unsigned int db,tovabb,B1,B2,I,i,j,L,p,q,elso,felso,elozo,isprime[65536],prime[6542],pr[6542];// pi(2^16)=6542
    time_t sec;
    mpz_t n,G,a,a1,szorzat,tomb[293];  // nagy prímvariánshoz kell a tömb
    mpz_init(n);
    mpz_init(G);
    mpz_init(a);
    mpz_init(a1);
    mpz_init(szorzat);
    
    for(i=0;i<293;i++)  mpz_init(tomb[i]);
    
    // 2^16-ig a prímek generálása Erathoszteneszi szitával
    for(i=0;i<65536;i++)  isprime[i]=(i>1);
    for(i=0;i<256;i++)  {
        if(isprime[i])  {
           for(j=i*i;j<65536;j+=i)  isprime[j]=0;
        }
    }
    db=0;
    for(i=0;i<65536;i++)
        if(isprime[i])  prime[db]=i,db++;

    printf("Pollard p-1 faktorizacios modszere tetszolegesen nagy szamokra\n");
    printf("Nagy primvarianssal.\n");
    printf("Lnkot csak az elso/masodik lepes utan szamol a program!\n");
    printf("p (prim)osztot talal a Pollard, ha p-1=Q*R, ahol Q minden\n");
    printf("primhatvanyosztoja <=B1, tovabba R=1 vagy R prim es (B1<) R<=B2\n");
    printf("Ha lnko=n lenne, akkor probalj kisebb B1/B2 ertekkel!\n");
    printf("B2<=B1 eseten masodik lepest mar nem szamol a program.\n");
    printf("Korlatok: 2<=B1,B2<2^31\n\n");
    printf("Kerem a szamot: (0-ra kilep a program)\n");
    
    while(gmp_scanf("%Zd",&n)!=EOF)  {
           if(mpz_sgn(n)==0)  return 0;
           
           mpz_abs(n,n);
           
           printf("B1 korlat: ");
           scanf("%u",&B1);
           printf("B2 korlat: ");
           scanf("%u",&B2);
           assert(B1<=2147483647&&B2<=2147483647);
           
           sec=time(NULL);
           
           mpz_set_ui(a,2);  // az alap a Pollard módszerben, nem lényeges
           for(I=0;I<=B1;I+=65536)  {
               felso=I+65536;
               if(felso>B1)  felso=B1+1;
               L=felso-I;
               // szita a következõ prímek megtalálására
               for(i=0;i<L;i++)  isprime[i]=1;
               if(I==0)  isprime[0]=0,isprime[1]=0;
               for(i=0;(i<6542)&&(prime[i]*prime[i]<=felso);i++)  {
                   p=prime[i];
                   elso=((I+p-1)/p)*p;
                   if(elso<=p)  elso=2*p;  // kis primeket ne huzzuk ki
                   for(j=elso-I;j<L;j+=p)  isprime[j]=0;
               }
               db=0;
               for(i=0;i<L;i++)
                   if(isprime[i])  {
                      p=I+i;
                      q=p;
                      while(q<=B1/p)  q*=p;  // B1-ig hatvanyozzuk p-t
                      pr[db]=q;
                      db++;
                   }
               // a Pollard módszer lelke
               for(i=0;i<db;i++)  mpz_powm_ui(a,a,pr[i],n);
           }
           
           mpz_sub_ui(a1,a,1);
           mpz_gcd(G,a1,n);
           if(mpz_cmp_ui(G,1)!=0&&mpz_cmp(G,n)!=0)  {
               printf("Nemtrivialis osztot talalt a p-1 modszer az elso lepesben:\n");
               gmp_printf("n osztoja: %Zd\n",G);
               tovabb=0;
           }
           else  if(mpz_cmp(G,n)==0)  {
               printf("lnko=n, probalj kisebb B1-el!\n");
               tovabb=0;
           }
           else  {
               printf("lnko=1 az elso lepesben.\n");
               tovabb=1;
           }
           if(tovabb&&(B2>B1))  {
               
               printf("Nagy prim varians inditasa.\n");
               
               // tomb[i]=a^i mod n
               mpz_set_ui(tomb[0],1);
               for(i=1;i<=292;i++)  // maximális gap=292 a legfeljebb 31 bites prímekre
                   {mpz_mul(tomb[i],tomb[i-1],a);mpz_mod(tomb[i],tomb[i],n);}
               mpz_powm_ui(a,a,B1,n);
               elozo=B1;
               mpz_set_ui(szorzat,1);
               
           for(I=B1+1;I<=B2;I+=65536)  {
               felso=I+65536;
               if(felso>B2)  felso=B2+1;
               L=felso-I;
               // szita a következõ prímek megtalálására
               for(i=0;i<L;i++)  isprime[i]=1;
               for(i=0;(i<6542)&&(prime[i]*prime[i]<=felso);i++)  {
                    p=prime[i];
                    elso=((I+p-1)/p)*p;
                    if(elso<=p)  elso=2*p;
                    for(j=elso-I;j<L;j+=p)  isprime[j]=0;
               }
               db=0;
               for(i=0;i<L;i++)
                   if(isprime[i])  pr[db]=I+i,db++;
               
               // Pollard módszer második lépése
               for(i=0;i<db;i++)  {
                   mpz_mul(a,a,tomb[pr[i]-elozo]);
                   mpz_mod(a,a,n);
                   
                   mpz_sub_ui(a1,a,1);
                   mpz_mul(szorzat,szorzat,a1);
                   mpz_mod(szorzat,szorzat,n);

                   elozo=pr[i];
               }
           }
           mpz_gcd(G,szorzat,n);
           if(mpz_cmp_ui(G,1)!=0&&mpz_cmp(G,n)!=0)  {
               printf("Nemtrivialis osztot talalt a p-1 modszer a masodik lepesben:\n");
               gmp_printf("n osztoja: %Zd\n",G);
           }
           else  if(mpz_cmp(G,n)==0)  {
               printf("lnko=n, probalj kisebb B2-el!\n");
           }
           else  {
               printf("lnko=1 a masodik lepesben is.\n");
           }
        }
      printf("ido=%ld sec.\n\n",time(NULL)-sec);
      printf("Kerem a szamot: (0-ra kilep a program)\n");
   }
   
   mpz_clear(n);
   mpz_clear(G);
   mpz_clear(a);
   mpz_clear(a1);
   mpz_clear(szorzat);
   for(i=0;i<293;i++)  mpz_clear(tomb[i]);
   
   return 0;
}
