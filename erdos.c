#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "gmp.h" // vagy "mpir.h"

// f(n) is the total number of nonnegative integers for that n*n+1<2^e or 10^e
// s(n) is the reciprocial sum

#define E 34  // x^2+1<2^(2*E)-ig szitalunk
#define blokk 17179869184LL  // 64-gyel oszthatonak kell lennie

int rem8[8]={0,0,0,1,0,1,0,0}; // jacobi-hoz

long long int mulmod(long long int a, long long int b, long long int p)  {
    long long int y=(long long int)((double)a*(double)b/p+0.5);
    long long int  r=a*b-y*p;
    if(r<0) r=r+p;

    return r;
}

long long int hatvanymod(long long int alap,long long int kitevo,long long int mod)  {
    long long int maradek=1,h=alap;
         
    while(kitevo)  {
          if(kitevo&1)  maradek=mulmod(maradek,h,mod);
          h=mulmod(h,h,mod);
          kitevo>>=1;
    }
    return maradek;
}

int jacobi(long long int a,long long int b)
{
       if(b<0) b=-b;
       a%=b;
       if(a<0) a+=b;
       if(a==0)  return 0;
       int sym=0;
       long long int c;
       while(a>1) {
               if((a&3)==0)  a>>=2;
               else if((a&3)==2) {
                       sym^=rem8[b&7],a>>=1;
                       if((a&b&3)==3) sym^=1;
                               c=a,a=b,b=c;
                               if(a>=b)  {
                                       a-=b;
                                       if(a>=b)  {
                                          a-=b;
                                            if(a>=b)  {
                                                    a-=b;
                                                    if(a>=b)  {
                                                            a-=b;
                                                            if(a>=b)  {
                                                                    a-=b;
                                                                    if(a>=b)  a%=b;
                                                                       }
                                                              }
                                                      }
                                              }
                                      }
                    }
       else{
             if((a&b&3)==3) sym^=1;
             c=a,a=b,b=c;
             if(a>=b)  {
                     a-=b;
                     if(a>=b)  {
                         a-=b;
                         if(a>=b)  {
                                 a-=b;
                                 if(a>=b)  {
                                         a-=b;
                                         if(a>=b)  {
                                                 a-=b;
                                                 if(a>=b)  {
                                                         a-=b;
                                                         if(a>=b)  {
                                                                 a-=b;
                                                                 if(a>=b)  {
                                                                         a-=b;
                                                                         if(a>=b)  {
                                                                                 a-=b;
                                                                                 if(a>=b)  a%=b;
                                                                                    }
                                                                           }
                                                                   }
                                                           }
                                                   }
                                           }
                                   }
                           }
                   }
               }
     }

    return 1-(sym<<1);
}

int kvadratikusnemmaradekkeres(long long int p)  {
    int a=2;
    while(jacobi(a,p)!=-1)  a++;
    // jacobi(a,p)=-1
    return a;
}

int main()  {
    
    int lo,hi,seconds=time(NULL);
    unsigned int *A,*isprime,Bit[32],g,e2,e10;
    long long int *prim,*primtomb,c,L,U,i,j,k,a,s[2],n,ct,temp;
    long long int st,en,p,primepisq,N,SQ,db,hatar2,hatar10,blokk2,blokk64;
    double p2,p10;
    
    mpz_t gmp_n;
    mpf_t tot,rec;
    
    mpf_set_default_prec(128);
    mpz_init(gmp_n);
    mpf_init(tot);
    mpf_init(rec);
    
    mpf_set_d(tot,0.5);

    
    N=1LL<<E;
    SQ=(long long int) sqrt(N)+1;
    blokk2=blokk/2;
    blokk64=blokk/64;
    e2=4;
    e10=1;
    p2=4;
    p10=sqrt(10);
    hatar2=4;
    hatar10=4;
    db=0;
    
    A=(unsigned int*)(malloc)(blokk64*sizeof(unsigned int));
    isprime=(unsigned int*)(malloc)(SQ*sizeof(unsigned int));
    primtomb=(long long int*)(malloc)(SQ*sizeof(long long int));
    prim=(long long int*)(malloc)(SQ*sizeof(long long int));
    
    Bit[0]=1;
    for(n=1;n<32;n++)  Bit[n]=2*Bit[n-1];

    for(n=0;n<SQ;n++)  isprime[n]=(n>1);
    for(n=2;n*n<SQ;n++)  {
        if(isprime[n])  {
           for(i=n*n;i<SQ;i+=n)  isprime[i]=0;
        }
    }
    primepisq=0;
    for(n=0;n<SQ;n++)
        if(isprime[n])  prim[primepisq]=n,primepisq++;

    for(st=0;st<N;st+=blokk)  {
    en=st+blokk;
    for(n=0;n<blokk64;n++)  A[n]=0xffffffff;
    if(st==0)  A[0]&=~Bit[0],db++;

    for(L=0;L<en;L+=SQ)  {
        if(L==0)  {
           for(i=0;i<primepisq;i++)  primtomb[i]=prim[i];
           ct=primepisq;
        }
        else {
           U=L+SQ;
           for(j=0;j<SQ;j++)  isprime[j]=1;
           for(i=0;(i<primepisq)&&(prim[i]*prim[i]<=U);i++)  {
                p=prim[i];
                for(j=((p+L-1)/p)*p-L;j<SQ;j+=p)  isprime[j]=0;
           }
           ct=0;
           for(j=0;j<SQ;j++)
               if(isprime[j])  primtomb[ct]=L+j,ct++;
        }
        for(i=0;i<ct;i++)  {
            p=primtomb[i];
            if(p%4==1)  {
               a=kvadratikusnemmaradekkeres(p);
               s[0]=hatvanymod(a,(p-1)/4,p);
               s[1]=p-s[0];
               if(s[0]%2==1)  s[0]+=p;
               if(s[1]%2==1)  s[1]+=p;                  
               for(k=0;k<2;k++)  {
                   c=((st-s[k]+2*p-1)/(2*p))*2*p+s[k];
                   if(((p-1)/c==c)&&(p==c*c+1))  c+=2*p;
                   for(j=(c-st)/2;j<blokk2;j+=p)   A[j>>5]&=~Bit[j&31];
               }
            }
        }
   }

   for(i=0;i<blokk64;i++)  {
       g=A[i];
       n=st+64*i;
       for(j=0;j<32;j++)  {
           if(g&1)  {
              db++;
              temp=n;
              lo=n&((1<<30)-1);
              temp>>=30;
              hi=temp;
              mpz_set_ui(gmp_n,hi);
              mpz_mul_2exp(gmp_n,gmp_n,30);
              mpz_add_ui(gmp_n,gmp_n,lo);
              
              mpz_mul(gmp_n,gmp_n,gmp_n);
              mpz_add_ui(gmp_n,gmp_n,1);
              
              mpf_set_z(rec,gmp_n);
              mpf_ui_div(rec,1,rec);
              mpf_add(tot,tot,rec);
           }
           g>>=1;         
           n+=2;
           if(n>=hatar2)  {
              gmp_printf("f(2%c%d)=%lld,s(2%c%d)=%.*Ff\n",'^',e2,db,'^',e2,30,tot);
              e2++;
              p2*=sqrt(2);
              if(e2&1) hatar2=(long long int) (p2+1.0);
              else     hatar2=(long long int) (p2+0.4);
           }
           if(n>=hatar10)  {
              gmp_printf("f(10%c%d)=%lld,s(10%c%d)=%.*Ff\n",'^',e10,db,'^',e10,30,tot);
              e10++;
              p10*=sqrt(10);
              if(e10&1)  hatar10=(long long int) (p10+1.0);
              else       hatar10=(long long int) (p10+0.4);
           }
       }
   }
   }
   free(isprime);
   free(primtomb);
   free(A);
   printf("ido=%ld masodperc.\n",time(NULL)-seconds);
   return 0;
}
