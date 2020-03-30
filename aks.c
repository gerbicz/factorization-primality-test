// Gerbicz Robert
#include <stdio.h>
#include <stdlib.h>
#include <math.h> // log, sqrt fuggvenyekhez
#include <time.h>
#include "gmp.h" // vagy "mpir.h"

// Ha n<2^40000, akkor varhatoan r<2^31 az AKS modszernel
// igy r az int-be is belefer, ezt hasznalni is fogjuk.

double dlog2(mpz_t n)  {// log_2(n) kozelitese, kb. 10^(-9) hibaval
    
    double d=0.0;
    mpz_t temp;
    mpz_init(temp);
    
    mpz_set(temp,n);
    while(mpz_cmp_ui(temp,1073741824)>=0)  {
         d+=1.0;
         mpz_fdiv_q_2exp(temp,temp,1);
    }
    return d+log((double) mpz_get_ui(temp))/log(2);
}

int keres(mpz_t n,double d2)  {// d2=log_2(n)
    // r keresese AKS-hez, melyre o_r(n)>log_2(n)^2
    // itt valojaban a rendet nem szamitjuk ki,
    // kulonben nem lenne ez polinomialis algoritmus
    int f,hatar,v,r,res,mult;
    
    mpz_t temp;
    mpz_init(temp);
    
    hatar=(int) ((double) d2*d2);
    
    r=1;
    while(1)  {
         r++;
         res=1;
         mult=mpz_mod_ui(temp,n,r);
         v=1;
         for(f=1;f<=hatar;f++)  {
             res=((long long int) mult*res)%r;// res=n^f mod r
             if(res==1)  {v=0;break;}
         }
         if(v)  return r;// megtalaltuk az r erteket
    }
    return -1;// ide valojaban nem jutunk el
}

int eulerphi(int n)  {// gyorsasaga nem lenyeges

    int e=1,p;
    
    for(p=2;p*p<=n;p++)
        if(n%p==0)  {
           n/=p;
           e*=p-1;
           while(n%p==0)  n/=p,e*=p-1;
        }
    
    if(n>1)  e*=n-1;
    
    return e;
}

mpz_t *coeff;
void F(mpz_ptr result,int bithossz,int st,int en)  {
    // result=sum(i=st,en,coeff[i]*b^i) kiszamitasa
    // oszd meg es uralkodj elvevel
    // b=2^bithossz
    
    if(en>st)  {
       mpz_set_ui(result,0);
    }
    
    if(st==en)  {
       mpz_mul_2exp(result,coeff[st],st*bithossz);
       return;
    }
    
    int mid=(st+en)>>1;
    F(result,bithossz,st,mid);
    mpz_t temp;
    mpz_init(temp);
    F(temp,bithossz,mid+1,en);
    mpz_add(result,result,temp);
    mpz_clear(temp);
    return;
}

void invF(mpz_t T,int bithossz,int st,int en)  {
    // T=sum(i=st,en,coeff[i]*b^i) osszefuggesbol
    // coeff[i] meghatarozasa oszd meg es uralkodj elvevel
    // b=2^bithossz itt
    if(st>en)  return;
    
    if(st==en)  {
       mpz_fdiv_q_2exp(coeff[st],T,bithossz*st);
       return;
    }
    
    int mid=(st+en)>>1;// kozepso tag indexe
    mpz_t temp;
    mpz_init(temp);
    
    mpz_fdiv_r_2exp(temp,T,bithossz*(mid+1));
    invF(temp,bithossz,st,mid);
    mpz_sub(T,T,temp);
    invF(T,bithossz,mid+1,en);
    mpz_clear(temp);

    return;
}

void egyszerusit(mpz_t n,int r)  {
    // mod (x^r-1,n) redukcio (2*r-1-edfoku a polinomunk az abrazolasban)
    int i;
    for(i=0;i<r;i++)  {
        mpz_add(coeff[i],coeff[i],coeff[i+r]);
        mpz_mod(coeff[i],coeff[i],n);
    }
    return;
}

int gcd(int a,int b)  {// lnko(a,b) eukleideszi algoritmussal
    
    if(b==0)  return 0;
    return gcd(b,a%b);
}

int AKS(mpz_t n)  {// AKS(n)=isprime(n)

    if(mpz_cmp_ui(n,2)<0)  {printf("n>1-nek kell lennie!\n");return -1;}

    if(mpz_perfect_power_p(n)!=0)  return 0; // trivi eset, n teljes hatvany

    int a,G,i,j,L,r,t,hatar,*bit,bithossz,v,res;
    double d2;
    mpz_t T,temp;
    mpz_init(T);
    mpz_init(temp);
    
    d2=dlog2(n);
    r=keres(n,d2); // r az AKS modszerhez
    for(a=1;a<=r;a++)  {
        res=mpz_mod_ui(temp,n,a);
        G=gcd(a,res);
        if(1<G&&mpz_cmp_ui(n,G)>0)  return 0; // trivi eset, n osszetett
    }
    if(mpz_cmp_ui(n,r)<=0)  return 1;// trivi eset, n prim

    bit=(int*)(malloc)(((int) d2+2)*sizeof(int));
    
    // n binaris jegyeinek meghatarozasa
    mpz_set(temp,n);
    L=0;
    while(L==0||mpz_sgn(temp)==1)  {
        bit[L++]=(mpz_odd_p(temp)!=0);
        mpz_fdiv_q_2exp(temp,temp,1);
    }
    
    // b=2^bithossz meghatarozasa, hogy a negyzetre emelesnel
    // ne legyen tulcsordulas
    bithossz=2*((int) d2+1);
    bithossz++;
    t=r;
    while(t)  t/=2,bithossz++;
    
    coeff=(mpz_t*)(malloc)(2*r*sizeof(mpz_t));
    for(i=0;i<2*r;i++)
        mpz_init(coeff[i]);
    
    hatar=(int) ((double) sqrt(eulerphi(r))*d2);
    for(a=1;a<=hatar;a++)  {
        // konstans 1 polinom felepitese
        for(i=0;i<r;i++)
            mpz_set_ui(coeff[i],(i==0));
        
        for(i=L-1;i>=0;i--)  {
            // polinom negyzetre emelese
            F(T,bithossz,0,r-1); // T meghatarozasa a polinom egyutthatoibol
            mpz_mul(T,T,T);      // negyzetre emeles
            invF(T,bithossz,0,2*r-1); // polinom egyutthatoinak visszanyerese
            egyszerusit(n,r);    // mod (x^r-1,n)
            
            if(bit[i])  {// (x+a)-val szorzunk, ez egyszerubb sokkal, mint 
                         // a negyzetreemeles, helyben elvegezheto
               mpz_set(temp,coeff[r-1]);
               for(j=r-1;j>0;j--)  {
                   mpz_mul_ui(coeff[j],coeff[j],a);
                   mpz_add(coeff[j],coeff[j],coeff[j-1]);
                   mpz_mod(coeff[j],coeff[j],n);
               }
               mpz_mul_ui(coeff[0],coeff[0],a);
               mpz_add(coeff[0],coeff[0],temp);
               mpz_mod(coeff[0],coeff[0],n);  
            }
        }
        
        // ellenorzes, hogy a hatvany egyenlo-e (x^n+a)-val mod (x^r-1,n)
        res=mpz_mod_ui(temp,n,r);
        v=1;
        for(i=0;i<r;i++)  {
            if(i==0)  mpz_set_ui(temp,a);
            else if(i==res)  mpz_set_ui(temp,1);
            else mpz_set_ui(temp,0);
            
            mpz_mod(temp,temp,n);            
            if(mpz_cmp(coeff[i],temp)!=0)  {v=0;break;}
        }
        if(!v)  { // n osszetett
           for(i=0;i<2*r;i++)  mpz_clear(coeff[i]);
           free(coeff);
           mpz_clear(temp);
           mpz_clear(T);
           return 0;
        }
    }
    for(i=0;i<2*r;i++)
        mpz_clear(coeff[i]);
    free(coeff);
    mpz_clear(temp);
    mpz_clear(T);
    return 1; // n prim
}

int main()  {
    
    mpz_t n; 
    mpz_init(n);
    time_t sec;
    
    printf("AKS teszt, kerem a szamokat!\n");

    while(gmp_scanf("%Zd",&n)!=EOF)  {
        time_t sec=time(NULL);
        int v=AKS(n);
        if(v==1)  printf("n prim\n");
        else if(v==0)  printf("n osszetett\n");
        printf("%ld sec.\n",time(NULL)-sec);
    }
    mpz_clear(n);
    
    return 0;
}
