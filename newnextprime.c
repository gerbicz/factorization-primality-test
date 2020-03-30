/*

This file is part of the MPIR Library.

The MPIR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 2.1 of the License, or (at
your option) any later version.

The MPIR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the MPIR Library; see the file COPYING.LIB.  If not, write
to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
Boston, MA 02110-1301, USA.
*/


#include <stdio.h>
#include <stdlib.h>
#include "mpir.h" // or use "gmp.h"

// Written by Robert Gerbicz
// fast next (pseudo)prime code, it is using product/remainder tree to sieve up to about log2(n)^2

void mpz_nextprime(mpz_ptr result, mpz_srcptr n)  {
     
     if(mpz_cmp_ui(n,30)<=0)  {// handle the small cases
        if(mpz_cmp_ui(n,2)<0)  {
           mpz_set_ui(result,2);
           return;
        }
        unsigned int N=mpz_get_ui(n),i;
        bool isprime=false;
        while(!isprime)  {
              isprime=true;
              N++;
              for(i=2;i*i<=N;i++)
                  if(N%i==0)  {
                     isprime=false;
                     break;
                  }
        }
        mpz_set_ui(result,N);        
        return;
     }
     
     unsigned int Bits[32],*isprime,*prm,*test;
     unsigned long i,i2,j,k,p,res,st,sq,primepi,ppi,primepisq,interval,interval2,interval64,as,as2,as64;
     bool foundpseudoprime=false;
     
     mpz_t B,B2,SQ,RES,ST,P,two,N;
     mpz_t *tree;
     mpz_init(B);
     mpz_init(B2);
     mpz_init(SQ);
     mpz_init(RES);
     mpz_init(ST);
     mpz_init(P);
     mpz_init(two);
     mpz_init(N);
     
     interval=mpz_sizeinbase(n,2);
     if(interval&1)  interval++; // now interval is about log2(n) and even number
     interval2=interval>>1;
     interval64=interval/64+1;
     isprime=(unsigned int*)(malloc)(interval64*sizeof(unsigned int));

     as=8*interval; // at once we sieve 8*interval numbers, so about 11.5*log(n) numbers
     as2=as>>1;
     as64=as/64+1;
     test=(unsigned int*)(malloc)(as64*sizeof(unsigned int));
     
     Bits[0]=1;
     for(i=1;i<32;i++)  Bits[i]=Bits[i-1]<<1;
     
     // sieve primes in [0,interval]
     for(i=0;i<interval64;i++)  isprime[i]=0xffffffff;
     isprime[0]&=~Bits[0];
     for(i=3;i*i<interval;i+=2)  {
         i2=(i-1)>>1;
         if(isprime[i2>>5]&Bits[i2&31])  {
            for(j=(i*i-1)>>1;j<interval2;j+=i)  isprime[j>>5]&=~Bits[j&31];
         }
     }
     
     primepisq=1; // this will be the number of primes in [0,interval)
     for(i=0;i<interval2;i++)
         if(isprime[i>>5]&Bits[i&31])  primepisq++;
     prm=(unsigned int*)(malloc)(primepisq*sizeof(unsigned int));// to store the primes in [0,interval)
     prm[0]=2;
     primepisq=1;
     for(i=0;i<interval2;i++)
         if(isprime[i>>5]&Bits[i&31])  {
            prm[primepisq]=2*i+1;
            primepisq++;
         }

     mpz_add_ui(result,n,1);
     if(mpz_even_p(result)!=0)
        mpz_add_ui(result,result,1);// the prime should be odd

     while(!foundpseudoprime)  {

     for(k=0;k<as64;k++)  test[k]=0xffffffff;
     mpz_set_ui(B,0);
     for(i=0;i<interval;i++)  {
        mpz_add_ui(B2,B,interval);
        mpz_sqrt(SQ,B2);
        sq=mpz_get_ui(SQ); // sq=sqrt(B2)
         
        // find odd primes in [B,B2)
        for(k=0;k<interval64;k++)  isprime[k]=0xffffffff;
        if(i==0)  isprime[0]&=~Bits[0];
        for(j=1;(j<primepisq)&&(prm[j]<=sq);j++)  {
            p=prm[j];
            if(i==0)  st=p*p;
            else  {
               res=mpz_mod_ui(RES,B,p);
               st=p-res;
               if((st&1)==0)  st+=p;
            }
            for(k=(st-1)>>1;k<interval2;k+=p)  isprime[k>>5]&=~Bits[k&31];
        }
        // primepi=number of odd primes in [B,B2)
        primepi=0;
        for(k=0;k<interval2;k++)
            if(isprime[k>>5]&Bits[k&31])  primepi++;
        
        // build up the product/remainder tree
        tree=(mpz_t*)(malloc)(2*primepi*sizeof(mpz_t));
        for(k=0;k<2*primepi;k++)
            mpz_init(tree[k]);
        
        ppi=0;
        for(k=0;k<interval2;k++)
            if(isprime[k>>5]&Bits[k&31])  {
               mpz_add_ui(tree[primepi+ppi],B,2*k+1);
               ppi++;
            }
        
        for(k=primepi-1;k>0;k--)
            mpz_mul(tree[k],tree[2*k],tree[2*k+1]);

        mpz_mod(tree[1],result,tree[1]);
        for(k=2;k<2*primepi;k++)
            mpz_mod(tree[k],tree[k>>1],tree[k]);
        
        // sieve by primes from [B,B2)
        ppi=primepi;
        primepi=0;
        for(k=0;k<interval2;k++)
            if(isprime[k>>5]&Bits[k&31])  {
               mpz_add_ui(P,B,2*k+1);// sieve by P (P is prime)
               mpz_sub(ST,P,tree[ppi+primepi]);
               if(mpz_cmp(ST,P)==0)  mpz_set_ui(ST,0);
               if(mpz_odd_p(ST)!=0)  mpz_add(ST,ST,P);
               mpz_divexact_ui(ST,ST,2);

               while(mpz_cmp_ui(ST,as2)<0)  {
                     st=mpz_get_ui(ST);
                     test[st>>5]&=~Bits[st&31];  
                     mpz_add(ST,ST,P);
               }
               primepi++;
            }
        
        for(k=0;k<2*primepi;k++)
            mpz_clear(tree[k]);
        free(tree);
        
        mpz_add_ui(B,B,interval); // B->B2
   }
       for(k=0;k<as2;k++)
            if(test[k>>5]&Bits[k&31])  {
               mpz_add_ui(N,result,2*k);
               // N is not divisible by primes from [0,interval^2), so primes up to about log2(n)^2
               // check if N is a Fermat pseudoprime for base=2 or not
               mpz_set_ui(two,2);
               mpz_powm(N,two,N,N);
               if(mpz_cmp_ui(N,2)==0)  {
                  // found pseudoprime
                  foundpseudoprime=true;
                  mpz_add_ui(result,result,2*k);
                  break;              
               }
            }
       if(foundpseudoprime) break;
       mpz_add_ui(result,result,as);
   }
   mpz_clear(B);
   mpz_clear(B2);
   mpz_clear(SQ);
   mpz_clear(RES);
   mpz_clear(ST);
   mpz_clear(P);
   mpz_clear(two);
   mpz_clear(N);

   free(isprime);
   free(prm);
   free(test);
   return;
}
