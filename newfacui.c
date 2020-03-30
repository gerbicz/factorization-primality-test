#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"

// Computation of n factorial by computing the prime factoriation of n!,
// using iterated squaring and multiplication by a "small" number idea and binary splitting
// written by Robert Gerbicz


/* mpz_fac_ui(result, n) -- Set RESULT to N!.

Copyright 1991, 1993, 1994, 1995, 2000, 2001, 2002, 2003 Free Software
Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA. */


void binary_splitting (mpz_ptr result,unsigned long int *a,unsigned long int L)  {
// mulptiplication by binary splitting
  unsigned long int i,L0,L1;
  mpz_t temp;

  if(L==0)  {
     mpz_set_ui (result, 1);
     return;
  }

  if(L<=3)  {
     mpz_set_ui(result,a[0]);
     for(i=1;i<L;i++)
          mpz_mul_ui(result,result,a[i]);
     return;
  }

  L0=L/2;
  L1=L-L0;
  binary_splitting(result,a,L1);
  mpz_init(temp);
  binary_splitting(temp,a+L1,L0);
  mpz_mul(result,result,temp);
  mpz_clear(temp);
  return;
}


void mpz_fac_ui(mpz_ptr result, unsigned long int n)  {
   
   if(n<2)  {
      mpz_set_ui(result,1);
      return;
   }
   
   unsigned long int Bit[32],e,N,g,p2,*S,*exponent,*isprime,*primes,count,i,p,primepi,sq,n2=(n-1)>>1,n64=(n>>6)+1;
   int h,expo;
   mpz_t temp;
   mpz_init(temp);
   isprime=(unsigned long int*)(malloc)(n64*sizeof(unsigned long int));
   
   Bit[0]=1;
   for(i=1;i<32;i++)  Bit[i]=Bit[i-1]<<1; // Bit[i]=2^i
   
   // determine all odd primes up to n by sieve
   for(i=0;i<n64;i++) isprime[i]=0xffffffff;
   
   sq=(unsigned long int) sqrt(n)+1;
   for(p=3;p<=sq;p+=2)  {
       if(isprime[p>>6]&Bit[(p>>1)&31])  {
          for(i=(p*p-1)>>1;i<=n2;i+=p)  isprime[i>>5]&=~Bit[i&31];
       }
   }
   
   primepi=0;
   for(i=0;i<=n2;i++)
       primepi+=((isprime[i>>5]&Bit[i&31])>0);
   
   primes=(unsigned long int*)(malloc)(primepi*sizeof(unsigned long int));
   S=(unsigned long int*)(malloc)(primepi*sizeof(unsigned long int));
   exponent=(unsigned long int*)(malloc)(primepi*sizeof(unsigned long int));

   primepi=0;
   // 1 is not prime and hasn't cancelled, so start from 1*2+1=3
   for(i=1;i<=n2;i++)
      if((isprime[i>>5]&Bit[i&31])>0)  {
          p=2*i+1;
          N=n;
          e=0;
          while(N)  N/=p,e+=N;
          primes[primepi]=p;      // store prime
          exponent[primepi]=e;    // exponent of p in the factorization of n!      
          primepi++;
      }
          
   free(isprime);
  
   mpz_set_ui(result,1);
   
   expo=0,p2=1;
   N=n;
   while(N)  N>>=1,p2<<=1,expo++;
   for(h=expo;h>=0;h--)  {
       // collect all primes for which in the factorization of n! the primes[g]'s exponent's h-th bit is 1
       count=0;
       // note that exponent[] is a decreasing array
       for(g=0;(g<primepi)&&(exponent[g]>=p2);g++)
           if(((exponent[g]>>h)&1)==1)  {
                S[count]=primes[g];
                count++;
           }
       binary_splitting(temp,S,count); // build the product by binary splitting
       mpz_pow_ui(result,result,2); // squaring
       mpz_mul(result,result,temp); // multiplcation by a not so large number
       p2>>=1;
   }
   
   N=n;
   e=0;
   while(N)  N>>=1,e+=N;
   mpz_mul_2exp(result,result,e);  // shift the number to finally get n!
   
   free(S);
   free(primes);
   free(exponent);
   mpz_clear(temp);
   return;
}
