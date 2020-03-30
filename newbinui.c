#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"

/* This file is part of the GNU MP Library.

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

// Major speedup for large input by binary splitting.
// Written by Robert Gerbicz from Hungary.

void binary_splitting (mpz_ptr result,mpz_t *a,unsigned long int L)  {
  unsigned long int i,L0,L1;
  mpz_t temp;

  if(L==0)  {
     mpz_set_ui (result, 1);
     return;
  }

  if(L<=3)  {
     mpz_set(result,a[0]);
     for(i=1;i<L;i++)
          mpz_mul(result,result,a[i]);
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

void mpz_bin_ui(mpz_ptr result,mpz_t n, unsigned long int k)  {
   signed int sg=1; 

   if(mpz_sgn(n)==-1)  {
      // property:
      // binomial(n,k)=(-1)^k*binomial(k-n-1,k)
      mpz_add_ui(n,n,1);
      mpz_neg(n,n);
      mpz_add_ui(n,n,k);
      if(k&1)  sg=-1;
   }

   if(mpz_cmp_ui(n,k)<0)   {
      mpz_set_ui(result,0);
      return;
   }
   
   // binomial(n,k)=binomial(n,n-k)
   // set k=min(k,n-k)
   mpz_t k2;
   mpz_init(k2);
   mpz_sub_ui(k2,n,k);
   if(mpz_cmp_ui(k2,k)<0)  k=mpz_get_ui(k2);

   // two trivial case:
   if(k==0)  {
      mpz_set_si(result,sg);
      mpz_clear(k2);
      return;
   }

   if(k==1)  {
      mpz_set(result,n);
      if(sg==-1)  mpz_neg(result,result);
      mpz_clear(k2);
      return;
   }

   unsigned long int i,j,m,p,q,sq,st,*A,*smallprimes;
   unsigned long int primepi_sq,pos,count,num,prod,mul,up; 
   sq=(unsigned long int) sqrt(k)+1;

   mpz_t res;
   mpz_init(res);
   mpz_t *S;
   S=(mpz_t*)(malloc)(k*sizeof(mpz_t));
   smallprimes=(unsigned long int*)(malloc)(sq*sizeof(unsigned long int));
   A=(unsigned long int*)(malloc)(sq*sizeof(unsigned long int));
   
   // determine primes up to sqrt(k) by sieve
   for(i=0;i<sq;i++)  A[i]=(i>1);
   for(i=0;i*i<sq;i++)  {
       if(A[i])  {
          for(j=i*i;j<sq;j+=i)  A[j]=0;
       }
   }
   primepi_sq=0;
   for(i=0;i<sq;i++)
       if(A[i])  smallprimes[primepi_sq]=i,primepi_sq++;
   
   for(i=0;i<k;i++)  mpz_init(S[i]);
   for(i=0;i<k;i++)  mpz_set_ui(S[i],1);
   
   // after the multiplication process binomial(n,k)=prod(i=0,k-1,(n-i)/S[i])
   // (multiplication is faster than division), that is the reason why we
   // not choose division and get binomial(n,k)=prod(i=0,k-1,S[i])
   // (obviously we need one division (n-i)/S[i] per each element in S array)
   for(i=0;;i+=sq)  {
       for(j=0;j<sq;j++)  A[j]=1;
       if(i==0)  A[0]=0,A[1]=0;
       
       if(i>k-sq)  up=k;
       else        up=i+sq-1;
       // sieve the primes in [i,up] interval
       for(j=0;(j<primepi_sq)&&(smallprimes[j]*smallprimes[j]<=up);j++)  {
            p=smallprimes[j];
            st=p-(i%p);
            if(st==p)  st=0;
            if((i<=p)&&(i+st<=p))  st=2*p-i;
            for(m=st;m<sq;m+=p)    A[m]=0;
       }
       // sieve done
       
       // sieve out the factors in the S array
       for(j=0;(j<=up-i);j++)  {
           if(A[j])  {
              p=i+j;
              q=1;
              while(q<=k/p)  {
                    q*=p;
                    pos=mpz_fdiv_r_ui(res,n,q); // n-pos is divisible by q
                    count=k/q;  // number of numbers up to k divisible by q
                    for(m=0;m<count;m++)  {
                        mpz_mul_ui(S[pos],S[pos],p);
                        pos+=q;
                    }
              }
          }
      }
      if(i>k-sq)  break;
  }
  free(A);
  free(smallprimes);
  
  num=0;
  for(i=0;i<k;i++)  {
      if(mpz_cmp(S[i],n)!=0)  {
         mpz_divexact(S[num],n,S[i]); // this is an exact division
         num++;
      }
      mpz_sub_ui(n,n,1);
  }

  binary_splitting(result,S,num);
  if(sg==-1)  mpz_neg(result,result);

  mpz_clear(k2);
  mpz_clear(res);
  for(i=0;i<k;i++)  mpz_clear(S[i]);
  free(S);
  return;
}
