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
// Solve the x==r[i] mod q[i] system (0<=i<n)

void mpz_crt(mpz_ptr result, int n, mpz_t *r, mpz_t *q)  {
// Here q[i]!=0, return by -1 if this is not true or the system is unsolvable or n<0
// using the Heindel-Horowitz algorithm

    if(n<0)  {
       mpz_set_si(result,-1);
       return;
    }
    if(n==0)  {
       mpz_set_ui(result,0);
       return;
    }

    int i;
    for(i=0;i<n;i++)
        if(mpz_cmp_ui(q[i],0)==0)  {
           mpz_set_si(result,-1);
           return;
        }
    
    int last=2;
    bool solvable=true;
    
    mpz_t *prodtree,*res,r1,r2,p1,p2,temp,mul,G;
    prodtree=(mpz_t*)(malloc)(2*n*sizeof(mpz_t));
    res=(mpz_t*)(malloc)(2*n*sizeof(mpz_t));
    for(i=0;i<2*n;i++)  mpz_init(prodtree[i]);
    for(i=0;i<2*n;i++)  mpz_init(res[i]);
    mpz_init(r1);
    mpz_init(r2);
    mpz_init(p1);
    mpz_init(p2);
    mpz_init(temp);
    mpz_init(mul);
    mpz_init(G);

    for(i=n;i<2*n;i++)  {
        mpz_abs(prodtree[i],q[i-n]);
        mpz_mod(res[i],r[i-n],prodtree[i]);
    }
    
    for(i=n-1;i>0;i--)  {
       // x==r1 mod p1
       // x==r2 mod p2
       // it is solvable if and only if gcd(p1,p2) divides r2-r1
    
       // G=gcd(p1,p2)
       // x=r1+k*p1
       // k*p1==r2-r1 mod p2
       // k*p1/G==(r2-r1)/G mod p2/G
       // k===(r2-r1)/G*(p1/G)^(-1) mod p2/G

       mpz_set(p1,prodtree[2*i]);
       mpz_clear(prodtree[2*i]);
       mpz_set(p2,prodtree[2*i+1]);
       mpz_clear(prodtree[2*i+1]);
    
       mpz_set(r1,res[2*i]);
       mpz_clear(res[2*i]);
       mpz_set(r2,res[2*i+1]);
       mpz_clear(res[2*i+1]);
        
       mpz_gcd(G,p1,p2);
       mpz_sub(temp,r2,r1);
       if(mpz_divisible_p(temp,G)==0)  { // there is no solution
          solvable=false;
          last=2*i;
          break;
       }
    
       mpz_divexact(p2,p2,G);  // p2=p2/G
       mpz_divexact(temp,temp,G);  // temp=(r2-r1)/G
       mpz_divexact(mul,p1,G); // mul=p1/G
       mpz_invert(mul,mul,p2); // mul=(p1/G)^(-1) mod p2/G
       mpz_mul(temp,temp,mul); // temp=(r2-r1)/G*(p1/G)^(-1)
       mpz_mod(temp,temp,p2);  // temp=(r2-r1)/G*(p1/G)^(-1) mod p2/G
    
       mpz_set(res[i],r1);
       mpz_addmul(res[i],temp,p1);
       mpz_mul(prodtree[i],p1,p2); // product=p1*p2/G=lcm(p1,p2)
    }
    
    
    if(solvable)  mpz_set(result,res[1]);
    else          mpz_set_si(result,-1);

    mpz_clear(r1);
    mpz_clear(r2);
    mpz_clear(p1);
    mpz_clear(p2);
    mpz_clear(temp);
    mpz_clear(mul);
    mpz_clear(G);
    for(i=0;i<last;i++)  {
       mpz_clear(res[i]);
       mpz_clear(prodtree[i]);
    }
    free(res);
    free(prodtree);
    
    return;
}
