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

// for reference see: http://cr.yp.to/lineartime/multapps-20080515.pdf
// section 23, page 370-372


void mpz_crt_coprime(mpz_ptr result, int n, mpz_t *r, mpz_t *q)  {
// Here q[i]!=0 are pairwise coprime integers, return by -1 if this is not true or n<0
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
     
     int last;
     bool coprime=true;
     
     mpz_t temp,*prodtree,*remaindertree,*T;
     prodtree=(mpz_t*)(malloc)(2*n*sizeof(mpz_t));
     remaindertree=(mpz_t*)(malloc)(2*n*sizeof(mpz_t));
     T=(mpz_t*)(malloc)(2*n*sizeof(mpz_t));
     for(i=0;i<2*n;i++)  mpz_init(prodtree[i]);
     for(i=0;i<2*n;i++)  mpz_init(remaindertree[i]);
     for(i=0;i<2*n;i++)  mpz_init(T[i]);
     mpz_init(temp);
     
     // build up the product tree
     for(i=n;i<2*n;i++)
         mpz_abs(prodtree[i],q[i-n]);

     for(i=n-1;i>0;i--)
         mpz_mul(prodtree[i],prodtree[2*i],prodtree[2*i+1]);
     
     mpz_clear(remaindertree[0]);
     // build up the remainder tree
     mpz_set(remaindertree[1],prodtree[1]);
     for(i=2;i<2*n;i++)  {
         mpz_mul(temp,prodtree[i],prodtree[i]);
         mpz_mod(remaindertree[i],remaindertree[i>>1],temp);
         if(i&1)  mpz_clear(remaindertree[i>>1]);
     }
     
     // this is an exact division
     for(i=n;i<2*n;i++)
         mpz_divexact(remaindertree[i],remaindertree[i],q[i-n]);

     for(i=n;i<2*n;i++)  {
         mpz_invert(T[i],remaindertree[i],q[i-n]);
         mpz_clear(remaindertree[i]);
         if(mpz_cmp_ui(T[i],0)==0)  {
            coprime=false;
            mpz_set_si(result,-1);
            last=i;
            break;
         }
         mpz_mul(T[i],T[i],r[i-n]);
         mpz_mod(T[i],T[i],prodtree[i]);  // prodtree[i]=abs(q[i-n])
     }
     if(coprime)  {
         for(i=n-1;i>0;i--)  {
            mpz_mul(T[i],T[2*i],prodtree[2*i+1]);
            mpz_clear(T[2*i]);
            mpz_clear(prodtree[2*i+1]);
            
            mpz_addmul(T[i],T[2*i+1],prodtree[2*i]);
            mpz_clear(T[2*i+1]);
            mpz_clear(prodtree[2*i]);
         }
         mpz_mod(result,T[1],prodtree[1]);
         
         for(i=0;i<2;i++)  {
            mpz_clear(T[i]);
            mpz_clear(prodtree[i]);
         }
     }
     else  {
         for(i=0;i<2*n;i++)      mpz_clear(prodtree[i]);
         for(i=0;i<2*n;i++)      mpz_clear(T[i]);
         for(i=last+1;i<2*n;i++) mpz_clear(remaindertree[i]);
     }

     mpz_clear(temp);
     
     free(prodtree);
     free(remaindertree);
     free(T);
     return;
}
