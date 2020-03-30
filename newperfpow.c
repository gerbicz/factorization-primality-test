/* mpz_perfect_power_p(arg) -- Return non-zero if ARG is a perfect power,
   zero otherwise.

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

// Major speedup by using power residues, remainder tree, strong Lehmer test
// written by Robert Gerbicz

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"


// for trial division up to 2^10=1024

static const unsigned short primes[] =
{  2,  3,  5,  7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
  59, 61, 67, 71, 73, 79, 83, 89, 97,101,103,107,109,113,127,131,
 137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,
 227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,
 313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,
 419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,
 509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,
 617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,
 727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,
 829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,
 947,953,967,971,977,983,991,997,1009,1013,1019,1021,0};


int single_strong_lehmer_test(mpz_t n,unsigned int p)  {
    // Here n-1 is divisible by p, where p>2 is prime and 1<n<p*p
    // we do only one test ( so trying one base ), so it is possible that we reject a prime
    // but in every case if it says n is prime then n is really prime
    mpz_t a;
    mpz_t a2;
    mpz_t e;
    mpz_init(a);
    mpz_init(a2);
    mpz_init(e);

    int answer;
    int res=mpz_mod_ui(a,n,15015);
    if((res%3==0)||(res%5==0)||(res%7==0)||(res%11==0)||(res%13==0))  {
        answer=0;
    }
    else  {
    mpz_sub_ui(e,n,1);
    mpz_divexact_ui(e,e,p);
           
    mpz_set_ui(a,2);
    mpz_powm(a,a,e,n);  // now a=2^((n-1)/p)  mod n
    mpz_powm_ui(a2,a,p,n);  // now a2=2^(n-1)  mod n
    if(mpz_cmp_ui(a2,1)!=0)  answer=0;
    else  {
         mpz_sub_ui(a,a,1);
         mpz_gcd(a,a,n);
         answer=(mpz_cmp_ui(a,1)==0);
    }}
    mpz_clear(a);
    mpz_clear(a2);
    mpz_clear(e);
    return answer;
}

static unsigned int
gcd (unsigned int a, unsigned int b)
{
 if(a==0)  return b;
 if(b==0)  return a;

 unsigned int c;

 while(b>0)  {
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
             }}}}}}}}
    c=a,a=b,b=c;
 }
 return a;
}



static int
isprime (unsigned int t)
{
  unsigned int q, r, d;

  if (t < 3 || (t & 1) == 0)
    return t == 2;

  for (d = 3, r = 1; r != 0; d += 2)
    {
      q = t / d;
      r = t - q * d;
      if (q < d)
	return 1;
    }
  return 0;
}

int
mpz_perfect_power_p (mpz_srcptr n)  {
// returns by non zero if n is a perfect power, otherwise by 0    
    unsigned long int b,c,i,numbits,test,seconds,found,p,bound_for_p,rem,exponent,special;
    unsigned long int primepi,count,ct,rootres,iq,prod,pos,ppi,exact,n2,*prime;
    unsigned long int answer=0,foundanswer,i1,i2;
    int sgn;

    mpz_t *pprime; // array for the 2kp+1 primes
    mpz_t *tree;  // for remainder tree
    mpz_t u;
    mpz_t p2;
    mpz_t q;
    mpz_t root;
    mpz_init(u);
    mpz_init(p2);
    mpz_init(q);
    mpz_init(root);
    
          if(mpz_cmpabs_ui(n,1)<=0)  {
             answer=1;  // 1,0,-1 are pefect powers
             goto smalldel;
          }
          
          sgn=mpz_sgn(n);
          mpz_abs(u,n);

          n2=mpz_scan1(u,0);
          mpz_tdiv_q_2exp(u,u,n2);

          if((sgn==-1)&&(n2>0))  {
             while((n2&1)==0)  n2>>=1;
          }
          if(n2==1)  {
             answer=0;			
             goto smalldel;
          }

          if(mpz_cmp_ui(u,1)==0)  {  // n is power of two
             answer=1;
             goto smalldel;
           }
          
          special=0;
          if(isprime(n2))  special=1;
          // special=1 if there is remaining only one exponent that we need to test by n-th root test algorithm
          // otherwise it's value is 0
          
          numbits=mpz_sizeinbase(u,2);
          for(i=1;(special==0)&&(primes[i]!=0);i++)   {
              p=primes[i];
              if(mpz_cmp_ui(u,p*p)<0)  {
                 answer=0;
                 goto smalldel;
              }
              
              if(mpz_divisible_ui_p(u,p)!=0)	// divisible by this prime?
              {
	               rem=mpz_tdiv_q_ui(q,u,p*p);
	               if(rem!=0)  {
	                  answer=0;  // prime divides exactly once, reject
	                  goto smalldel;
	               }
	               mpz_swap(q,u);
	               
	               mpz_set_ui(q,p);
	               exponent=2+mpz_remove(u,u,q);
	               
	               n2=gcd(n2,exponent);
	               if(sgn==-1)  {
	                  while((n2&1)==0)  n2>>=1;
	               }
	               
	               if(n2==1)   {
	                  answer=0;  // can't be a perfect power
	                  goto smalldel;
	               }
	               if(mpz_cmp_ui(u,1)==0)  {
	                  answer=1;  // fully factored
	                  goto smalldel;
	               }
	               if(isprime(n2))  special=1;
	               
	               numbits=mpz_sizeinbase(u,2);
	            }
	       }
         
	       if(special)  {
            if(n2==2)  {
               answer=mpz_perfect_square_p(u);
               goto smalldel;
            }
            answer=mpz_root(root,u,n2);
            goto smalldel;
         }

         if(10*2>numbits)  {
            answer=0;  // n isn't factorised fully, and n can't be a perfect power
            goto smalldel;
         }
         if(n2!=0)  {
            // factorize n2 and test its prime divisors as possible exponents
            for(p=2;p*p<=n2;p++)  {
                if(n2%p==0)  {
                   while(n2%p==0)  n2/=p;
                   if(p==2)  {
                      if((sgn==1)&&(mpz_perfect_square_p(u)!=0))  {
                         answer=1;
                         goto smalldel;
                      }
                   }
                   else  {
                      if(mpz_root(root,u,p)!=0)  {
                         answer=1;
                         goto smalldel;
                      }
                   }
                }
            }
            if(n2>1)  { // here n2 is prime
               answer=mpz_root(root,u,n2);
               goto smalldel;
            }
            answer=0;
            goto smalldel;
         }

         // we know no factor of n
          bound_for_p=numbits/10;  // previously we sieved up to 2^10=1024
          // primepi is an upper bound for number of primepowers up to bound_for_p
          primepi=(unsigned long int) ((double) 2.0*bound_for_p/log(bound_for_p));
          
          pprime=(mpz_t*)(malloc)(primepi*sizeof(mpz_t));  // for storing the 2kp+1 primes 
          for(i=0;i<primepi;i++)  mpz_init(pprime[i]);
          prime=(unsigned long int*)(malloc)(primepi*sizeof(unsigned long int));
          
          count=0;
          ppi=0;  // it will be the number of different odd primes up to bound_for_p
          for(p=3;p<=bound_for_p;p+=2)  {
              if(isprime(p))  {
                 ppi++;
                 mpz_set_ui(p2,p);
                 mpz_mul(p2,p2,p2);
                 // p2=p*p
                 mpz_set_ui(q,2*p);
                 mpz_add_ui(q,q,1);
                 // q=2*p+1
                 prod=1; 
                 // for very small primes we will use more than one primes for power-residue checking
                 // the idea behind this trick is (for non perfect numbers):
                 // if we use only one prime at power-residue checking
                 // then we expected to need O(log(log(bound_for_p)))=O(log(log(log(n)))) n-th root test
                 // using more than one prime we lower this number to O(1)
                 while(prod<=bound_for_p/p)  {
                       if(mpz_cmp(q,p2)>0)  {
                          // the number is too large for Lehmer test
                          iq=mpz_get_ui(q);
                          if(isprime(iq))  {
                             prime[count]=p;
                             mpz_set_ui(pprime[count],iq);
                             count++;
                             prod*=p;
                          }
                       }
                       else  {
                          // otherwise do a single Lehmer test
                          if(single_strong_lehmer_test(q,p))  {
                             prime[count]=p;
                             mpz_set(pprime[count],q);
                             count++;
                             prod*=p;
                          }
                       }
                       mpz_add_ui(q,q,2*p);
                 }
              }
          }

          tree=(mpz_t*)(malloc)(2*count*sizeof(mpz_t));
          if(count)  {    // if there is no prime, then do nothing
             for(i=0;i<2*count;i++)  mpz_init(tree[i]);

             // build up the product tree
             for(i=0;i<count;i++)    mpz_set(tree[count+i],pprime[i]);
             for(i=count-1;i>0;i--)  mpz_mul(tree[i],tree[2*i],tree[2*i+1]);

             // compute the remainders in the tree
             mpz_mod(tree[1],n,tree[1]);
             for(i=2;i<2*count;i++)
                 mpz_mod(tree[i],tree[i>>1],tree[i]);
          }

          pos=0;
          answer=0;
          foundanswer=0;
          for(i=0;(foundanswer==0)&&(i<ppi);i++)  {
              p=prime[pos];
              rootres=1;  // n is a power residue modulo p or divisible by p
              prod=1;
              while(prod<=bound_for_p/p)  {
                    prod*=p;
                    mpz_sub_ui(q,pprime[pos],1);
                    mpz_divexact_ui(q,q,p);
                    mpz_powm(q,tree[pos+count],q,pprime[pos]);
                    pos++;
                    // if n isn't divisible by q and is not a power residue then n is not a perfect power
                    if((mpz_cmp_ui(q,0)!=0)&&(mpz_cmp_ui(q,1)!=0))  rootres=0;
              }
              if(rootres)  {
                 if(mpz_root(root,u,p)!=0)  {
                    foundanswer=1;
                    answer=1;
                 }
              }
          }
          
          if((foundanswer==0)&&(sgn==1))  {
             if(mpz_perfect_square_p(u)!=0)  foundanswer=1,answer=1;
          }
          goto bigdel;

    smalldel:
        mpz_clear(u);
        mpz_clear(p2);
        mpz_clear(q);
        mpz_clear(root);
        return answer;

    bigdel:
          for(i=0;i<primepi;i++)  mpz_clear(pprime[i]);
          for(i=0;i<2*count;i++)  mpz_clear(tree[i]);
          free(tree);
          free(pprime);
          free(prime);
          goto smalldel;

    return 0;  // never reached
}
