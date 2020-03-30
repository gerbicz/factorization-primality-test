#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "gmp.h"

#define reps 5


int issquarefree(int n)  {
   int i;
   for(i=2;i*i<=n;i++)
       if(n%(i*i)==0)  return 0;
   return 1;
}

int isprime(int n)  {
   int i;
   if(n<2)  return 0;
   for(i=2;i*i<=n;i++)
       if(n%i==0)  return 0;
   return 1;
}

double F(mpz_t n,int mult,int B)  {// Knuth-Schroeppel szorzohoz
  
  mpz_t N,temp;
  mpz_init(N);
  mpz_init(temp);
  mpz_mul_ui(N,n,mult);
  
  int p,res;
  double s=-0.5*log(mult),v;
  for(p=2;p<B;p++)
      if(isprime(p)&&(p==2||mult%p==0||mpz_ui_kronecker(p,N)==1))  {
         if(p==2)  {
            res=mpz_mod_ui(temp,N,8);
            if(mult%2==0||res%4==3)  v=1.0/3.0;
            else if(res==5)          v=2.0/3.0;
            else                     v=4.0/3.0;
         }
         else  {
            if(mult%p==0)            v=(double) 1/(p+1);
            else                     v=(double) 2*p/(p*p-1);
        }
        s+=v*log(p);
     }

  mpz_clear(N);
  mpz_clear(temp);
  return s;
}
  
int szorzo(mpz_t n,int B)  {// legjobb szorzo megkeresese B-ig, 256-ig szitalva
  
  int mult,ret;
  double s,rec=-1000.0;
  
  for(mult=1;mult<256;mult++)
      if(issquarefree(mult))  {
         s=F(n,mult,B);
         if(s>rec)  {rec=s;ret=mult;}
      }
   return ret;
}


int CFRAC(mpz_t n)  {
	
	int primepi,k,l,ret,faktorbazis_meret,f32,omega,osztok_szama,*pr,*isp,ftlen,faktorizalt,len,mult,B1,EAS,MM[32];
	unsigned int *isprime,*kitevo,*pozicio,*om,**matrix,**offset,**pos2,**E,*off,*v,u,Bit[32];
	long long int i,j,pow2;
	EAS=3;
	
	mpz_t a0,a2,a,h0,h1,h2,h_st,k0,k1,k2,m,m2,d,d2,maradek,gmptemp,prod,prod2,R,G,GG,nn;
	mpz_t *oszto,*B,*eredeti,EASbound[EAS];
	
	mpz_init(a0);
	mpz_init(a2);
	mpz_init(a);
	mpz_init(h0);
	mpz_init(h1);
	mpz_init(h2);
	mpz_init(h_st);
	mpz_init(k0);
	mpz_init(k1);
	mpz_init(k2);
	mpz_init(m);
	mpz_init(m2);
	mpz_init(d);
	mpz_init(d2);
	mpz_init(maradek);
	mpz_init(gmptemp);
	mpz_init(prod);
	mpz_init(prod2);
	mpz_init(R);
	mpz_init(G);
	mpz_init(GG);
	mpz_init(nn);

	Bit[0]=1;
	for(i=1;i<32;i++)  Bit[i]=2*Bit[i-1];// Bit[i]=2^i
	
	for(i=0;i<EAS;i++) mpz_init(EASbound[i]);
  oszto=(mpz_t*)(malloc)(1024*sizeof(mpz_t));
	for(i=0;i<1024;i++)  mpz_init(oszto[i]);
  isp=(int*)(malloc)(1024*sizeof(int));

//L(n)=exp((log(n)*log(log(n)))^0.5)
//a=1/sqrt((6+2/(k+1)))
//c[i]=4*i*a*a/(k+1)^2
//teta[i]=i/(k+1)

//B1=L(n)^a
//Ha L^(teta[i]*a)-ig faktorizalatlan>n^(1/2*(1-c[1]-..-c[i])), akkor kilep

  mpz_set(nn,n);
  osztok_szama=0;

  int n2=mpz_sizeinbase(n,2)+8;
  double dlog2=(double) n2*log(2);
  double dc,de=0.5,da=(double) 1.0/sqrt(6+(double) 2/(EAS+1));
  
  de=sqrt(dlog2*log(dlog2))*da/log(2);
  B1=1<<((int) de);
  if(B1<8)  B1=8;
  for(i=2;i<=B1;i++)
      if(mpz_divisible_ui_p(n,(int)i)!=0)  {
         mpz_set_ui(gmptemp,(int)i);
         mpz_set_ui(oszto[osztok_szama],(int)i);
         isp[osztok_szama++]=1;
         j=mpz_remove(n,n,gmptemp);
         if(mpz_cmp_ui(n,1)==0)  goto del;
      }

  mpz_set(oszto[osztok_szama],n);
  isp[osztok_szama]=mpz_probab_prime_p(n,reps);
	osztok_szama++;

  mult=szorzo(n,256);
  mpz_mul_ui(n,n,mult);

  n2=mpz_sizeinbase(n,2);
  dlog2=(double) n2*log(2);
  de=0.5,da=(double) 1.0/sqrt(6+(double) 2/(EAS+1));
  
  de=sqrt(dlog2*log(dlog2))*da/log(2);
  B1=1<<((int) de);
  if(B1<8)  B1=8;
  
  de=0.5;
  for(i=0;i<EAS;i++)  {
      dc=(double) 4*(i+1)*da*da/((EAS+1)*(EAS+1));
      mpz_ui_pow_ui(gmptemp,B1,(int) i+1);
      mpz_root(gmptemp,gmptemp,EAS+1);
      MM[i]=mpz_get_ui(gmptemp);
      
      de-=0.5*dc;
      mpz_pow_ui(gmptemp,n,(int) ((double) 64.0*de));
      mpz_root(gmptemp,gmptemp,64);
      mpz_set(EASbound[i],gmptemp);
  }
  
  len=(B1+31)/32;
  isprime=(unsigned int*)(malloc)(len*sizeof(unsigned int));
  for(i=0;i<len;i++)  isprime[i]=0xffffffff;
  isprime[0]&=~Bit[0];
  isprime[0]&=~Bit[1];
  for(i=0;i*i<B1;i++)
      if((isprime[i>>5]&Bit[i&31])>0)  {
         for(j=i*i;j<B1;j+=i)  isprime[j>>5]&=~Bit[j&31];
      }
      
  primepi=1;
  for(i=0;i<len;i++)  {
      u=isprime[i];
      while(u)  primepi+=u&1,u>>=1;
  }
  pr=(int*)(malloc)(primepi*sizeof(int));

  pr[0]=-1;// -1 mindig benne van a faktorbazisban
	faktorbazis_meret=1;
	for(i=0;i<B1;i++)// faktorbazisban levo primek megkeresese B1-ig
	    if((isprime[i>>5]&Bit[i&31])>0&&(i==2||mult%i==0||mpz_kronecker_ui(n,(int)i)==1))
	        pr[faktorbazis_meret++]=i;
  f32=(faktorbazis_meret+31)>>5;

  B=(mpz_t*)(malloc)((faktorbazis_meret+1)*sizeof(mpz_t));
	for(i=0;i<=faktorbazis_meret;i++)  mpz_init(B[i]);

  eredeti=(mpz_t*)(malloc)((faktorbazis_meret+1)*sizeof(mpz_t));
	for(i=0;i<=faktorbazis_meret;i++)  mpz_init(eredeti[i]);

  matrix=(unsigned int**)(malloc)((faktorbazis_meret+1)*sizeof(unsigned int*));
  for(i=0;i<=faktorbazis_meret;i++)  matrix[i]=(unsigned int*)(malloc)(f32*sizeof(unsigned int));
  offset=(unsigned int**)(malloc)((faktorbazis_meret+1)*sizeof(unsigned int*));
  for(i=0;i<=faktorbazis_meret;i++)  offset[i]=(unsigned int*)(malloc)(f32*sizeof(unsigned int));
  pos2=(unsigned int**)(malloc)((faktorbazis_meret+1)*sizeof(unsigned int*));
  for(i=0;i<=faktorbazis_meret;i++)  pos2[i]=(unsigned int*)(malloc)((faktorbazis_meret+1)*sizeof(unsigned int));
  E=(unsigned int**)(malloc)((faktorbazis_meret+1)*sizeof(unsigned int*));
  for(i=0;i<=faktorbazis_meret;i++)  E[i]=(unsigned int*)(malloc)((faktorbazis_meret+1)*sizeof(unsigned int));

  kitevo=(unsigned int*)(malloc)((faktorbazis_meret+1)*sizeof(unsigned int));
  off=(unsigned int*)(malloc)(f32*sizeof(unsigned int));
  v=(unsigned int*)(malloc)(f32*sizeof(unsigned int));
  om=(unsigned int*)(malloc)((faktorbazis_meret+1)*sizeof(unsigned int));
  pozicio=(unsigned int*)(malloc)((faktorbazis_meret+1)*sizeof(unsigned int));

		mpz_set_ui(m,0);
		mpz_set_ui(d,1);
		mpz_sqrt(a,n);
		mpz_set(a0,a);
		
		mpz_set_ui(h0,1);
		mpz_set(h1,a0);
		mpz_set_ui(k0,0);
		mpz_set_ui(k1,1);
		
		mpz_set_ui(h_st,1);
	  pow2=1;// ciklizálás deketálásához
		
		for(i=0;i<faktorbazis_meret;i++)
		    for(j=0;j<f32;j++)  matrix[i][j]=0;
		
		if(!isp[osztok_szama-1])  {
    faktorizalt=0;
    
		for(i=0;!faktorizalt;i++)  {
		   	// lanctort kovetkezo tagjanak eloallitasa
		   	mpz_mul(m2,a,d);
		   	mpz_mod(m2,m2,n);
		   	mpz_sub(m2,m2,m);
		   	mpz_set(m,m2);
		   	
        mpz_mul(d2,m,m);
        mpz_mod(d2,d2,n);
        mpz_sub(d2,n,d2);
        mpz_divexact(d,d2,d); // pontos osztas
        
        mpz_add(a,a0,m2);
        mpz_fdiv_q(a,a,d);
        
        mpz_mul(h2,a,h1);
        mpz_add(h2,h2,h0);
        mpz_mod(h2,h2,n);
        
        mpz_mul(k2,a,k1);
        mpz_add(k2,k2,k0);
        mpz_mod(k2,k2,n);

			if(mpz_cmp(h2,h_st)==0)  {ret=-1;goto del;}// ciklizalt a modszer
			if(i==pow2)  {pow2<<=1;mpz_set(h_st,h2);}
			
			mpz_set(h0,h1);
			mpz_set(h1,h2);
			
			mpz_set(k0,k1);
			mpz_set(k1,k2);
			
			mpz_mul(maradek,h2,h2);
			mpz_mod(maradek,maradek,n);
			
			if(i&1)  {
			   mpz_sub(maradek,n,maradek);
			   omega=1;
			   kitevo[0]=1;
			   pozicio[0]=0;
			}
			else  omega=0;
			
			int pos=0;
			for(j=1;j<faktorbazis_meret;j++)  {// trial division
			    if(mpz_divisible_ui_p(maradek,pr[j])!=0)  {
			       mpz_set_ui(gmptemp,pr[j]);
			       kitevo[omega]=mpz_remove(maradek,maradek,gmptemp);
			       pozicio[omega]=j;
			       omega++;
			       if(mpz_cmp_ui(maradek,1)==0)  break;
			    }
			    if(pos<EAS&&pr[j]>MM[pos])  {// early abort strategy (EAS)
			       if(mpz_cmp(maradek,EASbound[pos])>0)  break;
			       pos++;
			    }
			}
			
			if(mpz_cmp_ui(maradek,1)==0)  {// B1-beli primek szorzata volt
			   
			   mpz_set(prod,h2);
			   
			   for(k=0;k<f32;k++)  off[k]=0;
			   for(k=0;k<f32;k++)  v[k]=0;
			   for(j=0;j<omega;j++)
			       if(kitevo[j]&1)  v[pozicio[j]>>5]+=Bit[pozicio[j]&31];
			   
			   ftlen=0;
			   for(j=pozicio[0];j<faktorbazis_meret;j++)
			       if((v[j>>5]&Bit[j&31])>0)  {
			          if((matrix[j][j>>5]&Bit[j&31])==0)  {ftlen=1;break;}
			          mpz_mul(prod,prod,B[j]);
			          mpz_mod(prod,prod,n);
			          for(k=0;k<f32;k++)  off[k]^=offset[j][k];
			          for(k=j>>5;k<f32;k++) v[k]^=matrix[j][k];
			      }
			      
			 if(ftlen)  {
			     off[j>>5]|=Bit[j&31];
			     for(k=j>>5;k<f32;k++)  matrix[j][k]=v[k];
			     for(k=0;k<f32;k++)     offset[j][k]=off[k];
			     mpz_set(B[j],prod);
			     mpz_set(eredeti[j],h2);
			     
			     om[j]=omega;
			     for(l=0;l<omega;l++)  {
			         pos2[j][l]=pozicio[l];
			         E[j][l]=kitevo[l];
			     }
			 }
			 else  {
			     om[faktorbazis_meret]=omega;
			     for(l=0;l<omega;l++)  {
			         pos2[faktorbazis_meret][l]=pozicio[l];
			         E[faktorbazis_meret][l]=kitevo[l];
			     }
			     
			     for(l=0;l<faktorbazis_meret;l++)  kitevo[l]=0;
			     
			     mpz_set_ui(prod2,1);
			     for(k=0;k<=faktorbazis_meret;k++)
			         if(k==faktorbazis_meret||(off[k>>5]&Bit[k&31])>0)  {
			            if(k!=faktorbazis_meret)  {
			               mpz_set(R,eredeti[k]);
			               mpz_mul(prod2,prod2,R);
			               mpz_mod(prod2,prod2,n);
			            }
			            else  {
			               mpz_set(R,h2);
			               mpz_mul(prod2,prod2,R);
			               mpz_mod(prod2,prod2,n);
			            }
			            for(l=0;l<(int) om[k];l++)
			               kitevo[pos2[k][l]]+=E[k][l];
			     }
			   
			     mpz_set_ui(prod,1);
			     for(l=0;l<faktorbazis_meret;l++)  {
			         mpz_set_si(gmptemp,pr[l]);
			         mpz_powm_ui(gmptemp,gmptemp,kitevo[l]>>1,n);
			         mpz_mul(prod,prod,gmptemp);
			         mpz_mod(prod,prod,n);
			     }
			     
			     mpz_sub(G,prod,prod2);
			     mpz_gcd(G,G,n);
			     
			     if(mpz_cmp_ui(G,1)>0&&mpz_cmp(G,n)<0)  {// osztot talaltunk
			        
			        for(j=0;j<osztok_szama;j++)
			            if(!isp[j])  {
			               mpz_gcd(GG,G,oszto[j]);
			               if(mpz_cmp_ui(GG,1)>0&&mpz_cmp(GG,oszto[j])<0)  {
			                  mpz_divexact(oszto[osztok_szama],oszto[j],GG);
			                  isp[osztok_szama]=mpz_probab_prime_p(oszto[osztok_szama],reps);
			                  mpz_set(oszto[j],GG);
			                  isp[j]=mpz_probab_prime_p(oszto[j],reps);
			                  osztok_szama++;
			               }
			            }
			        
			        faktorizalt=1;
			        for(j=0;j<osztok_szama;j++)
			            if(!isp[j])  {faktorizalt=0;break;}
			        
			        if(faktorizalt)  {ret=1;goto del;}// n-et faktirzaltuk
			     }
	}}}}
	
 del:
	
	for(i=0;i<osztok_szama;i++)
	    for(j=0;i+j+1<osztok_szama;j++)
	        if(mpz_cmp(oszto[j],oszto[j+1])>0)  {
	           mpz_set(gmptemp,oszto[j]);
	           mpz_set(oszto[j],oszto[j+1]);
	           mpz_set(oszto[j+1],gmptemp);
	        }
	
	FILE* out;
	out=fopen("cfrac.txt","a+");
	if(!faktorizalt)  fprintf(out,"nem teljes fakt. ");
	
	mpz_set(n,nn);
	gmp_fprintf(out,"%Zd=",n);
	int elso=1;
	for(i=0;i<osztok_szama;i++)  {
	    int e=mpz_remove(n,n,oszto[i]);
	    if(e>0)  {
	       if(!elso)  fprintf(out,"*");
	       elso=0;
	       gmp_fprintf(out,"%Zd",oszto[i]);
	       if(e>1)  fprintf(out,"%c%d",'^',e);
	    }
	}
	fprintf(out,"\n");
	fclose(out);

	 free(pr);
	 free(isp);
	 free(isprime);
	 free(kitevo);
	 free(pozicio);
	 free(om);
	 free(off);
	 free(v);
	 for(i=0;i<=faktorbazis_meret;i++)  free(matrix[i]);free(matrix);
	 for(i=0;i<=faktorbazis_meret;i++)  free(offset[i]);free(offset);
	 for(i=0;i<=faktorbazis_meret;i++)  free(pos2[i]);free(pos2);
	 for(i=0;i<=faktorbazis_meret;i++)  free(E[i]);free(E);

	 for(i=0;i<=faktorbazis_meret;i++)  mpz_clear(eredeti[i]);free(eredeti);
	 for(i=0;i<=faktorbazis_meret;i++)  mpz_clear(B[i]);free(B);

   for(i=0;i<1024;i++)  mpz_clear(oszto[i]);
	
	 mpz_clear(a0);
	 mpz_clear(a2);
	 mpz_clear(a);
	 mpz_clear(h0);
	 mpz_clear(h1);
	 mpz_clear(h2);
	 mpz_clear(h_st);
	 mpz_clear(k0);
	 mpz_clear(k1);
	 mpz_clear(k2);
	 mpz_clear(m);
	 mpz_clear(m2);
	 mpz_clear(d);
	 mpz_clear(d2);
	 mpz_clear(maradek);
	 mpz_clear(gmptemp);
	 mpz_clear(prod);
	 mpz_clear(prod2);
	 mpz_clear(R);
	 mpz_clear(G);
	 mpz_clear(GG);
	 mpz_clear(nn);
	
   return ret;
}


int main()  {

    mpz_t n;
    mpz_init(n);
    FILE* in;
    char filename[100];
    
    printf("Faktorizacio lanctortek segitsegevel (CFRAC).\nInput file neve? ");
    scanf("%s",filename);
    
    time_t sec=time(NULL);
    in=fopen(filename,"r");

    while(gmp_fscanf(in,"%Zd",&n)!=EOF) CFRAC(n);
    fclose(in);
    printf("Befejezve %ld masodperc alatt.\n",time(NULL)-sec);
    printf("cfrac.txt-ben lasd a faktorizaciokat.\n");
    
    mpz_clear(n);
    return 0;
}
