#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include "gmp.h"


int main()  {
	
	int B,N,power2,h,i,I,j,p,magassag,sq,*isprime,primepi,*prime,darab,darab2;
	int *count,**tomb,POWER2,MAGASSAG,*hely,pozicio,kitevo;
	time_t seconds;
	FILE* out;
	mpz_t F,*A,*szorzatfa,*SZORZATFA,maximum;
	gmp_randstate_t state;  // random generáláshoz kell
	gmp_randinit_default(state);
	gmp_randseed_ui(state,time(NULL));
	
	printf("Polinom ideju faktorizalas egyszerre, sok 'kis' szamra\n");
	printf("D. J. Bernstein modszere alapjan.\n");
	printf("Hany bites random szamokat szeretnel (<63 bit) ? ");
	scanf("%d",&B);
	assert(0<B&&B<63);
	printf("Hanyat generaljak? ");
	scanf("%d",&N);

  seconds=time(NULL);
	mpz_init(F);
	
	A=(mpz_t*)(malloc)(N*sizeof(mpz_t));
	for(i=0;i<N;i++)  mpz_init(A[i]);
	// N darab B bites szám generálása
	for(i=0;i<N;i++)  {
	    mpz_urandomb(A[i],state,B);
		mpz_add_ui(A[i],A[i],2);
	}
    
	mpz_init(maximum);
	mpz_set_ui(maximum,0);
	for(i=0;i<N;i++)  {
		if(mpz_cmp(A[i],maximum)>0)  mpz_set(maximum,A[i]);
	}
	gmp_printf("A random szamok maximuma=%Zd\n",maximum);
	
	mpz_sqrt(maximum,maximum);
	sq=mpz_get_ui(maximum);
	isprime=(int*)(malloc)((sq+1)*sizeof(int));
	for(i=0;i<=sq;i++)  isprime[i]=1;
	isprime[0]=0;
	isprime[1]=0;
	// Erathoszteneszi szita a prímekre
	for(i=0;i*i<=sq;i++)  {
		if(isprime[i])  {
			for(j=i*i;j<=sq;j+=i)  isprime[j]=0;
		}
	}
	primepi=0;
	for(i=0;i<=sq;i++)  primepi+=isprime[i];
	prime=(int*)(malloc)(primepi*sizeof(int));
	primepi=0;
	for(i=0;i<=sq;i++)
	   if(isprime[i])  prime[primepi]=i,primepi++;	
        printf("Primek szama a maximum negyzetgyokeig=%d\n",primepi);
        free(isprime);
        hely=(int*)(malloc)(primepi*sizeof(int));

	power2=1;
	magassag=0;
	while(power2<N)  power2<<=1,magassag++;
	power2<<=1;
	// faktorizálandó számok szorzatfájának a felépítése
	szorzatfa=(mpz_t*)(malloc)(power2*sizeof(mpz_t));
	for(i=0;i<power2;i++)  mpz_init(szorzatfa[i]);
	for(i=0;i<power2;i++)  mpz_set_ui(szorzatfa[i],1);
	
	for(i=0;i<N;i++)  mpz_set(szorzatfa[i+power2/2],A[i]);
	
	for(h=magassag;h>=1;h--)
		for(i=(1<<(h-1));i<(1<<h);i++)  mpz_mul(szorzatfa[i],szorzatfa[2*i],szorzatfa[2*i+1]);
	
    
	tomb=(int**)(malloc)((N+power2/2)*sizeof(int*));
	count=(int*)(malloc)((N+power2/2)*sizeof(int));
	
	POWER2=1;
	while(POWER2<primepi)  POWER2<<=1;
	POWER2<<=1;
	SZORZATFA=(mpz_t*)(malloc)(POWER2*sizeof(mpz_t));
	for(i=0;i<POWER2;i++)  mpz_init(SZORZATFA[i]);
	
	for(i=1;i<N+power2/2;i++)  {
		if(i==1)  {  // fa gyokere
	       darab=primepi;
		   for(j=0;j<primepi;j++)  hely[j]=j;
		}
		else  {
			darab=count[i>>1];
			for(j=0;j<darab;j++)  hely[j]=tomb[i>>1][j];
		}
	
	// prímek szorzatfájának az elkészítése
	POWER2=1;
	MAGASSAG=0;
	while(POWER2<darab)  POWER2<<=1,MAGASSAG++;
	POWER2<<=1;
	
	for(I=0;I<darab;I++)  mpz_set_ui(SZORZATFA[I+POWER2/2],prime[hely[I]]);
	for(I=darab;I<POWER2/2;I++) mpz_set_ui(SZORZATFA[I+POWER2/2],1);
	
	for(h=MAGASSAG;h>=1;h--)
		for(I=(1<<(h-1));I<(1<<h);I++)  mpz_mul(SZORZATFA[I],SZORZATFA[2*I],SZORZATFA[2*I+1]);
	
	// maradék fa felépítése
	mpz_mod(SZORZATFA[1],szorzatfa[i],SZORZATFA[1]);
	for(j=2;j<darab+POWER2/2;j++)  mpz_mod(SZORZATFA[j],SZORZATFA[j>>1],SZORZATFA[j]);
	
	darab2=0;
	for(j=0;j<darab;j++)
	    if(mpz_cmp_ui(SZORZATFA[j+POWER2/2],0)==0)  darab2++;
	count[i]=darab2;
	tomb[i]=(int*)(malloc)(darab2*sizeof(int));
	darab2=0;
	for(j=0;j<darab;j++)  // ha a maradek nulla, akkor a prim osztja a szamot
	    if(mpz_cmp_ui(SZORZATFA[j+POWER2/2],0)==0)  tomb[i][darab2]=hely[j],darab2++;
    
	}
	
	out=fopen("ki.txt","a+");
	for(i=0;i<N;i++)  {
		gmp_fprintf(out,"%Zd=",A[i]);
		pozicio=i+power2/2;
		darab=count[pozicio];
		mpz_set(F,A[i]);
		for(j=0;j<darab;j++)  { // talált prímekkel beosztom
			p=prime[tomb[pozicio][j]];
			kitevo=0;
			while(mpz_divisible_ui_p(F,p)!=0)  {
				  kitevo++;
				  mpz_divexact_ui(F,F,p);
			}
			if(j)  fprintf(out,"*");
			if(kitevo==1)  fprintf(out,"%d",p);
			else  fprintf(out,"%d%c%d",p,'^',kitevo);
		}
		if(mpz_cmp_ui(F,1)>0)  {  // egy nagy prímosztó
			if(darab!=0)  fprintf(out,"*");
			gmp_fprintf(out,"%Zd",F);
		}
		fprintf(out,"\n");
	}
	fclose(out);	
	printf("Kesz %ld masodperc alatt.\n",(int) time(NULL)-seconds);
	printf("Faktorizalas eredmenye a ki.txt fileban talalhato.\n");
	
	return 0;
}
