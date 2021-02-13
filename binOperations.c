#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//gcc binOperations.c -o bop -lm -O3 -fsanitize=address -static-libasan -g -Wall && time ./bop

void printbin(int n, int len){
	int k;
	for (int i = len-1; i >= 0; i--){
			k = n >> i;
			if (k & 1) printf("1");
			else       printf("0");
	}
	printf("\n");
}

void seedinit(void){
	time_t clock;
	srand((unsigned) time(&clock));
}

double uniform(void){
	// Returns a number in [0,1)
	return (double)rand()/RAND_MAX;
}

void crossover(unsigned int *a, unsigned int *b, unsigned int ndigits) {
	unsigned int ax, bx;
	unsigned int cut = uniform()*(ndigits-1) + 1;
	unsigned int mask = pow(2,cut)-1;
	ax = (~mask & *a) | ( mask & *b);
	bx = ( mask & *a) | (~mask & *b);
	*a = ax;
	*b = bx;
}

void mutate(unsigned int *a, double mu, unsigned int ndigits){
	//alters a single digit in sequence
	double p = uniform();
	if( p < mu){
		//printf("%.4lf < %.4lf\n", p, mu);
		unsigned int cut = uniform()*ndigits;
		unsigned int point = pow(2,cut);
		*a = (point ^ *a);
	}
}

int main(int argc, char const *argv[])
{
	seedinit();

	unsigned int ndigits = 10, maxbits, d, e;
	double mu = 0.8f;	

	maxbits = pow(2,ndigits);
	 
	d = uniform()*maxbits;
	e = uniform()*maxbits;

	//printf("%d,%d\n", d,e);	
	printbin(d, ndigits);
	printbin(e, ndigits);

	printf("\ncrossover:\n");
	crossover(&d,&e,ndigits);
	printbin(d, ndigits);
	printbin(e, ndigits);
	//printf("%d,%d\n", d,e); 	
 	printf("\nmutation(?):\n");
	printbin(d, ndigits);
	//mutate(&d, mu, ndigits);
	mutate(&d, 1, ndigits);
	printbin(d, ndigits);	

	//printf("%d\n", d);
	//r = d & e;
	//printbin(r, ndigits);
	return 0;
}
