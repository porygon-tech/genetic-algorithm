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
	return (double)rand()/RAND_MAX;
}; // Returns a number in [0,1)

/*
	
void OnePointCrossover( unsigned int p1, 
	                    unsigned int p2,
                        unsigned int *f1, 
                        unsigned int *f2){
	unsigned char d = uniform()*( 8* sizeof( unsigned int)- 1) + 1;

	unsigned int mask = 0xFFFFFFFFU << d;
	*f1 = ( p1 & mask) | ( p2 & ~mask);
	*f2 = ( p2 & mask) | ( p1 & ~mask);
}
*/
/*
void OnePointCrossover( unsigned int p1, 
	                    unsigned int p2,
                        unsigned int *f1, 
                        unsigned int *f2){
	unsigned char len = 8*sizeof(unsigned int);
	unsigned char d = uniform()*(len- 1) + 1, di = len - d;
	*f1 = ((p1>>d)<<d) | ((p2<<di)>>di);
	*f2 = ((p2>>d)<<d) | ((p1<<di)>>di);
}
*/


int main(int argc, char const *argv[])
{

	seedinit();

	unsigned int ndigits = 10, maxbits, d, e, r;
	//unsigned char d, e, r;
	maxbits = pow(2,ndigits);
	//maxbits = ndigits*sizeof(unsigned int);
	
	d = uniform()*maxbits;
	e = uniform()*maxbits;
	printf("length = %d (%d bits)\n%d, %d\n", ndigits, maxbits, d, e);
	printbin(d, ndigits);
	printbin(e, ndigits); 
	r = d & e;
	printbin(r, ndigits);
	return 0;
}
