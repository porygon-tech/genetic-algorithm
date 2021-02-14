#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*DEBUGGING:  
gcc evolver_base.c RKF78-2.2.c/RKF78.c -o evo -lm -O3 -fsanitize=address -static-libasan -g -Wall && time ./evo
*/

//antidiscretising functions
#define crom2IC(c) (((double) c)/1000)
#define crom2HSPar(c) (((double) c)/1099511627776UL)
#define crom2Par(c) (((double) c)/1048576U)
#define crom2LSPar(c) (((double) c)/1024U)

//discretising functions
#define Par2HSDiscPar(c) (c * 1099511627776UL)
#define Par2DiscPar(c) (c * 1048576U)
#define Par2LSDiscPar(c) (c * 1024U)

typedef struct {
	unsigned long IC[IC_GENES_NUMBER];
	unsigned long Pars[PARAMETERS_GENES_NUMBER];
	double DeltaPars[PARAMETERS_GENES_NUMBER];
	double fitness;
} individual;

//-----------------------------------
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

char *itobin(int n, int len){
	int k, t;
	char *p;
	t = 0;
	p = (char*)malloc(len+1);
	if (p == NULL) exit(EXIT_FAILURE);

	for (int i = len-1; i >= 0; i--){
			k = n >> i;
			if (k & 1) *(p+t) = 1 + '0';
			else       *(p+t) = 0 + '0';
			t++;
	}
	printf("\n");
	*(p+t) = '\0';

	return  p;
}

//-----------------------------------



int main(int argc, char const *argv[])
{
	double d;


	d = 0.023;
	printf("%.3lf\n", d);
	printbin(Par2LSDiscPar(d),10);

	d = 0.000001;
	printf("%.6lf\n", d);
	printbin(Par2DiscPar(d),20);

	d = 0.000000000003;
	printf("%.12lf\n", d);
	printbin(Par2HSDiscPar(d),40);

	char *seq;
	seq = itobin(Par2LSDiscPar(0.023),10);
	printf("%s\n", seq);



	return 0;
}



