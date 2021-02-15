#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*DEBUGGING:  
gcc evolver_base.c RKF78-2.2.c/RKF78.c -o evo -lm -O3 -fsanitize=address -static-libasan -g -Wall && time ./evo
*/

#define IC_GENES_NUMBER 3
#define PARAMETERS_GENES_NUMBER 11

#define CoreModelDIM 8
#define HMAX 1.0
#define HMIN 1.e-3
#define RKTOL 1.e-5

#define n_days 101 // Number of days in time series
#define n_vars 5   // Number of variables in time series

#define PopSize 1000000

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
	//printf("\n");
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

void encode(individual *ind){
	ind->Pars[0]  = Par2HSDiscPar(ind->DeltaPars[0]);
	ind->Pars[9]  = Par2HSDiscPar(ind->DeltaPars[9]);
	ind->Pars[10] = Par2HSDiscPar(ind->DeltaPars[10]);
	ind->Pars[1]  = Par2DiscPar  (ind->DeltaPars[1]);
	ind->Pars[4]  = Par2DiscPar  (ind->DeltaPars[4]);
	ind->Pars[5]  = Par2DiscPar  (ind->DeltaPars[5]);
	ind->Pars[6]  = Par2DiscPar  (ind->DeltaPars[6]);
	ind->Pars[8]  = Par2DiscPar  (ind->DeltaPars[8]);
	ind->Pars[2]  = Par2LSDiscPar(ind->DeltaPars[2]);
	ind->Pars[3]  = Par2LSDiscPar(ind->DeltaPars[3]);
	ind->Pars[7]  = Par2LSDiscPar(ind->DeltaPars[7]);
}

void indiv_init(individual *ind){
	for (int i = 0; i < PARAMETERS_GENES_NUMBER; ++i){
		ind->DeltaPars[i] = uniform();
	}
	encode(ind);
}

void printind(individual p){
	printf("%.12f  ",         p.DeltaPars[0]); printbin(p.Pars[0]  , 40); printf("\n");
	printf("%.12f  ",         p.DeltaPars[9]); printbin(p.Pars[9]  , 40); printf("\n");
	printf("%.12f  ",         p.DeltaPars[10]);printbin(p.Pars[10] , 40); printf("\n");
	printf("%.6f        ",    p.DeltaPars[1]); printbin(p.Pars[1]  , 20); printf("\n");
	printf("%.6f        ",    p.DeltaPars[4]); printbin(p.Pars[4]  , 20); printf("\n");
	printf("%.6f        ",    p.DeltaPars[5]); printbin(p.Pars[5]  , 20); printf("\n");
	printf("%.6f        ",    p.DeltaPars[6]); printbin(p.Pars[6]  , 20); printf("\n");
	printf("%.6f        ",    p.DeltaPars[8]); printbin(p.Pars[8]  , 20); printf("\n");
	printf("%.3f           ", p.DeltaPars[2]); printbin(p.Pars[2]  , 10); printf("\n");
	printf("%.3f           ", p.DeltaPars[3]); printbin(p.Pars[3]  , 10); printf("\n");
	printf("%.3f           ", p.DeltaPars[7]); printbin(p.Pars[7]  , 10); printf("\n");
}

//-----------------------------------



int main(int argc, char const *argv[])
{
	/*
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
	*/
	seedinit();

	individual p;
	indiv_init(&p);

	printind(p);


	return 0;
}



