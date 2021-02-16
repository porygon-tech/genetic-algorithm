# Genetic algorithms applied to 2020 SARS-CoV2 pandemic
## Theory
This kind of optimisation algorithm tries to emulate evolution by generating many different candidate solutions and then trying to improve them separately through applying a selective-pressure-like force on loci (parameters). The effect of this selective process is that beneficial alleles survive along time, while deletereous ones are renewed with new candidates.

New alleles are introduced in a population through the following processes:
- Mutation
- Immigration (from a different population)

Existing alleles are removed from population through:
- Selection
- Genetic drift
- (Emigration is negligible)

New combinations of alleles can arise by crossing-over their corresponding chromosomes.
Remember that natural selection do not act on alleles themselves, but on their phenotypic effect.

## Implementation

Search space must be discretised (split) into 2^n parts (where n is a positive integer). Thus, a \"phenotypic effect\" (a value belonging this discretised search space) can be encoded as a \"gene\" by rewritting it as a binary number of n digits corresponding to this value.

Each \"organism\" has a bunch of binary sequences associated that act as its \"chromosomes\".





More exploitatory variants (only useful when used on well-behaved functions!):
Elitist - takes only the fittest individuals.
Steady state - when selected parents generate their two children (with crossover and mutation), the two best ones out of this four-member family are selected. (¿Removes individuals at random (biased by fitness)?).

## Structures
Three different structures are created. The first one to store the individual, the second one to store the ODE 11 parameters and the last one to store the data.
```c
typedef struct {
	unsigned long IC[IC_GENES_NUMBER];          
	unsigned long Pars[PARAMETERS_GENES_NUMBER]; 
	double DeltaPars[PARAMETERS_GENES_NUMBER];  
	double fitness;
} individual;

typedef struct {
	double beta, phi, epsI, epsY, sigma, gamma1, gamma2, kappa, p, alpha, delta;
	unsigned PopSize;
} ODE_Parameters;

typedef struct {
        double PopSize;
	unsigned n_days;
	double Data_Time_Series[N_DAYS][N_VARS];	
} DataForFitting;
```

## Functions
### Random seed initialisation
`seedinit` initialises the random seed for random number generation using the system clock. Requires `#include <time.h>`
```c
void seedinit(void){
	time_t clock;
	srand((unsigned) time(&clock));
}
```
### Random uniform double generation
`uniform` generates random uniform decimal numbers between 0 and 1.
```c
double uniform(void){
	// Returns a number in [0,1)
	return (double)rand()/RAND_MAX;
}
```
### Binary sequence string conversion
here, n is a positive integer lower than 2^len, and len is the generated sequence length. It can be used as follows: 
```c
char *seq;
seq = itobin(Par2LSDiscPar(0.023), 10);
```
```c
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
```
this version only prints the result, it does not store it.
```c
void printbin(int n, int len){
	int k;
	for (int i = len-1; i >= 0; i--){
			k = n >> i;
			if (k & 1) printf("1");
			else       printf("0");
	}
	printf("\n");
}
```






