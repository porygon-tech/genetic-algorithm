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
### Fitness function computation
The function to compute the Vector Field, the function to generate predictions from individuals, the function to obtain the norm, and the function to merge all this functions and return the fitness of each individual inside a population; are the ones provided in the assignment with a few changes:

```c
void CoreModel(double t, double *x, unsigned CoreModelDIM, double *der, void *Params){
	ODE_Parameters *par = (ODE_Parameters *) Params; // To simplify the usage of Params (void pointer)
	
	double sigmae   = par -> sigma * x[1],
	       gamma1i1 = par -> gamma1 * x[2],
	       kappaA   = par -> kappa*x[3],
	       alphai2  = par -> alpha*x[5];

	der[0] = par -> phi*x[2] + x[3] + (1-par -> epsI)*(x[4]+x[5]) + (1-par -> epsY)*x[6];
	der[0] = - par -> beta * (x[0] * der[0])/par -> PopSize;
	der[1] = - der[0] - sigmae;
	der[2] = sigmae - gamma1i1;
	der[3] = (1-par -> p)*gamma1i1 - kappaA - par -> gamma2*x[3] ;
	der[4] = kappaA - par -> gamma2*x[4];
	der[5] = par -> p*gamma1i1 - par -> gamma2*x[5] - alphai2;
	der[6] = alphai2 - (par -> gamma2+par -> delta)*x[6];
	der[7] = par -> gamma2*(x[3] + x[4] + x[5] + x[6]);
}

int GeneratePredictionFromIndividual(double *xt, void *ODE_pars, DataForFitting *Pred) {
	register unsigned ndays;
	double t = 0.0, err, h = 1.e-3;
	for (ndays=1; ndays <= Pred -> n_days; ndays++) { 
		int status;
		while(t+h < ndays) {
			status = RKF78Sys(&t, xt, CORE_MODEL_DIM, &h, &err, HMIN, HMAX, RKTOL, ODE_pars, CoreModel);
			if(status) return status;
		}
		h = ndays - t;
		status = RKF78Sys(&t, xt, CORE_MODEL_DIM, &h, &err, HMIN, HMAX, RKTOL, ODE_pars, CoreModel);
		if(status) return status;
		Pred -> Data_Time_Series[ndays][0] = xt[4];
		Pred -> Data_Time_Series[ndays][1] = xt[5];
		Pred -> Data_Time_Series[ndays][2] = xt[6];
		Pred -> Data_Time_Series[ndays][3] = xt[7];
		Pred -> Data_Time_Series[ndays][4] = Pred -> PopSize - (xt[0]+xt[1]+xt[2]+xt[3]+xt[4]+xt[5]+xt[6]+xt[7]);
	}
	return 0;
}

double Parameters2norm(double *parameters, double *x0, void *TheData_void){
      register unsigned ndays;
      DataForFitting *TheData = (DataForFitting *) TheData_void;
      double t=0.0, err, h=1.e-3, norm=0.0;
      ODE_Parameters ODE_pars = { parameters[0], parameters[1], parameters[2], parameters[3],
                                  parameters[4], parameters[5], parameters[6], parameters[7],
                                  parameters[8], parameters[9], parameters[10], TheData->PopSize };
      double xt[CORE_MODEL_DIM]={x0[0], x0[1], x0[2], x0[3], x0[4], x0[5], x0[6], x0[7] };
      
      for (ndays=1; ndays<=TheData->n_days; ndays++) {
          int status;
          while (t+h < ndays) {
                status=RKF78Sys(&t, xt, CORE_MODEL_DIM, &h, &err, HMIN, HMAX, RKTOL, &ODE_pars, CoreModel);
                if (status) return MAXDOUBLE;
          }
          h=ndays-t;
          status=RKF78Sys(&t, xt, CORE_MODEL_DIM, &h, &err, HMIN, HMAX, RKTOL, &ODE_pars, CoreModel);
          if (status) return MAXDOUBLE;
          
          double d0=xt[4]-TheData -> Data_Time_Series[ndays][0];
          double d1=xt[5]-TheData -> Data_Time_Series[ndays][1];
          double d2=xt[6]-TheData -> Data_Time_Series[ndays][2];
          double d3=xt[7]-TheData -> Data_Time_Series[ndays][3];
          double d4=TheData -> PopSize - (xt[0]+xt[1]+xt[2]+xt[3]+xt[4]+xt[5]+xt[6]+xt[7]) - TheData -> Data_Time_Series[ndays][4];
          norm += ndays * (d0*d0+d1*d1+d2*d2+d3*d3+d4*d4);
      }
      return norm;      
}

void CoreModelVersusDataQuadraticError(individual *ind, void *TheData) {
	DataForFitting *TDfF = (DataForFitting *) TheData;
	DataForFitting ThePrediction = { 
		TDfF -> PopSize, 
		TDfF -> n_days, {
		{ TDfF -> Data_Time_Series[0][0], 
		  TDfF -> Data_Time_Series[0][1],
		  TDfF -> Data_Time_Series[0][2],
		  TDfF -> Data_Time_Series[0][3],
		  TDfF -> Data_Time_Series[0][4]
		}}
	};

	double xt[CORE_MODEL_DIM] = { 
		TDfF -> PopSize, 
		crom2IC(ind -> IC[0]), 
		crom2IC(ind -> IC[1]),
		crom2IC(ind -> IC[2]), 
		TDfF -> Data_Time_Series[0][0],
		TDfF -> Data_Time_Series[0][1], 
		TDfF -> Data_Time_Series[0][2],
		TDfF -> Data_Time_Series[0][3] 
	};

	xt[0] -= (xt[1]+xt[2]+xt[3]+xt[4]+xt[5]+xt[6]);

	ODE_Parameters ODE_pars = { 
		ind -> DeltaPars[0], 
		ind -> DeltaPars[1],
		ind -> DeltaPars[2], 
		ind -> DeltaPars[3],
		ind -> DeltaPars[4], 
		ind -> DeltaPars[5],
		ind -> DeltaPars[6], 
		ind -> DeltaPars[7],
		ind -> DeltaPars[8], 
		ind -> DeltaPars[9],
		ind -> DeltaPars[10], 
		TDfF -> PopSize 
	};

	if(GeneratePredictionFromIndividual(xt, &ODE_pars, &ThePrediction)) {
		ind -> fitness = MAXDOUBLE;
		return;
	};
	ind -> fitness = Parameters2norm(ind->DeltaPars, xt, &ThePrediction);
}
```




