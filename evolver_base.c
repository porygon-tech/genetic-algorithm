#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "RKF78.h"

/*DEBUGGING:  pra
gcc evolver_base.c RKF78-2.2.c/RKF78.c -o evo -lm -O3 -fsanitize=address -static-libasan -g -Wall && time ./evo
*/

#define IC_GENES_NUMBER 3
#define PARAMETERS_GENES_NUMBER 11

#define CORE_MODEL_DIM 8
#define HMAX 1.0
#define HMIN 1.e-3
#define RKTOL 1.e-5

#define MAXDOUBLE 1.7976931348623157e+308

#define N_DAYS 101 // Number of days in time series
#define N_VARS 5   // Number of variables in time series

#define POPSIZE 2000 // Lo que recomienda Alseda


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

typedef struct {
	double beta, phi, epsI, epsY, sigma, gamma1, gamma2, kappa, p, alpha, delta;
	unsigned PopSize;
} ODE_Parameters;

typedef struct {
        double PopSize;
	unsigned n_days;
	double Data_Time_Series[N_DAYS][N_VARS];	
} DataForFitting;

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

void mutate2(unsigned int *a, double mu, unsigned int ndigits){
	//lifelike mutation effect
	double p;
	unsigned long int points = 0, add = 0;
	for (int i = 0; i < ndigits; ++i){
		p = uniform();
		if( p < mu){
			add = (int) pow(2, i);
			points = points ^ add;			
		}
		//printbin(points, ndigits); printf(" %d\n", points);
	}
	*a = (points ^ *a);
}

//-----------------------------------
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

//-----------------------------------

int main(int argc, char const *argv[])
{
	seedinit();	
	//int PopSize = 2000;

     /* Real Data (DataForFittint, TheData_void) */
	float RD[N_DAYS][N_VARS] = {
	 	{1.000, 1.000, 0.000, 0.000, 0.000},	
		{1.841, 1.253, 0.056, 0.348, 0.000},
		{2.285, 1.607, 0.123, 0.733, 0.002},
		{2.571, 2.041, 0.203, 1.169, 0.004},
		{2.812, 2.541, 0.300, 1.670, 0.008},
		{3.059, 3.096, 0.414, 2.249, 0.013},
		{3.337, 3.702, 0.546, 2.915, 0.020},
		{3.655, 4.356, 0.698, 3.681, 0.029},
		{4.018, 5.060, 0.870, 4.555, 0.040},
		{4.428, 5.814, 1.062, 5.548, 0.054},
		{4.889, 6.626, 1.276, 6.673, 0.071},
		{5.402, 7.499, 1.513, 7.941, 0.091},
		{5.972, 8.442, 1.773, 9.364, 0.115},
		{6.603, 9.462, 2.058, 10.959, 0.143},
		{7.299, 10.570, 2.369, 12.739, 0.175},
		{8.067, 11.776, 2.710, 14.723, 0.212},
		{8.914, 13.091, 3.083, 16.931, 0.253},
		{9.847, 14.529, 3.490, 19.382, 0.301},
		{10.875, 16.103, 3.934, 22.102, 0.355},
		{12.008, 17.830, 4.420, 25.116, 0.415},
		{13.257, 19.725, 4.952, 28.454, 0.483},
		{14.633, 21.808, 5.534, 32.147, 0.559},
		{16.150, 24.099, 6.172, 36.231, 0.644},
		{17.823, 26.620, 6.871, 40.746, 0.738},
		{19.667, 29.395, 7.638, 45.735, 0.844},
		{21.701, 32.453, 8.480, 51.245, 0.960},
		{23.943, 35.822, 9.406, 57.331, 1.090},
		{26.416, 39.535, 10.422, 64.051, 1.233},
		{29.143, 43.628, 11.540, 71.469, 1.392},
		{32.151, 48.141, 12.770, 79.656, 1.569},
		{35.469, 53.116, 14.124, 88.693, 1.763},
		{39.128, 58.603, 15.615, 98.666, 1.979},
		{43.164, 64.654, 17.256, 109.670, 2.217},
		{47.616, 71.328, 19.064, 121.812, 2.480},
		{52.527, 78.688, 21.056, 135.210, 2.770},
		{57.944, 86.805, 23.252, 149.991, 3.091},
		{63.918, 95.758, 25.671, 166.298, 3.446},
		{70.508, 105.632, 28.338, 184.289, 3.837},
		{77.777, 116.523, 31.278, 204.136, 4.269},
		{85.795, 128.535, 34.520, 226.031, 4.745},
		{94.638, 141.782, 38.094, 250.183, 5.271},
		{104.391, 156.393, 42.035, 276.826, 5.851},
		{115.149, 172.506, 46.380, 306.216, 6.492},
		{127.013, 190.277, 51.172, 338.635, 7.198},
		{140.099, 209.874, 56.456, 374.394, 7.978},
		{154.529, 231.485, 62.282, 413.837, 8.838},
		{170.444, 255.316, 68.708, 457.342, 9.787},
		{187.994, 281.594, 75.793, 505.327, 10.833},
		{207.347, 310.569, 83.606, 558.252, 11.988},
		{228.687, 342.516, 92.222, 616.623, 13.261},
		{252.217, 377.737, 101.722, 681.000, 14.666},
		{278.160, 416.566, 112.197, 751.997, 16.215},
		{306.763, 459.370, 123.746, 830.294, 17.924},
		{338.296, 506.551, 136.480, 916.637, 19.809},
		{373.057, 558.552, 150.519, 1011.850, 21.887},
		{411.372, 615.862, 165.995, 1116.839, 24.180},
		{453.603, 679.015, 183.056, 1232.602, 26.708},
		{500.144, 748.598, 201.861, 1360.238, 29.496},
		{551.431, 825.258, 222.588, 1500.956, 32.570},
		{607.940, 909.702, 245.431, 1656.087, 35.960},
		{670.196, 1002.705, 270.602, 1827.096, 39.697},
		{738.773, 1105.119, 298.338, 2015.593, 43.818},
		{814.301, 1217.874, 328.894, 2223.347, 48.361},
		{897.473, 1341.989, 362.553, 2452.305, 53.369},
		{989.042, 1478.579, 399.625, 2704.605, 58.889},
		{1089.838, 1628.858, 440.449, 2982.593, 64.974},
		{1200.765, 1794.155, 485.396, 3288.848, 71.680},
		{1322.811, 1975.914, 534.874, 3626.195, 79.069},
		{1457.053, 2175.710, 589.326, 3997.735, 87.212},
		{1604.666, 2395.253, 649.239, 4406.864, 96.183},
		{1766.929, 2636.396, 715.142, 4857.302, 106.065},
		{1945.231, 2901.150, 787.612, 5353.119, 116.949},
		{2141.079, 3191.686, 867.280, 5898.765, 128.936},
		{2356.106, 3510.345, 954.827, 6499.101, 142.133},
		{2592.078, 3859.646, 1050.997, 7159.429, 156.661},
		{2850.898, 4242.291, 1156.593, 7885.531, 172.651},
		{3134.614, 4661.169, 1272.484, 8683.701, 190.246},
		{3445.426, 5119.362, 1399.610, 9560.781, 209.600},
		{3785.685, 5620.137, 1538.981, 10524.199, 230.885},
		{4157.899, 6166.949, 1691.682, 11582.004, 254.286},
		{4564.733, 6763.432, 1858.876, 12742.907, 280.004},
		{5009.003, 7413.384, 2041.802, 14016.312, 308.258},
		{5493.676, 8120.751, 2241.781, 15412.350, 339.286},
		{6021.857, 8889.602, 2460.210, 16941.914, 373.345},
		{6596.777, 9724.093, 2698.565, 18616.677, 410.714},
		{7221.773, 10628.429, 2958.390, 20449.122, 451.691},
		{7900.265, 11606.813, 3241.302, 22452.549, 496.601},
		{8635.723, 12663.380, 3548.971, 24641.082, 545.789},
		{9431.630, 13802.127, 3883.121, 27029.662, 599.628},
		{10291.435, 15026.829, 4245.506, 29634.031, 658.513},
		{11218.496, 16340.939, 4637.900, 32470.693, 722.867},
		{12216.019, 17747.481, 5062.074, 35556.871, 793.138},
		{13286.985, 19248.930, 5519.771, 38910.428, 869.799},
		{14434.065, 20847.086, 6012.679, 42549.783, 953.349},
		{15659.535, 22542.936, 6542.399, 46493.794, 1044.309},
		{16965.181, 24336.518, 7110.407, 50761.622, 1143.224},
		{18352.193, 26226.782, 7718.015, 55372.566, 1250.660},
		{19821.068, 28211.456, 8366.327, 60345.878, 1367.197},
		{21371.503, 30286.926, 9056.199, 65700.547, 1493.434},
		{23002.296, 32448.131, 9788.183, 71455.071, 1629.977},
		{24711.260, 34688.480, 10562.489, 77627.195, 1777.437}
	};
	
	DataForFitting realData;
        realData.PopSize = 1000000;
	memcpy(realData.Data_Time_Series, RD, sizeof(float)*N_DAYS*N_VARS);
	//printf("%f\n", realData.Data_Time_Series[2][1]);
	
     /* Initial Population */
        individual *population;
	population = (individual *) malloc(sizeof(individual) * POPSIZE);

	for (int i = 0; i < POPSIZE; ++i){
	     indiv_init(&population[i]); 	     
	}
	//printf("%f\n", population -> fitness);
     
        individual p;
	indiv_init(&p);
	//printind(p);
     
        CoreModelVersusDataQuadraticError(&p, &realData);
        printf("%.50f\n", p.fitness);
     /* Compute Fitness 
        for (int i=0; i < POPSIZE; ++i){
	     CoreModelVersusDataQuadraticError(&population[i], &realData);
	}
	printf("%f\n", population -> fitness);
     */

/*	//MUTATION TEST
	unsigned int ndigits = 20, chr;
	double mu = 0.1f;	

	double d;
	d = uniform(); printf("%.6lf\n", d);
	chr = Par2DiscPar(d);
	printbin(chr,ndigits);
	mutate2(&chr, mu, ndigits);
	printbin(chr,ndigits);
	d = crom2Par(chr); printf("%.6lf\n", d);
*/


/*	//INDIVIDUAL INITIALISATION TEST
	individual p;
	indiv_init(&p);
	printind(p);
*/


/*	//PHENOTYPE-GENOTYPE ENCODING TEST
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

	free(population);
	return 0;
}



