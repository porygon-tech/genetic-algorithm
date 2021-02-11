
double Parameters2norm(double 
	DataForFitting *TheData = (DataForFitting *) TheData_void;
	double t = 0.0, err, h = 1.e-3, norm = 0.0; 
	ODE_Parameters ODE_pars 

/*
double xt1CoreModelDIM1 = 
'parameters, double x0, void TheData_void 

  { Paraneters[01, parameters111, parameters parameters[61, parameters171. parameters ( x0[01, x0111, x0121, x0131, x0141, x0151 
for(ndays=1: ndays <= TheData>N Days;	
 ndays+4) ( int status;
  while(t+h < ndays) status = RKF78Sys(&t, xt, CoreModelDIM, &h. &err. MIN, MAX. RKTOL, &ODE pars. CoreModet);  
   if(status) return MAXDOUBLE: ) h  ndays t;    status  RKF78Syst&t, xt, CoreModelDIM, &h, &err, HMIN, HMAX, RKTOL, &ODE pars, CoreModel);    
    if(status) return HAXDOUBLE;     
1) register unsigned ndays; 
[2). parameters(31, parameters(41, parameters151. [81. parameters191. parameters1101, TheData>PopSize );
 , x0[61, x0171 };  
double de = xt141 - TheData->Data Time Series1ndays1101: double dl  xt151 TheData->Data Time Series1ndays1111;
double d2  xt(6)  TheData->Data Time Sertes1ndays1 (21 ;
double d3 = xt[71 - TheData->Data Time SeriesInday$1131;
double d4 = TheData->PopSize (xt[01+xt111+xt121+xt131+xt141+xt norm +. ndays  (0de + d2(12 + d3d3 + d4.14): Mora noted 
*define crom2iC(c) (((double) c)/10001 *define crom2HSPar(c) (((double) c)/1099511627776UL) *define crom2Par(c) (((double) c)/1048576U) *define crom2LSPar(c) (((double) c)/1024U) *define Par2HSDiscPar(c) (c  1099511627776UL) *define Per2DiscPar(c) (c  1048576U) *define Par2LSDiscPar(c) (c  1024U) 
*/





int Steepest_Descent_backtracking(double *, double *, double *, double (*) (double *, double *, void *), double *, void *);

void CoreModelOuadraticErrorFitness(individual *ind, void TheData_void) {DataForFitting TheData = (DataForFitting *) TheData_void;
	double ic[CoreModelDIM] = { TheData -> PopSize, crom2IC(ind->IC[0]), crom2IC(ind->IC[1]), crom2IC(ind->IC[2]), 
	                            TheData -> Data_Time_Series[0][0]
	                            TheData -> Data_Time_Series[0][3] }

	ind -> DeltaPars[0]  =  crom2HSPar (ind->Pars[0]  + ind->DeltaPars[0] );
	ind -> DeltaPars[9]  =  crom2HSPar (ind->pars[9]  + ind->DeltaPars[9] );
	ind -> DeltaPars[10] =  crom2HSPar (ind->pars[10] + ind->DeltaPars[10]);
	ind -> DeltaPars[1]  =  crom2Par   (ind->pars[1]  + ind->DeltaPars[1] );
	ind -> DeltaPars[4]  =  crom2Par   (ind->pars[4]  + ind->DeltaPars[4] );
	ind -> DeltaPars[5]  =  crom2Par   (ind->pars[5]  + ind->DeltaPars[5] );
	ind -> DeltaPars[6]  =  crom2Par   (ind->pars[6]  + ind->DeltaPars[6] );
	ind -> DeltaPars[8]  =  crom2Par   (ind->pars[8]  + ind->DeltaPars[8] );
	ind -> DeltaPars[2]  =  crom2LSPar (ind->pars[2]  + ind->DeltaPars[2] );
	ind -> DeltaPars[3]  =  crom2LSPar (ind->pars[3]  + ind->DeltaPars[3] );
	ind -> DeltaPars[7]  =  crom2LSPar (ind->pars[7]  + ind->DeltaPars[7] );

	int niters = Steepest_Descent_backtracking(ind->DeltaPars, &(ind->fitness), &fitness_not_improved, Parameters2norm, ic, TheData_void). 

/*
...........................
*/

}











