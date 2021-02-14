
#define n_days 101 // Number of days in time series
#define n_vars 5   // Number of variables in time series

double Params2norm(double *parameters, doble *x0, void *TheData void){
      register unsigned ndays;
      DataForFitting *TheData = (DataForFitting *) TheData void;
      double t=0.0, err, h=1.e-3, norm=0.0;
      ODE_Parameters ODE pars = { parameters[0], parameters[1], parameters[2], parameters[3],
                                  parameters[4], parameters[5], parameters[6], parameters[7],
                                  parameters[8], parameters[9], parameters[10], TheData->PopSize };
      double xt[CoreModelDim]={x0[0], x0[1], x0[2], x0[3], x0[4], x0[5], x0[6], x0[7] };
      
      for (ndays=1; ndays <= TheData->n_days; ndays++) {
          int status;
          while (t+h < ndays) {
                status=RKF78Sys(&t, xt, CoreModelDIM, &h, &err, HMIN, HMAX, RKTOL, &ODE pars, CoreModel);
                if (status) return MAXDOUBLE;
          }
          h=ndays-t;
          status=RKF78Sys(&t, xt, CoreModelDIM, &h, &err, HMIN, HMAX, RKTOL, &ODE pars, CoreModel);
          if (status) return MAXDOUBLE;
          
          double d0=xt[4]-TheData -> Data_time_Series[ndays][0];
          double d1=xt[5]-TheData -> Data_time_Series[ndays][1];
          double d2=xt[6]-TheData -> Data_time_Series[ndays][2];
          double d3=xt[7]-TheData -> Data_time_Series[ndays][3];
          double d4=TheData -> PopSize - (xt[0]+xt[1]+xt[2]+xt[3]+xt[4]+xt[5]+xt[6]+xt[7]) -TheData -> Data_time_Series[ndays][4];
          norm += ndays * (d0*d0+d1*d1+d2*d2+d3*d3+d4*d4);
      }
      return norm;      
}

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











