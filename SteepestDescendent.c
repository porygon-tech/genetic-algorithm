// Steepest Descendent Method Backtracking 
#include <stdio.h>
#include <math.h>
#include "GeneticAlgorithm.h"

#define TOL 1.0e-6
#define ZEROT 1.e-14
#define ITER_MAX 1000

void copy_vect (double *dest, double *source, unsigned int n) {
   register unsigned int i;
   for (i=0; i<n; i++, dest++, source++) *dest = *source;
}

void LinComb_vect (double *dest, double *v1, double *v2, double alpha, unsigned int n) {
   register unsigned int i;
   for (i=0; i<n; i++) dest[i] = v1[i] - alpha *v2[i];
} 

double dot_product_vect (double *v, unsigned int n) {
   register unsigned int i;
   double dot = 0.0;
   for (i=0; i<n; i++, v++) {
      dot += (*v) * (*v);
   }
   return dot;
}

double max_norm_vect (double *v, unsigned int n) {
   register unsigned int i;
   double x, norm = 0.0;
   for (i=0; i<n; i++){
      if (norm < (x = fabs(v[i]))) norm = x;
   }
   return norm;
}

double forward_gradient ( double *, unsigned int, double (*)(double *, double *, void *), double *, double *, void *, double, double *);


int Steepest_Descent_backtracking ( double *x, double *fx, double *fx_base, double (*f)(double *, double *, void *), double *ic, void *TheData) {
   
   register unsigned iter = 0U;
   double gradient[PARAMETERS_GENES_NUMBER], grad_abs_error_inf_norm, sigma = 0.1, rho = 0.5;
   double gradient_max_norm = forward_gradient (gradient, PARAMETERS_GENES_NUMBER, f, x, ic, TheData, 1.0e-6, &grad_abs_error_inf_norm);
   *fx_base = *fx = f(x, ic, TheData);
   
   while ( gradient_max_norm >= TOL ) {
       double x_new_try[PARAMETERS_GENES_NUMBER], fx_new_try, alpha = gradient_max_norm;
       
       for ( i=0; i<PARAMETERS_GENES_NUMBER; i++){
          gradient[i] /= gradient_max_norm;
          if (fabs(gradient[i]) < ZEROT || x[i] < ZEROT || x[i] > ONET) {
             gradient[i] = 0.0; continue;
          }
          double alpha_new = (gradient[i] > 0.0) ? x[i] / gradient[i] : (x[i] - 1.0) / gradient[i];
          if (alpha_new < alpha) alpha = alpha_new;
       }
       
       double sigma_gradient_L2norm = sigma * dot_product_vect(gradient, PARAMETERS_GENES_NUMBER);
       
       LinComb_vect(x_new_try, x, alpha, gradient, PARAMETERS_GENES_NUMBER); 
       fx_new_try = f(x_new_try, ic, TheData);
       
       while ( fx_new_try > *fx + alpha * sigma_gradient_L2norm){
          alpha *= rho; 
          LinComb_vect(x_new_try, x, alpha, gradient, PARAMETERS_GENES_NUMBER);
          fx_new_try = f(x_new_try, ic, TheData);
       } ;
       if (*fx < fx_new_try + TOL || alpha * max_norm_vect(gradient, PARAMETERS_GENES_NUMBER) < TOL) return iter;
       
       double relative_fitness_error = (*fx - fx_new_try) / (*fx);
       copy_vect(x, x_new_try, PARAMETERS_GENES_NUMBER); *fx = fx_new_try; iter ++;
       
       if ( iter >= ITER_MAX || relative_fitness_error < TOL ) return iter;
       
       gradient_max_norm = forward_gradient ( gradient, PARAMETERS_GENES_NUMBER, f, x, ic, TheData, TOL, &grad_abs_error_ind_norm);
   }
   return iter;
}

// Numerical computation of the gradient
#define GSL EPSILON 2.220446049250313e-16
#define GSL MAX(a, b) ((a)>(b) ? (a):(b))

//This is the important one, because we see the increments

double partial_derivative_forward (double (*f)(double *, double *, void *), unsigned int i,
                                   double *x, double *ic, void *TheData,
                                   double h, double *abserr round, double *abserr trunc) {

/* Compute the derivative using the 4point rule (x+h/4, x+h/2, x+3h/4, x+h)
   Compute the error using the difference between the 4point and the 2point rule (x+h/2, x+h) */
   double xi = x[i];
   x[i] = xi+h/4.0;
   double f1 = f(x, ic, TheData);
   x[i] = xi+h/2.0;
   double f2 = f(x, ic, TheData);
   x[i] = xi+h*(3.0/4.0);
   double f3 = f(x, ic, TheData);
   x[i] = xi+h;
   double f4 = f(x, ic, TheData);
   
   double r2 = 2.0*(f4-f2);
   double r4 = (22.0/3.0)*(f4-f3) - (62.0/3.0)*(f3-f2) + (52.0/3.0)*(f2-f1);
// Estimate the rounding error for r4
   double e4 = 2 * 20.67 * (fabs(f4)+fabs(f3)+fabs(f2)+fabs(f1)) * GSL EPSILON;

// The next term is due to finite precision in x+h=O(eps*x)
   double dy = GSL MAX(fabs(r2/h), fabs(r4/h)) * fabs(xi/h) * GSL EPSILON;

/* The truncation error in the r4 approximation itself is O(h³)
   However, we estimate the error from r4-r2, which is O(h)
   By scaling h we minimise this estimated error, not the actual truncation error in r4 */
   *abserr trunc = fabs((r4-r2)/h); //Estimated truncation error O(h)
   *abserr round = fabs(e4/h)+dy;
   return (r4/h);
}

// Computes the forward partial derivatives
double forward_partial_derivative (double (*f)(double *, double *, void *), unsigned int i,
                                   double *x, double *ic, void *TheData,
                                   double h, double *abserr) {
   double round, trunc;
   double r 0 = partial_derivative_forward (f, i, x, ic, TheData, h, &round, &trunc);
   *abserr = round + trunc;
   if (!(round < trunc && (round > 0 && trunc > 0))) retunr r 0;
   
/* Compute an optimised stepsize to minimize the total error, 
   using the scaling of the estimated truncation error O(h) and 
   rounding error O(1/h) */
   double r opt = partial_derivative_forward (f, i, x, ic, TheData, h*sqrt(round/trunc), &round, &trunc); round += trunc;
   
/* Check that the new error is smaller, and that the new derivative is consistent
   with the error bounds of the original estimate */
   if (!(round < *abserr && fabs(r opt - r 0) < 4.0 * (*abserr))) return r 0;
   
   *abserr = round;
   return r opt;                                 
}

//Computes the 11 derivatives for the gradient
double forward_gradient (double *grad, unsigned int n, double (*f)(double *, double *, void *),
                         double *x, double *ic, void *TheData,
                         double h, double *grad_abs_error_inf_norm){
   register unsigned int i;
   double max_norm = -1.0; *grad_abs_error_inf_norm = -1.0;
   for(i=0; i<n; i++) {
      double abserr;
      *(grad+i) ) forward_partial_derivative (f, i, x, xi, TheData, h, &abserr);
      if(abserr > *grad_abs_error_inf_norm) *grad_abs_error_inf_norm = abserr;
      if(max_norm < (abserr = fabs(*(grad+i)))) max_norm = abserr;
   }
   return max_norm;
}




