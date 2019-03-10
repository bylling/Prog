#include<stdio.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>
#include<math.h>

double integrand (double t, void * params) {
  double f = log(t)/sqrt(t);
	return f;
}

double integrandnorm (double t, void * params) {
  double alpha = *(double*)params;
  double f = exp(-alpha * t*t);
  return f;
}
double integrandham (double t, void * params) {
  double alpha = *(double*)params;
  double f = (-alpha * alpha * t*t/2 + alpha/2 + t * t/2) * exp(-alpha * t*t);
  return f;
}


int main(){
printf("Exercise 1\n");
    int limit=100;
  	gsl_integration_workspace * w;
  	w = gsl_integration_workspace_alloc (limit);
    gsl_function F;
  	 F.function = integrand;
    double result,error,acc=1e-8,eps=1e-8;
    int flag = gsl_integration_qags
    		(&F, 0, 1, acc, eps, limit, w, &result, &error);
    gsl_integration_workspace_free(w);
    if(flag!=GSL_SUCCESS) return printf("Integral_1 = %g\n",NAN);
    else printf("Integral_1 = %g\n",result);
printf("Exercise 1 is done\n");

printf("Exercise 2\n");
int limit2=100;
gsl_integration_workspace * w2;
w2 = gsl_integration_workspace_alloc (limit2);
gsl_function F2;
F2.function = integrandnorm;

int limit3=100;
gsl_integration_workspace * w3;
w3 = gsl_integration_workspace_alloc (limit3);
gsl_function F3;
F3.function = integrandham;


for(double alpha=0.1;alpha<20;alpha+=0.1){

  F2.params = (void*)&alpha;
	double result,error,acc2=1e-8,eps2=1e-8;
	int flag2 = gsl_integration_qagi
		(&F2, acc2, eps2, limit2, w2, &result, &error);
  	if(flag2!=GSL_SUCCESS) printf("Integration error in <>.");

  F3.params = (void*)&alpha;
  double result2,error2,acc3=1e-8,eps3=1e-8;
  int flag3 = gsl_integration_qagi
    (&F3, acc3, eps3, limit3, w3, &result2, &error2);
    if(flag3!=GSL_SUCCESS) printf("Integration error in <H>");

  fprintf(stderr,"%g %g %g %g\n",alpha,result,result2,result2/result);
}
gsl_integration_workspace_free(w2);
gsl_integration_workspace_free(w3);
printf("Exercise 2 is done and plotted\n");
printf("We know the ground state of the 1d harmonic oscillator is a gaussian with alpha = 1\n");
printf("This groundstate has excactly an energy of 1/2 hbar, so we have just found the ground state\n");
return 0;
}
