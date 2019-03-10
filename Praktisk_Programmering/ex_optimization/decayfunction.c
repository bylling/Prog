#include<math.h>
#include<gsl/gsl_vector.h>

struct expdata {int n; double *tim, *sig, *err;};
double decayfunction(double tim, double A, double T, double B){
	return A*exp(-tim/T)+B;
}


double master(const gsl_vector *x, void *params){
	double A = gsl_vector_get(x,0);
	double T = gsl_vector_get(x,1);
	double B = gsl_vector_get(x,2);
	struct expdata data = *(struct expdata *)params;
	double chi2=0;
	for(int i=0;i<data.n;i++){
		double tim=data.tim[i];
		double sig=data.sig[i];
		double err=data.err[i];
		chi2+=pow(decayfunction(tim,A,T,B)-sig,2)/err/err;
		}
	return chi2;
}
