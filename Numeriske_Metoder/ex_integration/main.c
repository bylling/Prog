#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_integration.h>
#include"integrationfunctions.h"


double my_function(double x){return x*x;}
int main (void)
{

  // Exercise A
  fprintf(stdout, "Exercise A has started\n\n" );
  fprintf(stdout, "An recursive adaptive integrator have been implemented. It can with a given accuracy-goal find the integral through a recursive calculation, with a error-calulation based on a lower order method. For each layer of recursion, points are saved and reused when sepperating in smaller steps. \n" );
  fprintf(stdout, "The integrator is tested on the following integrals, with returned values, errors and number of recursive iterations: \n\n\n" );

  fprintf(stdout, "Test 1: I = ∫ from 0 to 1 of √(x) dx = 2/3 = 0.66..\n" );

  // We initialze the integration:
  int calls1 = 0;
  double err1 = 0;
  double a=0, b=1, acc=1e-4, eps=1e-4;
  // We setup the funciton
  double f(double x){calls1++; return sqrt(x);};
  // We calculate the integral
  fprintf(stdout, "The calculation is started. \n" );
  double Q1 = adapt(f, a, b, acc, eps, &err1);
  fprintf(stdout, "The calculation is done. \n" );

  // We print the values
	fprintf(stdout,"The found values are: integral I=%lg, error=%lg, calls=%i\n\n",Q1,err1, calls1);


  // Second integral
  fprintf(stdout, "Test 2: I = ∫ from 0 to 1 of 1/√(x) dx = 2 \n" );

  // We initialze the integration:
  int calls2 = 0;
  double err2 = 0;
  a=0, b=1, acc=1e-4, eps=1e-4;
  // We setup the funciton
  double f2(double x){calls2++; return 1./sqrt(x);};
  // We calculate the integral
  fprintf(stdout, "The calculation is started. \n" );
  double Q2 = adapt(f2, a, b, acc, eps, &err2);
  fprintf(stdout, "The calculation is done. \n" );
  // We print the values
	fprintf(stdout,"The found values are: integral I=%lg, error=%lg, calls=%i\n\n",Q2,err2, calls2);


  // Third integral
  fprintf(stdout, "Test 3: I = ∫ from 0 to 1 of ln(x)/√(x) dx = -4 \n" );

  // We initialze the integration:
  int calls3 = 0;
  double err3 = 0;
  a=0, b=1, acc=1e-4, eps=1e-4;
  // We setup the funciton
  double f3(double x){calls3++; return log(x)/sqrt(x);};
  // We calculate the integral
  fprintf(stdout, "The calculation is started. \n" );
  double Q3 = adapt(f3, a, b, acc, eps, &err3);
  fprintf(stdout, "The calculation is done. \n" );
  // We print the values
	fprintf(stdout,"The found values are: integral I=%lg, error=%lg, calls=%i\n\n",Q3,err3, calls3);


  fprintf(stdout, "The tests have succeded, all with the same wanted accuracy, and as is clearly seen in the number of calls, some functions are more difficult to estimate than others, especially when dealing with diverging terms as the logarithm and 1/x. \n" );
  fprintf(stdout, "For serious calculations, we want to calculate the following integral with as many significant digits as possible, until we at some point reaches a divergence in calculation-cost of the recursive calculator, probably due to the restriction from the machine-epsilon. \n\n" );

  // Fourth integral
  fprintf(stdout, "Calculation of I = ∫ from 0 to 1 of 4 *√(1-(1-x)2) dx = π  \n" );

  // We initialze the integration:
  int calls4 = 0;
  double err4 = 0;
  a=0, b=1;
  acc=1e-19, eps=1e-19;
  // We setup the funciton
  double f4(double x){calls4++; return 4. * sqrt(1.-(1.-x)*(1.-x));};
  // We calculate the integral
  fprintf(stdout, "The calculation is started. \n" );
  double Q4 = adapt(f4, a, b, acc, eps, &err4);
  fprintf(stdout, "The calculation is done. \n" );
  // We print the values
  fprintf(stdout, "The found values are: \n Known solution: π = 3.141592653589793238462 \n Found solution I =%25.23g \n With error=%lg and calls=%i\n Beyond this point of accuracy, the calulation time increases dramatically. \n",Q4,err4, calls4);
  fprintf(stdout, "As is seen at this point we reach an accuracy smaller than the machine-epsilon, which results in the fact that the recursive relation will not be able to converge in a simple manner any more. \n" );
  fprintf(stdout, "With this successfull demonstration of the tecursive adaptive integrator, the first exercise is hereby finished.\n\n\n" );

  // Exercise B
  fprintf(stdout, "Exercise B has started\n\n" );
  fprintf(stdout, "An recursive adaptive integrator have been implemented using the Clenshaw-Curtis variable transform. This transformation can for some integrals improve accuracy, while for others reduce accuracy, dependent on whether the transformation suits the function, in order to direct a easier function for the recursive solver. We test this integrator on the same integrals as before, since a few of these diverges at the borders, and thereby can get an advantage when being transformed. Furthermore we calculate the same integrals by the GSL-procedure for comparison, here we use the Gauss-Kronrod 21 integration algorithm QAGS to handle some the the diverging regimes. The calculation and comparison proceeds: \n\n" );


  fprintf(stdout, "Test 1: I = ∫ from 0 to 1 of √(x) dx = 2/3 = 0.66..\n" );

  // We initialze the integration:
  int calls11 = 0;
  double err11 = 0;
  a=0, b=1, acc=1e-4, eps=1e-4;
  // We setup the funciton again, with the new calls
  double f11(double x){calls11++; return sqrt(x);};
  // We calculate the integral
  fprintf(stdout, "The calculation is started. \n" );
  double Q11 = adapt_clenshaw_curtis(f11, a, b, acc, eps, &err11);
  fprintf(stdout, "The calculation is done. \n" );

  // We calculate the GSL-integration regimes as done in practical programming
  int calls111 = 0;
  int limit=100;
  gsl_integration_workspace * w;
  w = gsl_integration_workspace_alloc (limit);
  // We setup the funciton again, with the new calls
  double f111(double x,void * params){calls111++; return sqrt(x);};
  gsl_function  F;
  F.function = f111;
  double Q111;
  double err111;
  gsl_integration_qags(&F, 0, 1, acc, eps, limit, w, &Q111, &err111);
  gsl_integration_workspace_free(w);

  // We print the values
	fprintf(stdout,"The previous values are:    integral I=%lg, error=%lg, calls=%i\n",Q1,err1, calls1);
  fprintf(stdout,"The transformed values are: integral I=%lg, error=%lg, calls=%i\n",Q11,err11, calls11);
  fprintf(stdout,"The GSL values are:         integral I=%lg, error=%lg, calls=%i\n\n",Q111,err111, calls111);


  // Second integral
  fprintf(stdout, "Test 2: I = ∫ from 0 to 1 of 1/√(x) dx = 2 \n" );

  // We initialze the integration:
  int calls22 = 0;
  double err22 = 0;
  a=0, b=1, acc=1e-4, eps=1e-4;
  // We setup the funciton again, with the new calls
  double f22(double x){calls22++; return 1./sqrt(x);};
  // We calculate the integral
  fprintf(stdout, "The calculation is started. \n" );
  double Q22 = adapt_clenshaw_curtis(f22, a, b, acc, eps, &err22);
  fprintf(stdout, "The calculation is done. \n" );

  // We calculate the GSL-integration regimes as done in practical programming
  int calls222 = 0;
  limit=100;
  gsl_integration_workspace * w2;
  w2 = gsl_integration_workspace_alloc (limit);
  // We setup the funciton again, with the new calls
  double f222(double x,void * params){calls222++; return 1./sqrt(x);};
  gsl_function F2;
  F2.function = f222;
  double Q222;
  double err222;
  gsl_integration_qags(&F2, 0, 1, acc, eps, limit, w2, &Q222, &err222);
  gsl_integration_workspace_free(w2);

  // We print the values
  fprintf(stdout,"The previous values are:    integral I=%lg, error=%lg, calls=%i\n",Q2,err2, calls2);
  fprintf(stdout,"The transformed values are: integral I=%lg, error=%lg, calls=%i\n",Q22,err22, calls22);
  fprintf(stdout,"The GSL values are:         integral I=%lg, error=%lg, calls=%i\n\n",Q222,err222, calls222);

  // Third integral
  fprintf(stdout, "Test 3: I = ∫ from 0 to 1 of ln(x)/√(x) dx = -4 \n" );

  // We initialze the integration:
  int calls33 = 0;
  double err33 = 0;
  a=0, b=1, acc=1e-4, eps=1e-4;
  // We setup the funciton again, with the new calls
  double f33(double x){calls33++; return log(x)/sqrt(x);};
  // We calculate the integral
  fprintf(stdout, "The calculation is started. \n" );
  double Q33 = adapt_clenshaw_curtis(f33, a, b, acc, eps, &err33);
  fprintf(stdout, "The calculation is done. \n" );

  // We calculate the GSL-integration regimes as done in practical programming
  int calls333 = 0;
  limit=100;
  gsl_integration_workspace * w3;
  w3 = gsl_integration_workspace_alloc (limit);
  // We setup the funciton again, with the new calls
  double f333(double x,void * params){calls333++; return log(x)/sqrt(x);};
  gsl_function F3;
  F3.function = f333;
  double Q333;
  double err333;
  gsl_integration_qags(&F3, 0, 1, acc, eps, limit, w3, &Q333, &err333);
  gsl_integration_workspace_free(w3);


  // We print the values
  fprintf(stdout,"The previous values are:    integral I=%lg, error=%lg, calls=%i\n",Q3,err3, calls3);
  fprintf(stdout,"The transformed values are: integral I=%lg, error=%lg, calls=%i\n",Q33,err33, calls33);
  fprintf(stdout,"The GSL values are:         integral I=%lg, error=%lg, calls=%i\n\n",Q333,err333, calls333);

  // Fourth integral
  fprintf(stdout, "Calculation of I = ∫ from 0 to 1 of 4 *√(1-(1-x)2) dx = π  \n" );

  // We initialze the integration:
  int calls44 = 0;
  double err44 = 0;
  a=0, b=1;
  acc=1e-19, eps=1e-19;
  // We setup the funciton
  double f44(double x){calls44++; return 4. * sqrt(1.-(1.-x)*(1.-x));};
  // We calculate the integral
  fprintf(stdout, "The calculation is started. \n" );
  double Q44 = adapt_clenshaw_curtis(f44, a, b, acc, eps, &err44);
  fprintf(stdout, "The calculation is done. \n" );

  // We calculate the GSL-integration regimes as done in practical programming
  int calls444 = 0;
  acc = 1e-13;
  eps = 1e-13;
  limit=10000;
  gsl_integration_workspace * w4;
  w4 = gsl_integration_workspace_alloc (limit);
  // We setup the funciton again, with the new calls
  double f444(double x,void * params){calls444++; return 4. * sqrt(1.-(1.-x)*(1.-x));};
  gsl_function F4;
  F4.function = f444;
  double Q444;
  double err444;
  gsl_integration_qags(&F4, 0, 1, acc, eps, limit, w4, &Q444, &err444);
  gsl_integration_workspace_free(w4);


  // We print the values
  fprintf(stdout, "The found values are: \n Known solution:      π =   3.141592653589793238462 \n Previous solution    I = %25.23g \n Transformed solution I = %25.23g \n GSL solution         I =  %25.23g \n With previous   error=%lg  and calls=%i\n New transformed error=%lg  and calls=%i\n and GSL         error=%lg   and calls=%i\n Beyond this point of accuracy, the calulation time increases dramatically. \n",Q4,Q44,Q444,err4, calls4,err44, calls44,err444,calls444);
  fprintf(stdout, "From the calculations, we see that for integral test 2 and 3, where the logarithm and 1/x will diverge in the interval, the transformation significantly improves the integration procedure. However for the others we recieve a similar error, but with more function-calls, which might be problematic if the function is significantly difficult to calculate. \n" );
  fprintf(stdout, "Comparing to the GSL-routines, when the function diverges, theese have significantly better accuracy than our routines but also require a bit more functioncalls. However this can be changed by changing the GSL-algoritm used in the routine, in which there are a dusin. We see at test 1, for a simple function, our minmalistic first approach still provides the best result, however in regions with singularities the others are superior. An interesting aspect is the fact that for the last calculation, the GSL-method will be performed way faster and in way less functioncalls. Thereby it seems to be superior at fast calculations with high accuracy. However the restrictions on the GSL-solver is a bit less, in order not to have errors when it compiles. \n" );
  fprintf(stdout, "With this successfull demonstration of the recursive adaptive integrator using a Clenshaw-Curtis transformation, and a comparison with the GSL-routines, the second exercise is hereby finished.\n\n\n" );


  // Exercise C
  fprintf(stdout, "Exercise C has started\n\n" );
  fprintf(stdout, "An recursive adaptive integrator have been implemented using a variable transform to be able to calculate from a value 'a' to infinity. Furthermore we calculate the same integrals by the GSL-procedure for comparison, here we use the gsl_integration_qagiy integration algorithm, which actually works by the same transformation. The calculation and comparison proceeds: \n\n" );


    fprintf(stdout, "As a test function we use: I = ∫ from 0 to inf of x/(exp(x)-1) dx = pi²/6.\n" );

    // We initialze the integration:
    int callsinf = 0;
    double errinf = 0;
    a=0, acc=1e-12, eps=1e-12;
    // We setup the funciton again, with the new calls
    double finf(double x){callsinf++; return x/(exp(x)-1);};
    // We calculate the integral
    fprintf(stdout, "The calculation is started. \n" );
    double Qinf = adapt_from_a_to_infinite(finf, a, acc, eps, &errinf);
    fprintf(stdout, "The calculation is done. \n" );

    // We calculate the GSL-integration regimes as done in practical programming
    int callsinf2 = 0;
    limit=10000;
    gsl_integration_workspace * winf;
    winf = gsl_integration_workspace_alloc (limit);
    // We setup the funciton again, with the new calls
    double finf2(double x,void * params){callsinf2++; return x/(exp(x)-1);};
    gsl_function  Finf;
    Finf.function = finf2;
    double Qinf2;
    double errinf2;
    gsl_integration_qagiu(&Finf, a, acc, eps, limit, winf, &Qinf2, &errinf2);
    gsl_integration_workspace_free(winf);

    // We print the values
    fprintf(stdout, "The found values are: \n Known solution: π²/6 =  1.6449340668482264364724 \n My solution        I = %25.23g \n  GSL solution      I = %25.23g \n With My found error=%lg  and calls=%i\n and GSL       error=%lg and calls=%i\n Beyond this point of accuracy, the calulation time increases drastically. \n",Qinf,Qinf2,errinf, callsinf,errinf2,callsinf2);
    fprintf(stdout, "We see that both in accuracy and number of function-calls, the GSL-routine is superior. However the implemented transformation works, and provides a rather good accuracy.\n" );
    fprintf(stdout, "With this successfull demonstration of the recursive adaptive integrator using a transformation to handle infinite limits, and a comparison with the GSL-routines, the third exercise is hereby finished.\n This also concludes the examination of adaptive integration. \n\n" );





	return 0;
}
