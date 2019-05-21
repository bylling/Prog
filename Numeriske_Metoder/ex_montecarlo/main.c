#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"montecarlofunctions.h"


double my_function(double x){return x*x;}
int main (void)
{

  // Exercise A
fprintf(stdout, "Exercise A has started\n\n" );
fprintf(stdout, "A plain monte carlo integrator have been implemented. It estimates an integral between two points and returns both the result and an error estimate. To test it, we use the three test integrals from the previous integration-exercise, since theese are all one dimensional and increase in difficulity due to singularities, they are well suited to test this integrator.\n" );
int N = 1e6;
fprintf(stdout, "All of the integrals will be made with the monte carlo method using %i points.\n", N );

// First integral
fprintf(stdout, "Test 1: I = ∫ from 0 to 1 of √(x) dx = 2/3 = 0.66..\n" );

// We initialze the integration:
int dim1 = 1;
gsl_vector* a1 = gsl_vector_alloc(dim1);
gsl_vector* b1 = gsl_vector_alloc(dim1);
double result1;
double error1;
// Start and end values:
gsl_vector_set(a1,0,0.);
gsl_vector_set(b1,0,1.);
// We setup the funciton
double f(gsl_vector* x){
  double x1 = gsl_vector_get(x,0);
  return sqrt(x1);};
// We calculate the integral
fprintf(stdout, "The calculation is started. \n" );
plainmc(a1, b1, f, N, &result1,&error1);
fprintf(stdout, "The calculation is done. \n" );
// We print the values
fprintf(stdout,"The found values are: integral I=%lg with error=%lg \n\n",result1,error1);

// Second integral
fprintf(stdout, "Test 2: I = ∫ from 0 to 1 of 1/√(x) dx = 2 \n" );

// We initialze the integration
gsl_vector* a2 = gsl_vector_alloc(dim1);
gsl_vector* b2 = gsl_vector_alloc(dim1);
double result2;
double error2;
// Start and end values:
gsl_vector_set(a2,0,0.);
gsl_vector_set(b2,0,1.);
// We setup the funciton
double f2(gsl_vector* x){
  double x1 = gsl_vector_get(x,0);
  return 1./sqrt(x1);};
// We calculate the integral
fprintf(stdout, "The calculation is started. \n" );
plainmc(a2, b2, f2, N, &result2,&error2);
fprintf(stdout, "The calculation is done. \n" );
// We print the values
fprintf(stdout,"The found values are: integral I=%lg with error=%lg \n\n",result2,error2);


// Third integral
fprintf(stdout, "Test 3: I = ∫ from 0 to 1 of ln(x)/√(x) dx = -4 \n" );

// We initialze the integration
gsl_vector* a3 = gsl_vector_alloc(dim1);
gsl_vector* b3 = gsl_vector_alloc(dim1);
double result3;
double error3;
// Start and end values:
gsl_vector_set(a3,0,0.);
gsl_vector_set(b3,0,1.);
// We setup the funciton
double f3(gsl_vector* x){
  double x1 = gsl_vector_get(x,0);
  return log(x1)/sqrt(x1);};
// We calculate the integral
fprintf(stdout, "The calculation is started. \n" );
plainmc(a3, b3, f3, N, &result3,&error3);
fprintf(stdout, "The calculation is done. \n" );
// We print the values
fprintf(stdout,"The found values are: integral I=%lg with error=%lg \n\n",result3,error3);


fprintf(stdout, "We see that as the difficulity of the integral increases, the error also increases. For theese calculations we have held all other parameters constant, but now we are going to challenge the integrator to the fullest, since we will calculate a singular integral in 3 Dimensions. \n" );
int N2 = 1e7;
fprintf(stdout, "This time, we use a higher number of points, since we are dealing with a 3 dimensional integral, we use %i points.\n", N2);

// Calculation
fprintf(stdout, "Calculation: I =∫ from 0 to π  dx/π ∫ from 0 to π  dy/π ∫ form 0 to π  dz/π of [1-cos(x)cos(y)cos(z)]^-1 = Γ(1/4)^4/(4π^3)  \n" );

// We initialze the integration
int dim3 = 3;
gsl_vector* a4 = gsl_vector_alloc(dim3);
gsl_vector* b4 = gsl_vector_alloc(dim3);
double result4;
double error4;
// Start and end values:
// X-values
gsl_vector_set(a4,0,0.);
gsl_vector_set(b4,0,M_PI);

// Y-values
gsl_vector_set(a4,1,0.);
gsl_vector_set(b4,1,M_PI);

// Z-values
gsl_vector_set(a4,2,0.);
gsl_vector_set(b4,2,M_PI);

// We setup the funciton
double f4(gsl_vector* x){
  double x1 = gsl_vector_get(x,0);
  double y = gsl_vector_get(x,1);
  double z = gsl_vector_get(x,2);
  return 1./(M_PI*M_PI*M_PI*(1. - cos(x1) * cos(y) * cos(z)));};
// We calculate the integral
fprintf(stdout, "The calculation is started. \n" );
plainmc(a4, b4, f4, N2, &result4,&error4);
fprintf(stdout, "The calculation is done. \n" );
// We print the values
fprintf(stdout,"The found values are:  \n Known solution:  Γ(1/4)^4/(4π^3) = 1.3932039296856768591842462603255  \n Found solution                 I =%25.23g \n With                       error = %lg \n\n",result4,error4);

fprintf(stdout, "With this successfull demonstration of the monte carlo integrator doing a multidimensional integral, the first exercise is hereby finished.\n\n\n" );

// Exercise B
fprintf(stdout, "Exercise B has started\n\n" );
fprintf(stdout, "To check that the error of the plain Monte Carlo method is of order O(1/sqrt(N)) we choose the particularly simple integral from Test 1, and evaluate it at multiple different values of N, to see the dependency of the error. \n\n" );

// First integral
fprintf(stdout, "We calculate I = ∫ from 0 to 1 of √(x) dx = 2/3 = 0.66..\n" );

// We initialze the integration procedure
gsl_vector* a5 = gsl_vector_alloc(dim1);
gsl_vector* b5 = gsl_vector_alloc(dim1);
double result5;
double error5;
FILE* file = fopen("plotdata1.txt", "w");
// Start and end values:
gsl_vector_set(a5,0,0.);
gsl_vector_set(b5,0,1.);
// We setup the funciton
double f5(gsl_vector* x){
  double x1 = gsl_vector_get(x,0);
  return sqrt(x1);};
// We calculate the integral with different orders of N
fprintf(stdout, "The calculation is started. \n" );
for (double N = 1; N < 7; N+=0.05) {
  plainmc(a5, b5, f5, (int)round(pow(10,N)), &result5,&error5);
  // We print the results to a file
  fprintf(file, "%i %g\n", (int)round(pow(10,N)),fabs(error5));
}
fprintf(stdout, "The calculation is done. \n" );
// We print the values
fprintf(stdout,"The found values from the calculation are plotted in the corresponding figure. The figure is plotted with the absolute error as a funciton of 1/sqrt(N), whereas a linear function will correspond to the expected dependency. \n");
fprintf(stdout, "We see in the figure that for few N the error-dependency varies a bit, but for large N we clearly see the linear dependency in the plot, corresponding to a 1/sqrt(N)-behaviour. \n" );
fprintf(stdout, "With this successfull demonstration of the error-behaviour as a function of N for the plain monte carlo integrator, through solution of a simple integral for different N and correlation of the results in the corresponding figure, the second exercise is hereby finished.\n\n\n" );


// Exercise C
fprintf(stdout, "Exercise C has started.\n\n" );

fprintf(stdout, "A recursive stratified sampling monte carlo integrator have been implemented. It estimates an integral between two points and returns both the result and an error, that will be below the requested error. To test it, we use a combination of the first and second test integrals from the previous integration-exercise, since the first is fairly simple, where the second will be a bit more complicated. This allows the sampling algorithm to focus on the dificult integral-dimension. The combination of the two integrals are herby well suited to test this integrator.\n\n" );

fprintf(stdout, "Test: I = ∫ from 0 to 1 dx ∫ from 0 to 1 dy of √(x)/√(y) = 2 * 2/3 = 4/3 \n" );
// We initialze the integration
double N3 = 60;
fprintf(stdout, "The test is calculated with only %g point at each level of recursion. \n",N3 );

dim3 = 2;
double result6;
double error6;
double acc = 1e-2;
double eps = 1e-2;
gsl_vector* a6 = gsl_vector_alloc(dim3);
gsl_vector* b6 = gsl_vector_alloc(dim3);

// Start and end values:
// X-values
gsl_vector_set(a6,0,0.);
gsl_vector_set(b6,0,1.);

// Y-values
gsl_vector_set(a6,1,0.);
gsl_vector_set(b6,1,1.);

// We setup the function
double f6(gsl_vector* x){
  double x1 = gsl_vector_get(x,0);
  double y = gsl_vector_get(x,1);
  return sqrt(x1)/(sqrt(y));};
// We calculate the integral
fprintf(stdout, "The calculation is started. \n" );
recursivemc(a6, b6, f6, N3, &result6,&error6, acc, eps);
fprintf(stdout, "The calculation is done. \n" );
fprintf(stdout,"The found values are:  \n Known solution:              4/3 = 1.333333  \n Found solution                 I =%25.23g \n With                       error = %lg \n\n",result6,error6);
fprintf(stdout, "With this successfull demonstration of the recursive stratified sampling monte carlo integrator, through solution of a well suited integral, the third exercise is hereby finished.\n" );
fprintf(stdout, "This ends the exercise on Monte Carlo Integration.\n" );
// We free the parameters
gsl_vector_free(a1);
gsl_vector_free(b1);
gsl_vector_free(a2);
gsl_vector_free(b2);
gsl_vector_free(a3);
gsl_vector_free(b3);
gsl_vector_free(a4);
gsl_vector_free(b4);
gsl_vector_free(a5);
gsl_vector_free(b5);
gsl_vector_free(a6);
gsl_vector_free(b6);
fclose(file);
return 0;
}
