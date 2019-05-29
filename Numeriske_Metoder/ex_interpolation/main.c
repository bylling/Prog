#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"interpfunctions.h"
#include<gsl/gsl_spline.h>

int main(int argc, char** argv){

  // Exercise A
fprintf(stdout, "Exercise A has started.\n\n" );
fprintf(stdout, "We have implemented a function, that makes the linear spline using binary search from a table at at given point z. Furthermore a function that calculates the integral of this linear spline has been made. As a test a linear interpolation and integration has been made from a dataset of a cosine funciton. \n\n" );
// We load in the data from the input
int n=atoi(argv[1]);
double xdata[n],cosdata[n],sindata[n];
 for(int i=0;i<n;i++){
  scanf("%lg %lg %lg",xdata+i,cosdata+i,sindata+i);}

// We initialize the Interpolation:
FILE* fileA = fopen("plot.data","w");
int interpolgridperpoint=20;
int interpoltop=(n)*interpolgridperpoint;
double interpoltopdiv = interpoltop;
double stepsize= (xdata[n-1]-xdata[0])/interpoltopdiv;
double xinter[interpoltop];
double zinter[interpoltop];
double zintegral[interpoltop];
// We call the interpolation function as well as the integral function
fprintf(stdout, "We start the calculation.\n");
for(int i=0;i<interpoltop;i++){
  xinter[i]=(stepsize*i);
  zinter[i]=mylinterp(n,xdata,cosdata,xinter[i]);
  zintegral[i]=myinteglin(n,xdata,cosdata,xinter[i]);
  // We print the solution
  fprintf(fileA,"%g %g %g\n",xinter[i],zinter[i],zintegral[i]);}
fprintf(stdout, "The found linear spline and integral spline is shown in the corresponding figure.\n");
fprintf(stdout, "With a successfull demonstration of the spline and intragral function the first exercise is hereby done.\n\n\n");

// Exercise B
fprintf(stdout, "Exercise B has started.\n\n" );
fprintf(stdout, "We have implemented a function, that makes the quadratic spline using binary search from a table at at given point z. Furthermore a function that calculates the integral and the derivative of this quadratic spline has been made. As a test a quadratic interpolation, derivation and integration has been made from a dataset of a constant, linear and quadratic funciton. \n\n" );

// The wanted data is made, a simple linear, constant and quadratic regression;
FILE* file1 = fopen("qinitial.data","w");
int n2=5;
double xq[n2],yq1[n2],yq2[n2],yq3[n2];
for (int i = 0; i < n2; i++) {
xq[i]=i+1;
yq1[i]=1;
yq2[i]=i+1;
yq3[i]=(i+1)*(i+1);
fprintf(file1,"%g %g %g %g\n",xq[i],yq1[i],yq2[i],yq3[i]);
}
// The structs are allocated for use
qspline *s1 = qspline_alloc(n2, xq, yq1);
qspline *s2 = qspline_alloc(n2, xq, yq2);
qspline *s3 = qspline_alloc(n2, xq, yq3);

// We initialize the values for the splines
int interpolgridperpoint2=20;
int interpoltop2=(n)*interpolgridperpoint2;
double interpoltopdiv2 = interpoltop2;
double stepsize2= (xq[n2-1]-xq[0])/interpoltopdiv2;
double xqinter[interpoltop2];
double zq1inter[interpoltop2];
double zq2inter[interpoltop2];
double zq3inter[interpoltop2];
double zq1interdiff[interpoltop2];
double zq2interdiff[interpoltop2];
double zq3interdiff[interpoltop2];
double zq1interintegrate[interpoltop2];
double zq2interintegrate[interpoltop2];
double zq3interintegrate[interpoltop2];
FILE* file = fopen("qspline.data","w");
// Here all of the three functions are splined, differentiated and integrated
fprintf(stdout, "We start the calculation.\n");
for(int i=1;i<interpoltop2;i++){
xqinter[i]=(stepsize2*i + 1);
zq1inter[i] = qspline_eval(s1,xqinter[i]);        /* evaluates the prebuilt spline at point z */
zq2inter[i] = qspline_eval(s2,xqinter[i]);        /* evaluates the prebuilt spline at point z */
zq3inter[i] = qspline_eval(s3,xqinter[i]);        /* evaluates the prebuilt spline at point z */
zq1interdiff[i] = qspline_derivative(s1,xqinter[i]);        /* evaluates the prebuilt spline at point z */
zq2interdiff[i] = qspline_derivative(s2,xqinter[i]);        /* evaluates the prebuilt spline at point z */
zq3interdiff[i] = qspline_derivative(s3,xqinter[i]);        /* evaluates the prebuilt spline at point z */
zq1interintegrate[i] = qspline_integral(s1,xqinter[i]);
zq2interintegrate[i] = qspline_integral(s2,xqinter[i]);
zq3interintegrate[i] = qspline_integral(s3,xqinter[i]);
// We print the values.
fprintf(file,"%g %g %g %g %g %g %g %g %g %g\n",xqinter[i],zq1inter[i],zq2inter[i],zq3inter[i],zq1interdiff[i],zq2interdiff[i],zq3interdiff[i], zq1interintegrate[i],zq2interintegrate[i],zq3interintegrate[i]);}
//THe functions are freed.
qspline_free(s1); /* free memory allocated in qspline_alloc */
qspline_free(s2); /* free memory allocated in qspline_alloc */
qspline_free(s3); /* free memory allocated in qspline_alloc */
fprintf(stdout, "The found quadratic spline with derivative and integral is shown in the corresponding figure, and gives the expected trivial values for the three functions.\n");
fprintf(stdout, "With a successfull demonstration of the quadratic spline, derivative and intregral, the second exercise is hereby done.\n\n\n");

//Exercise C will be started now.
fprintf(stdout, "Exercise C has started.\n\n" );
// We initialize the parameters
fprintf(stdout, "We have implemented a function, that makes the cubic spline using binary search from a table at at given point z. Furthermore a function that calculates the integral and the derivative of this cubic spline has been made. As a test a quadratic interpolation, derivation and integration has been made from a dataset of a cosine function, and a comparison is made with the GSL-library functions. \n\n" );
FILE* file2 = fopen("cspline.data","w");
cubic_spline *s4 = cubic_spline_alloc(n,xdata,cosdata); /* allocates and builds the quadratic spline */
double xintercub[interpoltop];
double zintercub[interpoltop];
double zintercubdiff[interpoltop];
double zintercubinteg[interpoltop];
double zintergsl[interpoltop];

//As in the examples, we define the spline
 gsl_interp_accel *acc = gsl_interp_accel_alloc();
 gsl_spline* csplinegsl = gsl_spline_alloc(gsl_interp_cspline, n);
 gsl_spline_init(csplinegsl, xdata, cosdata, n);
fprintf(stdout, "We start the calculation.\n");
// We calculate
for(int i=0;i<interpoltop;i++){
  xintercub[i]=(stepsize*i);
  zintercub[i]=cubic_spline_eval(s4,xintercub[i]);
  zintercubdiff[i]=cubic_spline_derivative(s4,xintercub[i]);
  zintercubinteg[i]=cubic_spline_integral(s4,xintercub[i]);
  zintergsl[i] = gsl_spline_eval(csplinegsl, xintercub[i], acc);
// We print
 fprintf(file2,"%g %g %g %g %g\n",xintercub[i],zintercub[i],zintercubdiff[i],zintercubinteg[i],zintergsl[i]);
}
// We free the parameters
gsl_spline_free(csplinegsl);
gsl_interp_accel_free(acc);
cubic_spline_free(s4); /* free memory allocated in qspline_alloc */
free(file1);
free(file2);
free(fileA);

fprintf(stdout, "The found cubic spline with derivative and integral is shown in the corresponding figure, and gives the expected sine-values for the integral and derivative. We also conclude that the cubic spline are well represented by the GSL-function.\n");
fprintf(stdout, "With a successfull demonstration of the cubic spline, derivative and intregral, the third exercise is hereby done.\n This ends the examination on interpolation.\n\n");

 return 0;
}
