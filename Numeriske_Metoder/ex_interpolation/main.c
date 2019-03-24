#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"mylinterp.c"
#include"myinteglin.c"
#include"qspline.c"
#include"cspline.c"
#include <gsl/gsl_spline.h>


double mylinterp(int n, double *x, double *y, double z);
double myinteglin(int n, double *x, double *y, double z);
int main(int argc, char** argv){

  // Exercise 1

// data load

int n=atoi(argv[1]);
double xdata[n],cosdata[n],sindata[n];
 for(int i=0;i<n;i++){
  scanf("%lg %lg %lg",xdata+i,cosdata+i,sindata+i);}

  // Interpolation:

int interpolgridperpoint=20;
int interpoltop=(n)*interpolgridperpoint;
double interpoltopdiv = interpoltop;
double stepsize= (xdata[n-1]-xdata[0])/interpoltopdiv;

double xinter[interpoltop];
double zinter[interpoltop];
double zintegral[interpoltop];

for(int i=0;i<interpoltop;i++){
  xinter[i]=(stepsize*i);
  zinter[i]=mylinterp(n,xdata,cosdata,xinter[i]);
  zintegral[i]=myinteglin(n,xdata,cosdata,xinter[i]);
  fprintf(stderr,"%g %g %g\n",xinter[i],zinter[i],zintegral[i]);}

// Exercise 2
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
// The struckts are allocated for use
qspline *s1 = qspline_alloc(n2, xq, yq1); /* allocates and builds the quadratic spline */
qspline *s2 = qspline_alloc(n2, xq, yq2); /* allocates and builds the quadratic spline */
qspline *s3 = qspline_alloc(n2, xq, yq3); /* allocates and builds the quadratic spline */

// We derine the values for the splines
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

fprintf(file,"%g %g %g %g %g %g %g %g %g %g\n",xqinter[i],zq1inter[i],zq2inter[i],zq3inter[i],zq1interdiff[i],zq2interdiff[i],zq3interdiff[i], zq1interintegrate[i],zq2interintegrate[i],zq3interintegrate[i]);}
//THe functions are freed.
qspline_free(s1); /* free memory allocated in qspline_alloc */
qspline_free(s2); /* free memory allocated in qspline_alloc */
qspline_free(s3); /* free memory allocated in qspline_alloc */

//Exercise 3 will be started now.
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

for(int i=0;i<interpoltop;i++){
  xintercub[i]=(stepsize*i);
  zintercub[i]=cubic_spline_eval(s4,xintercub[i]);
  zintercubdiff[i]=cubic_spline_derivative(s4,xintercub[i]);
  zintercubinteg[i]=cubic_spline_integral(s4,xintercub[i]);
  zintergsl[i] = gsl_spline_eval(csplinegsl, xintercub[i], acc);
 fprintf(file2,"%g %g %g %g %g\n",xintercub[i],zintercub[i],zintercubdiff[i],zintercubinteg[i],zintergsl[i]);
}
gsl_spline_free(csplinegsl);
gsl_interp_accel_free(acc);
cubic_spline_free(s4); /* free memory allocated in qspline_alloc */


fprintf(stdout, "\nExercise 1 is solved with input from a cosine function, which integrated is the sine, the results can be found in the first plot");
fprintf(stdout, "\nExercise 2 is solved, the coefficients are trivial for the three tested functions. \n");
fprintf(stdout, "Exercise 2 results are for a const, lin. and quad. function shown in plot ex2.  \n");
fprintf(stdout, "Exercise 3 results are for a cosine function shown in plot ex3.  \n");
fprintf(stdout, "Here the homemade cspline is well reproduceb by the GSL function.  \n");

 return 0;
}
