#include<stdlib.h>
#include<assert.h>
#include<stdio.h>
#include<math.h>
#include"interpfunctions.h"



double mylinterp(int n, double * x, double * y, double z){
  // This function do the linear interpolation using binary search. We insert the dimension, x and y-values and the point in which we want to find z.
  // It returns the linear interpolated value at this point z.

  // We initialize
  int i=0, j=n-1;
  // We do the binary search
  while(j-i>1){
    int m=floor((i+j)/2);
    // Until we find the interval in which z is
    if(z>x[m]){
      i = m;    /* We set the lower boarder to the middle*/
    }
    else{
      j = m;     /* We set the upper boarder to the middle*/
    }
  }
double solution = y[i]+(y[i + 1]- y[i])/(x[i+1] - x[i])*(z-x[i]); /* Using lin. interpolation, we find the value within the i and i+1 to return*/
return solution;
}

double myinteglin(int n, double *x, double *y, double z){
  // This function calculates the integral up to this point
  // We call the linear interpolation at the point
    double yofz = mylinterp(n, x, y, z);

    // Linear integration by average square between the ponts
    double integral = 0;
    int i = 0;
    // We integrate every value up to the border before z
    while (x[i+1] < z) {
        integral += 1/2.0 * (y[i + 1] + y[i])*(x[i + 1] - x[i]);
        i++;
    }  // We now have integral up to the point before z *
 // We will now add the last contribution up to z with value S_z
    integral += 1/2.0 * (yofz + y[i])*(z - x[i]);
    return integral;
}

cubic_spline* cubic_spline_alloc(int n,double* x,double* y){//We build the cspline as  given in the lecture notes
cubic_spline* s= (cubic_spline*)malloc(sizeof(cubic_spline)); // spline
s->b = (double*)malloc((n-1)*sizeof(double)); // b i
s->c = (double*)malloc((n-1)*sizeof(double)); // c i
s->d = (double*)malloc((n-1)*sizeof(double)); // c i
s->x = (double*)malloc(n*sizeof(double)); // x i
s->y = (double*)malloc(n*sizeof(double));// y i
s->n = n; // n
// We insert the values in x and y
for(int i=0; i<n; i++){
  s->x[i]=x[i];
  s->y[i]=y[i];}
//We insert the h and p points for the cubic spline
double h[n-1],p[n-1]; //VLA from C99
for(int i=0; i <n-1; i ++){
  // We find the difference between the points
  h[i]=x[i+1]-x[i]; assert(h[i]>0);
}
  // We find the derivative
for(int i=0; i <n-1; i ++){
  p[i]=(y[i+1]-y[i])/h[i];
}

// We initialize the values
double D[n], Q[n-1], B[n];
// We start building the tridiagonal system as done for the cubic spline in the lecture notes
D[0]=2;// Eq. 21 in the lecture notes
for(int i=0; i<n-2; i++) { // Eq. 21 in the lecture notes
  D[i+1]=2*h[i]/h[i+1]+2; D[n-1]=2;Q[0]=1;
}
for(int i=0; i<n-2; i++){ // Eq. 22 in the lecture notes
  Q[i+1]=h[i]/h[i+1];
}
for(int i=0; i<n-2; i++) {// Gauss elimination:
  B[i+1]=3*(p[i]+p[i+1]*h[i]/h[i+1]); // eq. 23
  B[0]=3*p[0];
  B[n-1]=3*p[n-2];
}
for(int i=1; i<n; i++){
  D[i]-=Q[i-1]/D[i-1];B[i]-=B[i-1]/D[i-1]; // eq. 25
}
s->b[n-1]=B[n-1]/D[n-1];// back−substitution:
for(int i=n-2; i>=0;i--) {
  s->b[i]=(B[i]-Q[i]*s->b[i+1])/D[i]; // eq. 27
}
for(int i=0; i<n-1; i++){ // we set the c and d vectors
  s->c[i]=(-2*s->b[i]-s->b[i+1]+3*p[i])/h[i];
  s->d[i]=(s->b[i]+s->b[i+1]-2*p[i])/h[i]/h[i];
}
return s;
}


double cubic_spline_eval(cubic_spline *s,double z){
// evaluate s(z) using the cubic spline method of eq 14 in the lecture notes
int i=0, j=s->n-1;
// First we do a binary search as previouslt :
while(j-i>1){
  int m=floor((i+j)/2);
  if(z>s->x[m]) i=m;
  else j=m;
}
double h=z-s->x[i]; // we find the distance
return s->y[i]+h*(s->b[i]+h*(s->c[i]+h*s->d[i]));} // inerpolating polynomial  written as given in the lecture notes eq. 14



void cubic_spline_free(cubic_spline *s) { // free the allocated memory  written as given in the lecture notes
free(s->x);
free(s->y);
free(s->b);
free(s->c);
free(s->d);
free(s);
}


double cubic_spline_derivative(cubic_spline *s, double z){
  // The derivative function of the cubic spline takes a point and approximates the derivative in the point
// Defining the constants
    int n = s->n;
    double *x = s->x;
    double *b = s->b;
    double *c = s->c;
    double *d = s->d;

    int i=0, j=n-1;

    // We use binary search again,
        while(j-i>1){
          int m=floor((i+j)/2);
          if(z>s->x[m]) i=m;
          else j=m;
        }

    //calculation of derivative from analytic cubic formula just as previously
    double h = z - x[i];
    return b[i] + 2*c[i]*h + 3*d[i]*h*h; // analytic derivatve of a cubic function
}


double cubic_spline_integral(cubic_spline *s, double z){
// This function calculates the integral of the cubic spline up to the point of z

// We take out the data
    double *x = s->x;
    double *a = s->y;
    double *b = s->b;
    double *c = s->c;
    double *d = s->d;
  // As done before, we use the cubic formula for the integral, and sum all sections, just as the linear and quadratic function
    double integral = 0;
    int i = 0;
    while (x[i+1] < z) {
        integral += a[i] * (x[i+1] - x[i]) + b[i] * pow(x[i+1] - x[i], 2) / 2 + c[i] * pow(x[i+1] - x[i], 3) / 3 + d[i] * pow(x[i+1] - x[i], 4) / 4;
        i++;
      }
//lastly we find the integral from the last known point and up to x
    integral += a[i]*(z - x[i]) + b[i] * pow(z - x[i], 2)/2 + c[i] * pow(z - x[i], 3) / 3 + d[i] * pow(z - x[i], 4) / 4;
    return integral;

}



qspline* qspline_alloc(int n,double* x,double* y){//We build a quadratic spline written as given in the lecturenotes
qspline* s= (qspline*)malloc(sizeof(qspline)); // spline
s->b = (double*)malloc((n-1)*sizeof(double)); // b i
s->c = (double*)malloc((n-1)*sizeof(double)); // c i
s->x = (double*)malloc(n*sizeof(double)); // x i
s->y = (double*)malloc(n*sizeof(double));// y i
s->n = n;
for(int i=0; i<n; i++){ // we insert the x and y values
  s->x[i]=x[i];
  s->y[i]=y[i];}
int i;
double p[n-1],h[n-1]; //VLA from C99
for(i=0; i <n-1; i ++){ // We calcualte the linear distances and derivatives
  h[i]=x[i+1]-x[i];
  p[i]=(y[i+1]-y[i])/h[i];}
s->c[0]=0;
//recursion up :
for(i=0; i<n-2; i++){
  s->c[i+1]=(p[i+1]-p[i]-s->c[i]*h[i])/h[i+1]; // eq 11
}
s->c[n-2]/=2;
// recursion down :
for(i=n-3; i>=0;i--){
  s->c[i]=(p[i+1]-p[i]-s->c[i+1]*h[i+1])/h[i]; // eq 12
}
for(i=0; i<n-1;i++){
  s->b[i]=p[i]-s->c[i]*h[i]; // eq 13
}
return s;}


double qspline_eval(qspline *s,double z){
  // This funciton evaluate the quadratic spline in point z
// evaluate s(z)

int i=0, j=s->n-1;
// We start by binary search :
while(j-i>1){
  int m=floor((i+j)/2);
  if(z>s->x[m]) i=m;
  else j=m;
}
double h=z-s->x[i];

return s->y[i]+h*(s->b[i]+h*s->c[i]);} // inerpolating polynomial  written as given in the lecturenotes by eq 7



void qspline_free (qspline *s) { // free the allocated memory  written as given in the lecture notes
free(s->x);
free(s->y);
free(s->b);
free(s->c);
free(s);
}

double qspline_derivative(qspline *s, double z){
  // this function calculates the derivatie of the quadratic spline in the points z
// We take out the data
    int n = s->n;
    double *x = s->x;
    double *b = s->b;
    double *c = s->c;

// We binary search for derivative,
    int i=0, j=n-1;
    while(j-i>1){
      int m=floor((i+j)/2);
      if(z>s->x[m]) i=m;
      else j=m;
    }

//We calculate the derivative from analytic quadratic formula
    return b[i] + 2*c[i]*(z - x[i]);
}



double qspline_integral(qspline *s, double z){
  // This funciton calculates the quadratic spline integral up to a given point z
// We take out the data
    double *x = s->x;
    double *a = s->y;
    double *b = s->b;
    double *c = s->c;

// As done before, we use the quadratic formula for the integral, and sum all sections, just as the linear function
    double integral = 0;
    int i = 0;
    while (x[i+1] < z) {
        integral += a[i] * (x[i+1] - x[i]) + b[i] * pow(x[i+1] - x[i], 2) / 2 + c[i] * pow(x[i+1] - x[i], 3) / 3;
        i++;
    }
//Lastly we find the last section up to z
    integral += a[i]*(z - x[i]) + b[i] * pow(z - x[i], 2)/2 + c[i] * pow(z - x[i], 2)/3;
    return integral;
}
