#include<stdlib.h>
#include<assert.h>
#include<stdio.h>
#include<math.h>

typedef struct {int n; double *x,*y,*b,*c,*d;} cubic_spline;

cubic_spline* cubic_spline_alloc(int n,double* x,double* y){//buildscspline written as given in the lecture notes
cubic_spline* s= (cubic_spline*)malloc(sizeof(cubic_spline)); // spline
s->b = (double*)malloc((n-1)*sizeof(double)); // b i
s->c = (double*)malloc((n-1)*sizeof(double)); // c i
s->d = (double*)malloc((n-1)*sizeof(double)); // c i
s->x = (double*)malloc(n*sizeof(double)); // x i
s->y = (double*)malloc(n*sizeof(double));// y i
s->n = n;
for(int i=0; i<n; i++){
  s->x[i]=x[i];
  s->y[i]=y[i];}
double h[n-1],p[n-1]; //VLA from C99
for(int i=0; i <n-1; i ++){
  h[i]=x[i+1]-x[i]; assert(h[i]>0);}
  for(int i=0; i <n-1; i ++){
  p[i]=(y[i+1]-y[i])/h[i];}

double D[n], Q[n-1], B[n];//   b u i l d i n g   t h e   t r i d i a g o n a l   s y s t e m :
D[0]=2;
for(int i=0; i<n-2; i++) D[i+1]=2*h[i]/h[i+1]+2; D[n-1]=2;Q[0]=1;
for(int i=0; i<n-2; i++) Q[i+1]=h[i]/h[i+1];
for(int i=0; i<n-2; i++) B[i+1]=3*(p[i]+p[i+1]*h[i]/h[i+1]);B[0]=3*p[0]; B[n-1]=3*p[n-2];// Gauss   e l i m i n a t i o n   :
for(int i=1; i<n; i++){D[i]-=Q[i-1]/D[i-1];B[i]-=B[i-1]/D[i-1];}
s->b[n-1]=B[n-1]/D[n-1];// back−s u b s t i t u t i o n   :
for(int i=n-2; i>=0;i--) s->b[i]=(B[i]-Q[i]*s->b[i+1])/D[i];
for(int i=0; i<n-1; i++){s->c[i]=(-2*s->b[i]-s->b[i+1]+3*p[i])/h[i];
  s->d[i]=(s->b[i]+s->b[i+1]-2*p[i])/h[i]/h[i];}
return s;
}


double cubic_spline_eval(cubic_spline *s,double z){
// evaluate s(z)
assert(z>=s->x[0] && z<=s->x[s->n-1]);
int i=0, j=s->n-1;
// b i n a r y s e a r c h :
while(j-i>1){
  int m=floor((i+j)/2);
  if(z>s->x[m]) i=m;
  else j=m;
}
double h=z-s->x[i];
return s->y[i]+h*(s->b[i]+h*(s->c[i]+h*s->d[i]));} // inerpolating polynomial  written as given in the lecture notes



void cubic_spline_free(cubic_spline *s) { // free the allocated memory  written as given in the lecture notes
free(s->x);
free(s->y);
free(s->b);
free(s->c);
free(s->d);
free(s);
}


double cubic_spline_derivative(cubic_spline *s, double z){
// Defining the constants
    int n = s->n;
    double *x = s->x;
    double *b = s->b;
    double *c = s->c;
    double *d = s->d;

    // We binary search for derivative,

        int i=0, j=n-1;

        while(j-i>1){
          int m=floor((i+j)/2);
          if(z>s->x[m]) i=m;
          else j=m;
        }

    //calculation of derivative from analytic cubic formula just as previously
    double h = z - x[i];
    return b[i] + 2*c[i]*h + 3*d[i]*h*h;
}


double cubic_spline_integral(cubic_spline *s, double z){
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
//lastly we find the integral up to z
    integral += a[i]*(z - x[i]) + b[i] * pow(z - x[i], 2)/2 + c[i] * pow(z - x[i], 3) / 3 + d[i] * pow(z - x[i], 4) / 4;
    return integral;

}
