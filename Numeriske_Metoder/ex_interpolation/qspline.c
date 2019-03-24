#include<stdlib.h>
#include<assert.h>
#include<math.h>

typedef struct {int n; double *x,*y,*b,*c;} qspline;

qspline* qspline_alloc(int n,double* x,double* y){//buildsqspline written as given in the lecturenotes
qspline* s= (qspline*)malloc(sizeof(qspline)); // spline
s->b = (double*)malloc((n-1)*sizeof(double)); // b i
s->c = (double*)malloc((n-1)*sizeof(double)); // c i
s->x = (double*)malloc(n*sizeof(double)); // x i
s->y = (double*)malloc(n*sizeof(double));// y i
s->n = n;
for(int i=0; i<n; i++){
  s->x[i]=x[i];
  s->y[i]=y[i];}
int i;
double p[n-1],h[n-1]; //VLA from C99
for(i=0; i <n-1; i ++){
  h[i]=x[i+1]-x[i];
  p[i]=(y[i+1]-y[i])/h[i];}
s->c[0]=0;
// r e c u r s i o n up :
for(i=0; i<n-2; i++)s->c[i+1]=(p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];
s->c[n-2]/=2;
// r e c u r s i o n down :
for(i=n-3; i>=0;i--)s->c[i]=(p[i+1]-p[i]-s->c[i+1]*h[i+1])/h[i];
for(i=0; i<n-1;i++)s->b[i]=p[i]-s->c[i]*h[i];
return s;}


double qspline_eval(qspline *s,double z){
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
return s->y[i]+h*(s->b[i]+h*s->c[i]);} // inerpolating polynomial  written as given in the lecturenots



void qspline_free (qspline *s) { // free the allocated memory  written as given in the lecture notes
free(s->x);
free(s->y);
free(s->b);
free(s->c);
free(s);
}

double qspline_derivative(qspline *s, double z){
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

//calculation of derivative from analytic quadratic formula
    return b[i] + 2*c[i]*(z - x[i]);
}



double qspline_integral(qspline *s, double z){
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
