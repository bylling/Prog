#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>
#include"hydrogen.c"
#define TYPE gsl_multiroot_fsolver_hybrids
#define EPS 1e-6


void root(double* z,double* z1);


int main(){
    printf("Exercise 1\n\n");
  double x = 1.5;
  double y = 1.5;

  double *px = &x;
  double *py = &y;

  printf("Guess for maximum is [%g, %g]\n",*px,*py);
  root(px,py);

    printf("Result [%g,%g], It is Actually = [1,1]\n",*px,*py);

    printf("Exercise 2\n\n");


double xmax = 10;
double estart = -0.5;


gsl_multiroot_function F;
F.f=hydrogen_equation;
F.n=1;
F.params=(void*)&xmax;

gsl_multiroot_fsolver * S;
S = gsl_multiroot_fsolver_alloc(TYPE,F.n);

gsl_vector* start = gsl_vector_alloc(F.n);
gsl_vector_set(start,0,estart);
gsl_multiroot_fsolver_set(S,&F,start);

int flag;
do{
	gsl_multiroot_fsolver_iterate(S);
	flag=gsl_multiroot_test_residual(S->f,EPS);
}while(flag==GSL_CONTINUE);

double result=gsl_vector_get(S->x,0);
gsl_multiroot_fsolver_free(S);
gsl_vector_free(start);
printf("Expected hydrogen energy is -0.5 Hartree\n");
printf("Energy found e=%g\n",result);

printf("We check by inserting e in SE \n");

for (size_t i = 1; i < 100; i++) {
  fprintf(stderr,"%g %g %g\n",xmax*i/100,fe(result,xmax*i/100),xmax*i/100*exp(-xmax*i/100));
}

printf("The solution is plotted with the analytical solution in plot.svg \n");
return 0;
}
