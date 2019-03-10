#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"rosen.c"
#include"decayfunction.c"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multimin.h>


int main(int argc, char** argv){
  printf("Exercise 1 Start\n");
  printf("We start at [0,0]\n");
  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  /* Starting point */
  x = gsl_vector_alloc (2);
  gsl_vector_set (x, 0, 0.0);
  gsl_vector_set (x, 1, 0.0);

  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (2);
  gsl_vector_set_all (ss, 1.0);

  /* Initialize method and iterate */
  minex_func.n = 2;
  minex_func.f = rosen;
  minex_func.params = NULL;

  s = gsl_multimin_fminimizer_alloc (T, 2);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status)
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-2);

      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }

      printf ("%5ld %10.3e %10.3e f() = %7.3f size = %.3f\n",
              iter,
              gsl_vector_get (s->x, 0),
              gsl_vector_get (s->x, 1),
              s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < 100);

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);


printf("Exercise 2 Start\n");

if(argc<2){
  fprintf(stderr,"usage: %s number_of_lines_to_read\n",argv[0]);
}
int n=atoi(argv[1]);
double tim[n],sig[n],err[n];
for(int i=0;i<n;i++)
  scanf("%lg %lg %lg",tim+i,sig+i,err+i);

int dim=3;
struct expdata data;
data.n=n;
data.tim=tim;
data.sig=sig;
data.err=err;
gsl_multimin_function F;
F.f = master;
F.n = dim;
F.params=(void*)&data;

gsl_multimin_fminimizer *M;
#define TYPE gsl_multimin_fminimizer_nmsimplex2
M = gsl_multimin_fminimizer_alloc(TYPE,dim);
gsl_vector* start=gsl_vector_alloc(dim);
gsl_vector* step=gsl_vector_alloc(dim);
gsl_vector_set(start,0,6);
gsl_vector_set(start,1,4);
gsl_vector_set(start,2,1);
gsl_vector_set(step,0,0.2);
gsl_vector_set(step,1,0.2);
gsl_vector_set(step,2,0.1);

gsl_multimin_fminimizer_set(M,&F,start,step);

int iter2=0,status2;
double size2;
do{
  iter2++;
  status2 = gsl_multimin_fminimizer_iterate(M);
  if (status2) break;

  size2 = gsl_multimin_fminimizer_size (M);
  status2 = gsl_multimin_test_size (size2, 1e-2);

  if (status2 == GSL_SUCCESS)
      {
        printf ("converged to minimum at\n");
      }

    printf ("iter=%5d A=%g T=%g B=%g master=%g size=%g\n",
            iter2,
            gsl_vector_get (M->x, 0),
            gsl_vector_get (M->x, 1),
            gsl_vector_get (M->x, 2),
            M->fval, size2);
  }
while (status2 == GSL_CONTINUE && iter2 < 100);

double A=gsl_vector_get(M->x,0);
double Tc=gsl_vector_get(M->x,1);
double B=gsl_vector_get(M->x,2);
for(int i=0;i<n;i++) {fprintf(stderr,"%g %g %g %g\n",
  *(tim+i),*(sig+i),*(err+i)
  ,decayfunction(tim[i],A,Tc,B)
);}


printf("Coefficients is found to be A = %g, T= %g and B= %g\n",A,Tc,B);

printf("Exercise 2 fit is plottet in plot.svg\n");

  return 0;
}
