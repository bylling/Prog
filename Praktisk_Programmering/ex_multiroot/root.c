#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#define TYPE gsl_multiroot_fsolver_hybrids
#define EPS 1e-12

int root_equation
(const gsl_vector * x, void * params, gsl_vector * f)
{
  const double r = gsl_vector_get(x,0);
  const double r2 = gsl_vector_get(x,1);

	double mismatch =  2* ( 200 * r*r*r - 200 * r * r2 + r -1 );
	const double mismatch2 =  200*(r2 - r*r);
	gsl_vector_set(f,0,mismatch);
	gsl_vector_set(f,1,mismatch2);
return GSL_SUCCESS;
}

int root(double * z,double * z2){
	gsl_multiroot_function F;
	F.f=root_equation;
	F.n=2;
	F.params=NULL;

	gsl_multiroot_fsolver * S;
	S = gsl_multiroot_fsolver_alloc(TYPE,F.n);

	gsl_vector* start = gsl_vector_alloc(F.n);

  gsl_vector_set (start, 0, *z);
  gsl_vector_set (start, 1, *z2);

	gsl_multiroot_fsolver_set(S,&F,start);

  fprintf(stdout, "\nThe root finding starts: \n");
	int flag,iter=0;
	do{
		iter++;
		gsl_multiroot_fsolver_iterate(S);
		flag=gsl_multiroot_test_residual(S->f,EPS);
  fprintf(stdout,"x=%g, y=%g iter=%i  gradx=%g grady=%g\n",gsl_vector_get(S->x,0),gsl_vector_get(S->x,1),iter,gsl_vector_get(S->f,0),gsl_vector_get(S->f,1));
	}while(flag==GSL_CONTINUE);
  fprintf(stdout, "\n The root finding has ended \n\n");

  *z = gsl_vector_get(S->x,0);
  *z2 = gsl_vector_get(S->x,1);


  gsl_multiroot_fsolver_free(S);
	gsl_vector_free(start);
	return 0;
}
