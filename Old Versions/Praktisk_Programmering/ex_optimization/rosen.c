#include<gsl/gsl_vector.h>

double rosen(const gsl_vector *v, void *params)
{
  double x, y;

  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);

  return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
}
