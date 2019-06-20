#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_odeiv2.h>
#include<math.h>
#include<assert.h>
int ode_H(double x, const double y[], double dydx[], void *params)
{
	double e=*(double*)params;
	dydx[0] = y[1];
	dydx[1] = 2*(-1/x-e)*y[0];
	return GSL_SUCCESS;
}

double fe(double e,double xmax)
{ assert(xmax>=0);
	gsl_odeiv2_system sys;
	sys.function = ode_H;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = (void*)&e;

	gsl_odeiv2_driver *driver;
	double hstart = 0.001, abs = 1e-8, eps = 1e-8;
	driver = gsl_odeiv2_driver_alloc_y_new(&sys,
					       gsl_odeiv2_step_rkf45,
					       hstart, abs, eps);

	double x0 = 1e-4;
	double y[] = { x0-x0*x0, 1-2*x0};
	gsl_odeiv2_driver_apply(driver, &x0, xmax, y);
	gsl_odeiv2_driver_free(driver);
	return y[0];
}

int hydrogen_equation
(const gsl_vector * x, void * params, gsl_vector * f)
{
	double xmax = *(double*)params;
	double e = gsl_vector_get(x,0);
	double mismatch = fe(e,xmax);
	gsl_vector_set(f,0,mismatch);
return GSL_SUCCESS;
}
