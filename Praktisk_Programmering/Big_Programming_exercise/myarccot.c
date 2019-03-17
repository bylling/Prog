#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>
#include<math.h>

int ode_arccot(double t, const double y[], double dydt[], void* params){
  dydt[0]= -1/(t*t+1);
return GSL_SUCCESS;
}

double myarccot(double t){
	if (t<0) {return -myarccot(-t);};
  gsl_odeiv2_system sys;
  sys.function=ode_arccot;
  sys.jacobian=NULL;
  sys.dimension=1;
  sys.params=NULL;
  gsl_odeiv2_driver* driver;
  double hstart=0.0001,abs=1e-8,eps=1e-8;
  driver = gsl_odeiv2_driver_alloc_y_new(
    &sys,
    gsl_odeiv2_step_rkf45,
    hstart,
    abs,
    eps);

    double t0=0;
    double y[1]={M_PI/2.0};
    gsl_odeiv2_driver_apply(driver,&t0,t,y);

    gsl_odeiv2_driver_free(driver);
    return y[0];
}
