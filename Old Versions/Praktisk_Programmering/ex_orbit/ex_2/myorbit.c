#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>

int ode_orbit(double t, const double y[], double dydt[], void *params){
  double epsilon = *(double *) params;
  dydt[0]=y[1];
  dydt[1]= 1 - y[0] + epsilon * y[0] * y[0];
return GSL_SUCCESS;
}

double myorbit(double t, double epsilon, double uprime){
  gsl_odeiv2_system sys;
  sys.function=ode_orbit;
  sys.jacobian=NULL;
  sys.dimension=2;
  sys.params= (void *) &epsilon;
  gsl_odeiv2_driver* driver;
  double hstart=1e-3,abs=1e-6,eps=1e-6;
  driver = gsl_odeiv2_driver_alloc_y_new(
    &sys,
    gsl_odeiv2_step_rkf45,
    hstart,
    abs,
    eps);

    double t0=0;
    double y[2]={1,uprime};
    gsl_odeiv2_driver_apply(driver,&t0,t,y);

    gsl_odeiv2_driver_free(driver);
    return y[0];
}
