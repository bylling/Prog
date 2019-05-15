void ode_orbit(double t, gsl_vector* y, gsl_vector*  dydx);
void rkstep12(
  double t,                                  /* the current value of the variable */
	double h,                                  /* the step to be taken */
	gsl_vector* yt,                                /* the current value y(t) of the sought function */
	void f(double t, gsl_vector* y, gsl_vector* dydt), /* the right-hand-side, dydt = f(t,y) */
	gsl_vector* yth,                               /* output: y(t+h) */
	gsl_vector* err                                /* output: error estimate dy */
);


void driver(
	double* t,                             /* the current value of the variable, takes startpoint a and returns endpoint b*/
	double b,                              /* the end-point of the integration */
	double* h,                             /* the current step-size, takes an estimate of the initial size and returns the last accepted stepsize*/
	gsl_vector*yt,                             /* the current y(t), so takes y(a) and returns y(b) */
	double acc,                            /* absolute accuracy goal */
	double eps,                            /* relative accuracy goal */
	void stepper(                          /* the stepper function to be used */
		double t, double h, gsl_vector*yt,
		void f(double t,gsl_vector*y,gsl_vector*dydt),
		gsl_vector*yth, gsl_vector*err
		),
	void f(double t,gsl_vector*y,gsl_vector*dydt) /* right-hand-side */
);

void driver_mod(
	double* t,                             /* the current value of the variable, takes startpoint a and returns endpoint b*/
	double b,                              /* the end-point of the integration */
	double* h,                             /* the current step-size, takes an estimate of the initial size and returns the last accepted stepsize*/
	gsl_vector*yt,                             /* the current y(t), so takes y(a) and returns y(b) */
	double acc,                            /* absolute accuracy goal */
	double eps,                            /* relative accuracy goal */
	void stepper(                          /* the stepper function to be used */
		double t, double h, gsl_vector*yt,
		void f(double t,gsl_vector*y,gsl_vector*dydt),
		gsl_vector*yth, gsl_vector*err
		),
	void f(double t,gsl_vector*y,gsl_vector*dydt), /* right-hand-side */
  gsl_matrix* path,                    /*  The matrix which will be storing the path of the ODE solver. */
  int *pathlength                        /*  An integer to indicate how large a part of the Path-matrix which has been ocupied. */
);

void integration_by_ode(
	double* t,                             /* the current value of the variable, takes startpoint a and returns endpoint b*/
	double b,                              /* the end-point of the integration */
	double* h,                             /* the current step-size, takes an estimate of the initial size and returns the last accepted stepsize*/
	double acc,                            /* absolute accuracy goal */
	double eps,                            /* relative accuracy goal */
	void stepper(                          /* the stepper function to be used */
		double t, double h, gsl_vector*yt,
		void f(double t,gsl_vector*y,gsl_vector*dydt),
		gsl_vector*yth, gsl_vector*err
		),
	double function(double t), /* The function in which we want to integrate.*/
  double* result /* The result, that are send back */
);

double testfunction(double t);
