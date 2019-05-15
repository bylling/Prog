#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include"odefunctions.h"


void ode_orbit(double t, gsl_vector* y, gsl_vector*  dydx){
  // The ordinary differential equaiton from Practical Programming exercise orbit - 2 - III
  double epsilon = 0.01;
  double y_1 = gsl_vector_get(y,0);
  double y_2 = gsl_vector_get(y,1);
  // We return the value at this point of the differential equation
    gsl_vector_set(dydx,0,y_2);
    gsl_vector_set(dydx,1,1-y_1 + epsilon * y_1 * y_1);
}

double testfunction(double t){
  // The testfunction for test by the integrator using ODE is declared.
  // The funciton x^-x from 0 to 1, which is numerically shown if you sum up over all integers to their minus own power will give approx = 1.2913
  double result = pow(t,-t);
  return result;
}

void rkstep12(
  double t,                                  /* the current value of the variable */
	double h,                                  /* the step to be taken */
	gsl_vector* yt,                                /* the current value y(t) of the sought function */
	void f(double t, gsl_vector* y, gsl_vector* dydt), /* the right-hand-side, dydt = f(t,y) */
	gsl_vector* yth,                               /* output: y(t+h) */
	gsl_vector* err                                /* output: error estimate dy */
){
// This function will using Runge-Kutta 1-2 method calculate the value of a differential equation, given a step to be taken. Furthermore it returns an error estimate.
// The function is made as written in the lecture notes, but updated to consider gsl_vectors.
// We call the function to find the value at the function at y(t) with t
f(t, yt, err); // err is used later and therefore acts as k_0 here
// We calculate the next step as yi+1 = yi + h/2 * k0
gsl_vector_scale (err, h/2.0);
gsl_vector_add(yt, err);

//vector_sum(yt, 1., err, h/2); //y = 1.*y + h/2*err

// We find the value of the function at y(t+h/2) to give k_1/2 as 1/2=f(x0+12h,y0+12hk0)
f(t + h/2, yt, yth); //yh being used later and therefore acting as k_1/2 here

// We calculate the backward step as y = yi+1 - h/2 * k0
gsl_vector_scale (err, -1);
gsl_vector_add(yt, err);

// We estimate the error by substraction of the two steps evaluated in stepsize h/2
gsl_vector_scale (err, -1);
gsl_vector_scale (yth, -h/2.0);
gsl_vector_add(err, yth);
gsl_vector_scale (yth, -2.0/h);

// We finally take all of the step yh = k_1/2 * h + yt
gsl_vector_scale (yth, h);
gsl_vector_add(yth,yt);

}


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
){
// The driver function is made similarly to the function in the lecture notes, but with gsl_vectors.

// We initialize the parameters
int step = 0;
double error, normy, tol;
int n = yt->size;
gsl_vector* ystep = gsl_vector_alloc(n);
gsl_vector* err = gsl_vector_alloc(n);
double start = *t;

// We set up a do-while loop to let the driver continiously jump forward until it has reached the goal
do {
  // We check if the function is beyond the point of interrest, then we will progress backwards:
  if(fabs(*t + *h) > fabs(b)) {*h = b - *t;}

  // We call the step funciton with the initial value and current estimate of stepsize and the differential equation to find the new value and a estimate of the error
  stepper(*t,*h,yt,f,ystep,err);
  // We find the norm of the error
  error = gsl_blas_dnrm2(err);
  // We find the norm of the y-function
  normy = gsl_blas_dnrm2(yt);
  // We find the linear tolerance as defined in the lecture notes. τi= (‖yi‖+δ)√hib−a, ei=‖δyi‖
  tol = (normy*eps + acc)*sqrt(*h/(b-start));

  // If we have a smaller error than tolerance we accept the step and continue
  if (tol > error) {
  step++;
  // We update the values and continue
  gsl_vector_memcpy(yt,ystep);
  *t = *t + *h;
  }

  // We estimate the new error
 *h *= pow(tol/error, 0.25)* 0.95;
// We set up the condition for the solver to stop, when we have reached the limit t=b
}while(fabs(*t - b) > 1e-12 && step < 1e6);

gsl_vector_free(ystep);
gsl_vector_free(err);
}



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
){
// The modified driver function is made similarly to the driver function, but by storing every path through the input vector with size given by the input vector length.

// We initialize the parameters
int step = 0;
double error, normy, tol;
int n = yt->size;
gsl_vector* ystep = gsl_vector_alloc(n);
gsl_vector* err = gsl_vector_alloc(n);
double start = *t;

// We set up a do-while loop to let the driver continiously jump forward until it has reached the goal
do {
  // We check if the function is beyond the point of interrest, then we will progress backwards:
  if(fabs(*t + *h) > fabs(b)) {*h = b - *t;}

  // We call the step funciton with the initial value and current estimate of stepsize and the differential equation to find the new value and a estimate of the error
  stepper(*t,*h,yt,f,ystep,err);
  // We find the norm of the error
  error = gsl_blas_dnrm2(err);
  // We find the norm of the y-function
  normy = gsl_blas_dnrm2(yt);
  // We find the linear tolerance as defined in the lecture notes. τi= (‖yi‖+δ)√hib−a, ei=‖δyi‖
  tol = (normy*eps + acc)*sqrt(*h/(b-start));

  // If we have a smaller error than tolerance we accept the step and continue
  if (tol > error) {
  step++;
  // We update the values and continue
  gsl_vector_memcpy(yt,ystep);
  *t = *t + *h;

  // If the step succeds we must store the path of the solver, if there are any more room in the solver we write the following
  if (*pathlength < path->size1) {
    // Then we set the next value in the matrix to the found path
    gsl_matrix_set(path, *pathlength, 0, *t-*h);
    for (int i = 0; i < (yt->size); ++i) {
    gsl_matrix_set(path, *pathlength, i+1, gsl_vector_get(yt, i));

   }
  // Then we remember to note, that the size of the path matrix have grown so:
  (*pathlength)++;
  }

  }

  // We estimate the new error
 *h *= pow(tol/error, 0.25)* 0.95;
// We set up the condition for the solver to stop, when we have reached the limit t=b
}while(fabs(*t - b) > 1e-12 && step < 1e6);

gsl_vector_free(ystep);
gsl_vector_free(err);
}

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
){

// Since we need a function to send into the stepper, we define a one dimensional function, returns the value of the function, as the point-to-point differential for the solution to the integral, since this returns the pont of the function which are integrated.
  void function_in_integral(double t, gsl_vector* y, gsl_vector* functionvalue){
        gsl_vector_set(functionvalue, 0, function(t)); // y'=f(x)
      }

// We allocate the needed parameters:
gsl_vector* yt = gsl_vector_alloc(1); /* the current y(t), so takes y(a) and returns y(b) */
// We set starting value
gsl_vector_set(yt, 0., 0.); // y(a)=0;

// We call the Runge- Kutta 1-2 driver as
driver(t, b, h, yt, acc, eps, stepper,function_in_integral);

// The integral value will now be collected as the end of the ode
*result = gsl_vector_get(yt, 0); // : I = y(b),
// We free the parameters
gsl_vector_free(yt);

}
