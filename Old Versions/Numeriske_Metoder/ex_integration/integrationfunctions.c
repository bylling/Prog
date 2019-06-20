#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"integrationfunctions.h"

double adapt24(double f(double), double a, double b, double acc, double eps, double *err, double f2, double f3, int nrec)
{ // The Recursive adaptive integrator from the lecture notes, it takes the function in which we want to integrate f, the beginning point a, the end point b, the wanted accuracy acc, and the wanted weighted average eps.
  // It calculates the integration during the step, seperated in four pieces (f(a)) f1, f2, f3 and f4 (and f(bb)), in which the interval is found in sixth-terms of the total interval.
  // The user supplies value f2 and f3, which is a third and two thirds midpoints of the interval, since theese are calculated in the upper layer of the recurse formula, and should be re-used.
  // At last the function takes the number of times the recursive formula has been run, in order to stop it if the wanted errors are too strict.

  // if the function has been called more than a million times, we stop the prodecure
  assert(nrec < 1e6);
  // We calculate the point between f(a) and f2
  double f1 = f(a+(b-a)/6.);
  // We calculate the point between f3 and f(b)
  double f4=f(a+5.*(b-a)/6.);
  // We calculate value of the integral in between a and b, using the trapezium rule for the weights
  double Q = (2*f1 + f2 + f3 + 2*f4)/6. * (b-a);
  // We find the value of the integral using the rectangle rule for the weights, since this is less average than the trapezium rule
  double q=(f1+f2+f3+f4)/4.*(b-a);
  // We calculate the tolerance from the given errors
  double tolerance = acc + eps*fabs(Q);
  // We estimate the error by the difference between the 4' and the 2' method.
  *err = fabs(Q-q);
  // If we can tolerate the error, we return the value of the integral
  if (*err < tolerance) {
    return Q;
  }
  else{ // If not, we need to further part the section in two, in order to call this recursive formula on two smaller sections, to make the errors less and continue on.
    // We call the function on the first half of the integral, and reuse the found points in the interval:
    double Q1 = adapt24(f,a,(a+b)/2.,acc/sqrt(2.),eps,err,f1,f2,nrec+1);
    // We call the function on the second half on the integral, and reuse the found points in the interval:
    double Q2 = adapt24(f,(a+b)/2.,b,acc/sqrt(2.),eps,err,f3,f4,nrec+1);
    // We sum up the values
    return Q1 + Q2;
  }
}

double adapt(double f(double), double a, double b, double acc, double eps, double* err)
{ // This function for adaptiv recursive integration, is based on the one written in the lecture notes.
  // It takes the function in which we want to integrate f, the beginning point a, the end point b, the wanted accuracy acc, and the wanted weighted average eps.
  // It calls the recursive formula, and thereby acts as a initializer for the recursive formula, and returns the funciton value.
  // In order to call the recursive function, we split up the funciton in thirds
  double f2 = f(a+2.*(b-a)/6.);
  double f3 = f(a+4.*(b-a)/6.);
  // We reset the number of recurse calls
  int nrec = 0;

  // We call the recursive function and return the result
  return adapt24(f,a,b,acc,eps, err,f2,f3,nrec);
}

double adapt_clenshaw_curtis(double f(double), double a, double b, double acc, double eps, double* err)
{ // This function for adaptiv recursive integration using a Clenshaw Curtis transformation, is based on equation 55 in the lecture notes,
  // it transfers a integral form -1 to 1 of f(x) dx to the integral from 0 to pi og f(cos(theta))*sin(theta) dtheta
  // This variable transform is done to optimise the integrand to better suit the recursive integration, and can for some function increase accuraCy, while for others reduce accuracy for the same computing time
  // The function take all of the same values as the adapt() function, and do also call the adapt() function, just after the variables have been transformed.

  // We do the transformation to a new nested function based on the input function, and to replace the input function in the further adapt() call.
  // the function transforms x -> cos(theta) , dx -> sin(theta) dtheta
  // This is done with rescaling using a and b:

  double function_clenshaw_curtis(double theta){
    double functioncostheta = f((a+b)/2+(a-b)/2*cos(theta));
    return functioncostheta * sin(theta)*(b-a)/2;
  }
  // We set the new points
  double anew = 0;
  double bnew = M_PI;
  // After the transformation, we call the adaptive algorithm using the new transformed function
  return adapt(function_clenshaw_curtis, anew, bnew, acc, eps, err);
}

double adapt_from_a_to_infinite(double f(double), double a, double acc, double eps, double* err)
{ // This function for adaptiv recursive integration from a to infinity through a coordinate transformation, is based on equation 60 in the lecture notes,
  // it transfers a integral from a to inf of f(x) dx to the integral from 0 to 1 of f(a+(1-t)/t)/t² dt
    // The function take all of the same values as the adapt() function exept the value b, which here is set to infinity, and do also call the adapt() function, just after the variables have been transformed.

  // We do the transformation to a new nested function based on the input function, and to replace the input function in the further adapt() call.
  // the function transforms x -> a+ (1-t)/t , dx -> 1/t² dtheta
  double function_eq60(double t){
    double functiont = f(a + (1-t)/t);
    return functiont * 1/t * 1/t;
  }
  // We set the values points
  double anew = 0;
  double bnew = 1;
  // After the transformation, we call the adaptive algorithm using the new transformed function
  return adapt(function_eq60, anew, bnew, acc, eps, err);
}
