#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include"lineqfunctions.h"
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>



void function_linear_equation_with_j(gsl_vector* x,gsl_vector* y, gsl_matrix* J, gsl_vector* functioncall){ // The system of linear equation A*x*y = 1 , exp(-x) + exp(-y) = 1 + 1/A, with A = 10000 is implemented
gsl_vector_set(functioncall,0,gsl_vector_get(functioncall,0)+1);
    double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1);
    double A = 10000;
		gsl_vector_set(y,0, A * x_1*x_2 - 1); // We Calculate the value of the function A*x*y -1 = 0 in the given point
		gsl_vector_set(y,1, exp(-x_1) + exp(-x_2) - 1 - 1.0/A); // We Calculate the value of the function exp(-x) + exp(-y) -( 1 + 1/A)= 0  in the given point

    double J_11 = A * x_2 ; // Analytically Calculated d_X f1(x,y)
    double J_12 = A * x_1; // Analytically Calculated d_y  f1(x,y)
    double J_21 = -exp(-x_1); // Analytically Calculated d_X  f2(x,y)
    double J_22 = -exp(-x_2); // Analytically Calculated d_y f2(x,y)
    gsl_matrix_set(J, 0, 0, J_11);
    gsl_matrix_set(J, 0, 1, J_12);
    gsl_matrix_set(J, 1, 0, J_21);
    gsl_matrix_set(J, 1, 1, J_22);
		}


void function_rosenbrock_with_j(gsl_vector* x,gsl_vector* y, gsl_matrix* J, gsl_vector* functioncall){ // The Rosenbrock function f(x,y) = (1-x)2+100(y-x2)2 is implemented = 0;
gsl_vector_set(functioncall,1,gsl_vector_get(functioncall,1)+1);
    double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1);
		gsl_vector_set(y,0, 2*(1-x_1)*(-1)+100*2*(x_2-x_1*x_1)*(-1)*2*x_1); // We Calculate the value of the function in the given point
		gsl_vector_set(y,1, 100*2*(x_2-x_1*x_1)); // We calculate the gradient of the function
    double J_11 = 2 - 400*(x_2) + 1200*x_1*x_1; // Analytically Calculated d_X d_X f(x,y)
    double J_12 = -400*x_1; // Analytically Calculated d_X d_y f(x,y)
    double J_21 = -400*x_1; // Analytically Calculated d_y d_X f(x,y)
    double J_22 = 200; // Analytically Calculated d_y d_y f(x,y)
    gsl_matrix_set(J, 0, 0, J_11);
    gsl_matrix_set(J, 0, 1, J_12);
    gsl_matrix_set(J, 1, 0, J_21);
    gsl_matrix_set(J, 1, 1, J_22);
		}


void function_himmel_with_j(gsl_vector* x,gsl_vector* y, gsl_matrix* J, gsl_vector* functioncall){ // The Himmelblau function f(x,y) = (x2+y-11)2+(x+y2-7)2 is implemented
gsl_vector_set(functioncall,2,gsl_vector_get(functioncall,2)+1);
		double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1);
    gsl_vector_set(y,0,  4*(x_1*x_1 + x_2 - 11)*x_1 + 2*(x_1 + x_2*x_2 - 7)); // We Calculate the value of the function in the given point
		gsl_vector_set(y,1, 2*(x_1*x_1 + x_2 - 11) + 4*(x_1 + x_2*x_2 - 7)*x_2); // We calculate the gradient of the function

    double J_11 = 12*x_1*x_1 + 4*x_2-42; // Analytically Calculated d_X d_X f(x,y)
    double J_12 = 4*(x_1+x_2); // Analytically Calculated d_X d_y f(x,y)
    double J_21 = 4*(x_1+x_2); // Analytically Calculated d_y d_X f(x,y)
    double J_22 = 4*x_1+12*x_2*x_2-26; // Analytically Calculated d_y d_y f(x,y)
    gsl_matrix_set(J, 0, 0, J_11);
    gsl_matrix_set(J, 0, 1, J_12);
    gsl_matrix_set(J, 1, 0, J_21);
    gsl_matrix_set(J, 1, 1, J_22);
		}

void function_linear_equation(gsl_vector* x,gsl_vector* y, gsl_vector* functioncall){ // The system of linear equation A*x*y = 1 , exp(-x) + exp(-y) = 1 + 1/A, with A = 10000 is implemented
gsl_vector_set(functioncall,3,gsl_vector_get(functioncall,3)+1);
    double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1);
    double A = 10000;
		gsl_vector_set(y,0, A * x_1*x_2 - 1); // We Calculate the value of the function A*x*y -1 = 0 in the given point
		gsl_vector_set(y,1, exp(-x_1) + exp(-x_2) - 1 - 1.0/A); // We Calculate the value of the function exp(-x) + exp(-y) -( 1 + 1/A)= 0  in the given point
	}

void function_rosenbrock(gsl_vector* x,gsl_vector* y,  gsl_vector* functioncall){ // The Rosenbrock function f(x,y) = (1-x)2+100(y-x2)2 is implemented = 0;
gsl_vector_set(functioncall,4,gsl_vector_get(functioncall,4)+1);
    double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1);
		gsl_vector_set(y,0, 2*(1-x_1)*(-1)+100*2*(x_2-x_1*x_1)*(-1)*2*x_1); // We Calculate the value of the function in the given point
		gsl_vector_set(y,1, 100*2*(x_2-x_1*x_1)); // We calculate the gradient of the function
	}

void function_himmel(gsl_vector* x,gsl_vector* y, gsl_vector* functioncall){ // The Himmelblau function f(x,y) = (x2+y-11)2+(x+y2-7)2 is implemented
gsl_vector_set(functioncall,5,gsl_vector_get(functioncall,5)+1);
		double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1);
    gsl_vector_set(y,0,  4*(x_1*x_1 + x_2 - 11)*x_1 + 2*(x_1 + x_2*x_2 - 7)); // We Calculate the value of the function in the given point
		gsl_vector_set(y,1, 2*(x_1*x_1 + x_2 - 11) + 4*(x_1 + x_2*x_2 - 7)*x_2); // We calculate the gradient of the function
	}


int newton_with_jacobian(void function(gsl_vector* x, gsl_vector* fx, gsl_matrix* J, gsl_vector* functioncall), gsl_vector* x, double epsilon,  gsl_vector* functioncall){
    //  This function calculates the root by Newton's method with a analytically jacobian
    // Takes in f: takes the input vector x, calculates analytically the vector f(x), the jacobian matrix, secong arguement is vector x: on input contains the starting point, on output becomes the latest approximation to the root;
    //  at last it takes double epsilon: the accuracy goal.

// The function is based on the python script example in the lecture notes:
// We allocate the needed parameters
int number_of_eq=x->size;
gsl_matrix* R = gsl_matrix_alloc(number_of_eq, number_of_eq);
gsl_matrix* J = gsl_matrix_alloc(number_of_eq,number_of_eq);
gsl_vector* fx = gsl_vector_alloc(number_of_eq);

gsl_vector* deltax = gsl_vector_alloc(number_of_eq);
gsl_vector* new_fx = gsl_vector_alloc(number_of_eq);
gsl_vector* cfx = gsl_vector_alloc(number_of_eq);
double lambda = 1.00;
int step=0;
double functionnorm, newfunctionnorm;
// We start the stepping procedure:
do {
  step++;
function(x,fx,J,functioncall); // We update the Jacobian and funciton values.
gsl_vector_memcpy(cfx, fx);
// Now we need to Solve system J∆x=−f(x) for stepping size delta x this is done through QR decomposition and multiplying with Q-invese as done in previous exercise on linear equations.
// First we QR decompose
 qr_gs_decomp(J, R);
 // Then we solve the system using backsubstitution
 qr_gs_solve(J, R, cfx, deltax);
 lambda = 1.00;
 // We find the new x by x + lambda* -deltax
gsl_vector_scale (deltax, -lambda);
gsl_vector_add(x, deltax);

// we evaluate the function at the new x
function(x,new_fx,J,functioncall);
// We find the norm of the function, being the distance to 0.0 for both the old and new function
functionnorm = gsl_blas_dnrm2(fx);
newfunctionnorm = gsl_blas_dnrm2(new_fx);

// We find the propper lambda factor by evluating the criteria
 while (newfunctionnorm > (1-lambda/2.00)*functionnorm && lambda > 1.00/64.00) {
   // If the lambda is to large, we half it and calculate new values
   lambda /= 2.00;
   gsl_vector_scale(deltax, lambda);
   gsl_vector_add(x, deltax);
   // we evaluate the function at the new x
   function(x,new_fx,J,functioncall);
   newfunctionnorm = gsl_blas_dnrm2(new_fx);
 }

// When lambda is rescaled we find the new value of the function at the new place
function(x,new_fx,J,functioncall);
newfunctionnorm = gsl_blas_dnrm2(new_fx);
}while(newfunctionnorm > epsilon);

// We free the parameters
gsl_matrix_free(R);
gsl_matrix_free(J);
gsl_vector_free(fx);
gsl_vector_free(deltax);
gsl_vector_free(cfx);
gsl_vector_free(new_fx);
return step;
}


void numjacobian(void function(gsl_vector* x, gsl_vector* fx,gsl_vector* functioncall), gsl_matrix* J, gsl_vector* x, gsl_vector* fx, gsl_vector* fx_dx, double dx,gsl_vector* functioncall){
// Calculates the jacobian of the function numerically by the formulas in the lecture notes, evaluated in steps of dx
for (int i = 0; i < x->size; ++i) {
  double x_i = gsl_vector_get(x, i); // We pull out an x value from the vector
  gsl_vector_set(x, i, x_i + dx); // In its place we set the value plus the stepsize
  function(x, fx_dx,functioncall);  // We evaluate the function
  gsl_vector_set(x, i, x_i);  // We reset the x-value to the initial point

  for (int j = 0; j < x->size; ++j) {
    double df_jdx_i = (gsl_vector_get(fx_dx, j) - gsl_vector_get(fx, j)) / dx; // We find the numerical jacobian by numerical integration of deltaf/dx in every dimension
    gsl_matrix_set(J, j, i, df_jdx_i); // we return this jacobian in the matrix
    }
  }
}



int newton(void function(gsl_vector* x, gsl_vector* fx, gsl_vector* functioncall), gsl_vector* x, double dx, double epsilon,gsl_vector* functioncall){
//  This function calculates the root by Newton's method with a numerical jacobian
// Takes in f: takes the input vector x, calculates analytically the vector f(x), the jacobian is calculated within here, secong arguement is vector x: on input contains the starting point, on output becomes the latest approximation to the root;
//  at last it takes double epsilon: the accuracy goal.

// The function is completely copy-paste of the previous function, unless that at eny point of the stepping procedure, one would need to find the numerical jacobian once again
// We allocate the needed parameters
int number_of_eq=x->size;
gsl_matrix* R = gsl_matrix_alloc(number_of_eq, number_of_eq);
gsl_matrix* J = gsl_matrix_alloc(number_of_eq,number_of_eq);
gsl_vector* fx = gsl_vector_alloc(number_of_eq);

gsl_vector* deltax = gsl_vector_alloc(number_of_eq);
gsl_vector* new_fx = gsl_vector_alloc(number_of_eq);
gsl_vector* cfx = gsl_vector_alloc(number_of_eq);
double lambda = 1.00;
int step=0;
double functionnorm, newfunctionnorm;
// We start the stepping procedure:
do {
  step++;
  function(x,fx,functioncall);
numjacobian(function, J,  x, fx, new_fx, dx,functioncall);; // We update the Jacobian and funciton values.
gsl_vector_memcpy(cfx, fx);
// Now we need to Solve system J∆x=−f(x) for stepping size delta x this is done through QR decomposition and multiplying with Q-invese as done in previous exercise on linear equations.
// First we QR decompose
 qr_gs_decomp(J, R);
 // Then we solve the system using backsubstitution
 qr_gs_solve(J, R, cfx, deltax);
 lambda = 1.00;
 // We find the new x by x + lambda* -deltax
gsl_vector_scale (deltax, -lambda);
gsl_vector_add(x, deltax);

// we evaluate the function at the new x
function(x,new_fx,functioncall);
// We find the norm of the function, being the distance to 0.0 for both the old and new function
functionnorm = gsl_blas_dnrm2(fx);
newfunctionnorm = gsl_blas_dnrm2(new_fx);

// We find the propper lambda factor by evluating the criteria
 while (newfunctionnorm > (1-lambda/2.00)*functionnorm && lambda > 1.00/64.00) {
   // If the lambda is to large, we half it and calculate new values
   lambda /= 2.00;
   gsl_vector_scale(deltax, lambda);
   gsl_vector_add(x, deltax);
   // we evaluate the function at the new x
   function(x,new_fx,functioncall);
   newfunctionnorm = gsl_blas_dnrm2(new_fx);
 }

// When lambda is rescaled we find the new value of the function at the new place
function(x,new_fx,functioncall);
newfunctionnorm = gsl_blas_dnrm2(new_fx);
}while(newfunctionnorm > epsilon);

// We free the parameters
gsl_matrix_free(R);
gsl_matrix_free(J);
gsl_vector_free(fx);
gsl_vector_free(deltax);
gsl_vector_free(cfx);
gsl_vector_free(new_fx);
return step;
}

// The following is directly a modification of the multiroot procedure as in the multiroot exercise for Prktiskm Programmering.
int gsl_root_equation_lin(const gsl_vector * x, void * params, gsl_vector * f)
{
      double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1);
      double A = 10000;
  		gsl_vector_set(f,0, A * x_1*x_2 - 1); // We Calculate the value of the function A*x*y -1 = 0 in the given point
  		gsl_vector_set(f,1, exp(-x_1) + exp(-x_2) - 1 - 1.0/A); // We Calculate the value of the function exp(-x) + exp(-y) -( 1 + 1/A)= 0  in the given point
  return GSL_SUCCESS;
}

int gsl_root_equation_rosen(const gsl_vector * x, void * params, gsl_vector * f)
{
      double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1);
  		gsl_vector_set(f,0, 2*(1-x_1)*(-1)+100*2*(x_2-x_1*x_1)*(-1)*2*x_1); // We Calculate the value of the function in the given point
  		gsl_vector_set(f,1, 100*2*(x_2-x_1*x_1)); // We calculate the gradient of the function
return GSL_SUCCESS;
}

int gsl_root_equation_himmel(const gsl_vector * x, void * params, gsl_vector * f)
{
  		double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1);
      gsl_vector_set(f,0,  4*(x_1*x_1 + x_2 - 11)*x_1 + 2*(x_1 + x_2*x_2 - 7)); // We Calculate the value of the function in the given point
  		gsl_vector_set(f,1, 2*(x_1*x_1 + x_2 - 11) + 4*(x_1 + x_2*x_2 - 7)*x_2); // We calculate the gradient of the function
return GSL_SUCCESS;
}



int gslroot_lin(double * z,double * z2){
  double EPS = 1e-6; // We use the same epsilon as before.
	gsl_multiroot_function F;
	F.f=gsl_root_equation_lin;
	F.n=2;
	F.params=NULL;

	gsl_multiroot_fsolver * S;
	S = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids,F.n);

	gsl_vector* start = gsl_vector_alloc(F.n);

  gsl_vector_set (start, 0, *z); // We use the same start-values as before
  gsl_vector_set (start, 1, *z2);

	gsl_multiroot_fsolver_set(S,&F,start);

  fprintf(stdout, "\nThe root finding starts: \n");
	int flag,iter=0;
	do{
		iter++;
		gsl_multiroot_fsolver_iterate(S);
		flag=gsl_multiroot_test_residual(S->f,EPS);
  //fprintf(stdout,"x=%g, y=%g iter=%i  gradx=%g grady=%g\n",gsl_vector_get(S->x,0),gsl_vector_get(S->x,1),iter,gsl_vector_get(S->f,0),gsl_vector_get(S->f,1));
	}while(flag==GSL_CONTINUE);

  *z = gsl_vector_get(S->x,0);
  *z2 = gsl_vector_get(S->x,1);

  fprintf(stdout, "\n The GSL_root finding for the linear function has ended, in %i steps at point [%g,%g] \n",iter,gsl_vector_get(S->x,0),gsl_vector_get(S->x,1));

  gsl_multiroot_fsolver_free(S);
	gsl_vector_free(start);
	return 0;
}



int gslroot_rosen(double * z,double * z2){
  double EPS = 1e-6; // We use the same epsilon as before.
	gsl_multiroot_function F;
	F.f=gsl_root_equation_rosen;
	F.n=2;
	F.params=NULL;

	gsl_multiroot_fsolver * S;
	S = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids,F.n);

	gsl_vector* start = gsl_vector_alloc(F.n);

  gsl_vector_set (start, 0, *z); // We use the same start-values as before
  gsl_vector_set (start, 0, *z2);

	gsl_multiroot_fsolver_set(S,&F,start);

  fprintf(stdout, "\nThe root finding starts: \n");
	int flag,iter=0;
	do{
		iter++;
		gsl_multiroot_fsolver_iterate(S);
		flag=gsl_multiroot_test_residual(S->f,EPS);
  //fprintf(stdout,"x=%g, y=%g iter=%i  gradx=%g grady=%g\n",gsl_vector_get(S->x,0),gsl_vector_get(S->x,1),iter,gsl_vector_get(S->f,0),gsl_vector_get(S->f,1));
	}while(flag==GSL_CONTINUE);
  //fprintf(stdout, "\n The root finding has ended \n\n");

  *z = gsl_vector_get(S->x,0);
  *z2 = gsl_vector_get(S->x,1);
  fprintf(stdout, "\n The GSL_root finding for the Rosenbrock function has ended, in %i steps at point [%g,%g] \n",iter,gsl_vector_get(S->x,0),gsl_vector_get(S->x,1));


  gsl_multiroot_fsolver_free(S);
	gsl_vector_free(start);
	return 0;
}

int gslroot_himmel(double * z,double * z2){
  double EPS = 1e-6; // We use the same epsilon as before.
	gsl_multiroot_function F;
	F.f=gsl_root_equation_himmel;
	F.n=2;
	F.params=NULL;

	gsl_multiroot_fsolver * S;
	S = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids,F.n);

	gsl_vector* start = gsl_vector_alloc(F.n);

  gsl_vector_set (start, 1, *z); // We use the same start-values as before
  gsl_vector_set (start, 1, *z2);

	gsl_multiroot_fsolver_set(S,&F,start);

  fprintf(stdout, "\nThe root finding starts: \n");
	int flag,iter=0;
	do{
		iter++;
		gsl_multiroot_fsolver_iterate(S);
		flag=gsl_multiroot_test_residual(S->f,EPS);
  //fprintf(stdout,"x=%g, y=%g iter=%i  gradx=%g grady=%g\n",gsl_vector_get(S->x,0),gsl_vector_get(S->x,1),iter,gsl_vector_get(S->f,0),gsl_vector_get(S->f,1));
	}while(flag==GSL_CONTINUE);
  //fprintf(stdout, "\n The root finding has ended \n\n");

  *z = gsl_vector_get(S->x,0);
  *z2 = gsl_vector_get(S->x,1);
  fprintf(stdout, "\n The GSL_root finding for the Himmelblau function has ended, in %i steps at point [%g,%g] \n",iter,gsl_vector_get(S->x,0),gsl_vector_get(S->x,1));


  gsl_multiroot_fsolver_free(S);
	gsl_vector_free(start);
	return 0;
}
