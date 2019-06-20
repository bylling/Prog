#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"lineqfunctions.h"

double function_rosenbrock_with_H(gsl_vector* x,gsl_vector* y, gsl_matrix* H){ // The Rosenbrock function f(x,y) = (1-x)2+100(y-x2)2 is implemented = 0;
    double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1); // We collect the x-values
		gsl_vector_set(y,0, 2*(x_1-1) - 400*x_1*(x_2-x_1*x_1)); // We Calculate the value of the gradiant in one direction
		gsl_vector_set(y,1, 200*x_2 - 200*x_1*x_1); // We calculate the gradient in the other direction
    double H_11 = 2 - 400*(x_2) + 1200*x_1*x_1; // Analytically Calculated d_X d_X f(x,y)
    double H_12 = -400*x_1; // Analytically Calculated d_X d_y f(x,y)
    double H_21 = -400*x_1; // Analytically Calculated d_y d_X f(x,y)
    double H_22 = 200; // Analytically Calculated d_y d_y f(x,y)
    gsl_matrix_set(H, 0, 0, H_11);
    gsl_matrix_set(H, 0, 1, H_12);
    gsl_matrix_set(H, 1, 0, H_21);
    gsl_matrix_set(H, 1, 1, H_22);
return (1-x_1)*(1-x_1)+100*(x_2-x_1*x_1)*(x_2-x_1*x_1);
    }


double function_himmel_with_H(gsl_vector* x,gsl_vector* y, gsl_matrix* H){ // The Himmelblau function f(x,y) = (x2+y-11)2+(x+y2-7)2 is implemented
		double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1);
    gsl_vector_set(y,0, 4*x_1*(x_1*x_1+x_2-11) + 2*(x_1+x_2*x_2-7)); // We Calculate the gradient in one direction
		gsl_vector_set(y,1, 2*(x_1*x_1+x_2-11) + 4*x_2*(x_1+x_2*x_2-7)); // We calculate the gradient n the other direction

    double H_11 = 12*x_1*x_1 + 4*x_2-42; // Analytically Calculated d_X d_X f(x,y)
    double H_12 = 4*(x_1+x_2); // Analytically Calculated d_X d_y f(x,y)
    double H_21 = 4*(x_1+x_2); // Analytically Calculated d_y d_X f(x,y)
    double H_22 = 4*x_1+12*x_2*x_2-26; // Analytically Calculated d_y d_y f(x,y)
    gsl_matrix_set(H, 0, 0, H_11);
    gsl_matrix_set(H, 0, 1, H_12);
    gsl_matrix_set(H, 1, 0, H_21);
    gsl_matrix_set(H, 1, 1, H_22);
return (x_1*x_1+x_2-11)*(x_1*x_1+x_2-11)+(x_1+x_2*x_2-7)*(x_1+x_2*x_2-7);
    }

double function_rosenbrock_without_H(gsl_vector* x,gsl_vector* y){ // The Rosenbrock function f(x,y) = (1-x)2+100(y-x2)2 is implemented ;
    double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1); // We collect the x-values
    gsl_vector_set(y,0, 2*(x_1-1) - 400*x_1*(x_2-x_1*x_1)); // We Calculate the value of the gradiant in one direction
		gsl_vector_set(y,1, 200*x_2 - 200*x_1*x_1); // We calculate the gradient in the other direction
    return (1-x_1)*(1-x_1)+100*(x_2-x_1*x_1)*(x_2-x_1*x_1);
    	}

double function_himmel_without_H(gsl_vector* x,gsl_vector* y){ // The Himmelblau function f(x,y) = (x2+y-11)2+(x+y2-7)2 is implemented
    double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1); // We collect the x-values
    gsl_vector_set(y,0, 4*x_1*(x_1*x_1+x_2-11) + 2*(x_1+x_2*x_2-7)); // We Calculate the gradient in one direction
		gsl_vector_set(y,1, 2*(x_1*x_1+x_2-11) + 4*x_2*(x_1+x_2*x_2-7)); // We calculate the gradient n the other direction
    return (x_1*x_1+x_2-11)*(x_1*x_1+x_2-11)+(x_1+x_2*x_2-7)*(x_1+x_2*x_2-7);
  }

  double function_rosenbrock_simple(double* x){// The Rosenbrock function f(x,y) = (1-x)2+100(y-x2)2 is implemented
    double x_1=x[0], x_2=x[1]; // We collect the x-values
    return (1-x_1)*(1-x_1)+100*(x_2-x_1*x_1)*(x_2-x_1*x_1);
  }
  double function_himmel_simple(double* x){
    double x_1=x[0], x_2=x[1]; // We collect the x-values
    return (x_1*x_1+x_2-11)*(x_1*x_1+x_2-11)+(x_1+x_2*x_2-7)*(x_1+x_2*x_2-7);
  }

  double function_fit_simple(double* x){ // The Fitting function F(A,T,B)=∑i(f(ti)-yi)²/σi² with f(t)=A*exp(-t/T)+B is implemented
      double A=x[0], T=x[1], B=x[2]; // We collect the x-values

      // We recall the data
      double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
      double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
      double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
      int N = sizeof(t)/sizeof(t[0]);
      // We calculate the fitted function values in theese points:

      double func = 0;
      double data = 0;
      double mismatch = 0;

      for(int i=0;i<N;i++){ // For each point we calculate
      func = A*exp(-t[i]/T)+B; // We calculate the theoretical fitting funciton
      data = y[i]; // We collect the data
      mismatch += pow((func-data)/(e[i]),2); // We find the weighted mismatch, and sum up, since this is the value we want to minimize
      }
    return mismatch;
  }



double function_fit_without_H(gsl_vector* x,gsl_vector* grad){ // The Fitting function F(A,T,B)=∑i(f(ti)-yi)²/σi² with f(t)=A*exp(-t/T)+B is implemented
    double A=gsl_vector_get(x,0), T=gsl_vector_get(x,1), B=gsl_vector_get(x,2); // We collect the x-values

    // We recall the data
    double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
    double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
    double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
    int N = sizeof(t)/sizeof(t[0]);
    // We calculate the fitted function values in theese points:

    double func = 0;
    double data = 0;
    double mismatch = 0;

    for(int i=0;i<N;i++){ // For each point we calculate
    func = A*exp(-t[i]/T)+B; // We calculate the theoretical fitting funciton
    data = y[i]; // We collect the data
    mismatch += pow((func-data)/(e[i]),2); // We find the weighted mismatch, and sum up, since this is the value we want to minimize
    }

    // Since this code is solely writted for a user-provided gradient, we analytically find the gradient of each value and return.
    double difffuncA = 0;
    double difffuncT = 0;
    double difffuncB = 0;

    double diffequationA = 0;
    double diffequationB = 0;
    double diffequationT = 0;

    for(int i=0;i<N;i++){ // For each point we calculate
    func = A*exp(-t[i]/T)+B;  // We calculate the theoretical fitting funciton
    data = y[i]; // We collect the data
    difffuncA = exp(-t[i]/T); // We calculate the theoretical gradient of the fitting funciton
    diffequationA = diffequationA + 2*(func-data)*difffuncA/(e[i]*e[i]); // We calculate the final gradient and sum
    }

    for(int i=0;i<N;i++){ // For each point we calculate
    func = A*exp(-t[i]/T)+B;  // We calculate the theoretical fitting funciton
    data = y[i]; // We collect the data
    difffuncT = A*t[i]*exp(-t[i]/T)/(T*T); // We calculate the theoretical gradient of the fitting funciton
    diffequationT = diffequationT + 2*(func-data)*difffuncT/(e[i]*e[i]); //We calculate the final gradient and sum
    }


    for(int i=0;i<N;i++){ // For each point we calculate
    func = A*exp(-t[i]/T)+B;  // We calculate the theoretical fitting funciton
    data = y[i]; // We collect the data
    difffuncB = 1; // We calculate the theoretical gradient of the fitting funciton
    diffequationB = diffequationB + 2*(func-data)*difffuncB/(e[i]*e[i]); // We calculate the final gradient and sum
    }

    gsl_vector_set(grad,0, diffequationA); // We Calculate the gradient in the first parameter
		gsl_vector_set(grad,1, diffequationT); // We calculate the gradient in the second parameter
    gsl_vector_set(grad,2, diffequationB); // We calculate the gradient in the third parameter

    return mismatch;
  }


int newton_with_hessian(double function(gsl_vector* x, gsl_vector* fx, gsl_matrix* H), gsl_vector* x, double eps){
  //  This function calculates the minimum or maximum by Newton's method with a analytically hessian
  //   Takes in f, which takes the input vector x, calculates analytically the vector f(x), the Hessian matrix, secong arguement is vector x: on input contains the starting point, on output becomes the latest approximation to the minima or maxima;
  //  at last it takes double epsilon: the accuracy goal.

// The function is based on the script from the root - exercise, but modified for the minimization:
// We allocate the needed parameters
    int n = x->size;
    int step = 0;

    gsl_matrix* H = gsl_matrix_alloc(n, n);
    gsl_matrix* R = gsl_matrix_alloc(n, n);
    gsl_vector* grad_f = gsl_vector_alloc(n);
    gsl_vector* Delta_x = gsl_vector_alloc(n);

    function(x, grad_f,H);// We update the Hessen and the gradient
    double norm_grad_f = gsl_blas_dnrm2(grad_f);
// We start the stepping procedure:
    do {
        function(x, grad_f, H); // We update the Hessen and funciton values.
        // Now we need to Solve system grad(phi) + H(x) deltax =0  for stepping size delta x this is done here
        // We solve the system by QR decomposition and multiplying with Q-invese as done in previous exercise on linear equations.
        // First we QR decompose
        qr_gs_decomp(H, R);
        // Then we solve the system of eq. 5 in the lecture notes using backsubstitution
        qr_gs_solve(H, R, grad_f, Delta_x);

        // We find the new x by x + lambda* -deltax
        double lambda = 1.;
        double fx = function(x,grad_f,H);
        gsl_vector_scale(Delta_x, -lambda);
        gsl_vector_add(x, Delta_x);
        gsl_vector_scale(Delta_x, -1.00/lambda);

        // We keep halfing lambda until we satisfy the Armijo condition of eq. 9
        double dotprod;
        gsl_blas_ddot(Delta_x, grad_f, &(dotprod));
        while (function(x,grad_f,H) > (fx + 1e-5*lambda*dotprod) && lambda > 1.00/64.00) {
            lambda /= 2.00;
             // If the lambda is to large, we half it and calculate new values
            gsl_vector_scale(Delta_x, lambda);
            gsl_vector_add(x, Delta_x);
            gsl_vector_scale(Delta_x, 1.00/lambda);
        }
        // When the step  is rescaled we find the new value of the function at the new place
        function(x, grad_f,H);
        norm_grad_f = gsl_blas_dnrm2(grad_f);
        step++;
        // We check if the solution is good enough
    }while(norm_grad_f > eps);

    gsl_matrix_free(H);
    gsl_matrix_free(R);
    gsl_vector_free(Delta_x);
    gsl_vector_free(grad_f);
    // We return the number of steps
    return step;
}

int hessian_broyden_update(gsl_matrix* Hinv, gsl_matrix* H_new, gsl_vector* y, gsl_vector* s, double lambda){ // This function updates the hessian function by use of Broydens update and also return a copy of the needed matrix from the y and s -vectors defined in the lecture notes.
  // We allocate the needed parameters
  int n = Hinv->size1;
  gsl_vector* product = gsl_vector_alloc(n);
  gsl_vector* Hinv_y = gsl_vector_alloc(n);
  gsl_vector* parant = gsl_vector_alloc(n);

    // We make the Broydens  step. Hinv -> Hinv + (s-Hinv y)*s^T*Hinv / (y^T*Hinv s)
    gsl_vector_scale(s, -lambda); // We scale deltax by the step factor lambda to find s

    double y_trans_Hinv_s = 0.;  // We calculate Y^T * Hinv *S
    gsl_blas_dgemv(CblasNoTrans,1.0,Hinv,s,0.0,product); // First Hinv * S
    gsl_blas_ddot(y,product,&(y_trans_Hinv_s)); // Then Y^T Hprod

    // If the update diverges, we reset the inverse Hessian matrix to unity and continue
    if (fabs(y_trans_Hinv_s) < 1e-3) {gsl_matrix_set_identity(Hinv); return 0;}


    gsl_blas_dgemv(CblasNoTrans,1.0,Hinv,y,0.0,Hinv_y); // Hinv * y -> Hinv_y
    gsl_vector_memcpy(parant,s);
    gsl_vector_sub(parant,Hinv_y);  // we find (s-Hinv y)
    gsl_blas_dger(1.00/y_trans_Hinv_s,parant,product,Hinv);//updating the Hinv -> Hinv + (s-Hinv y)*s^T*Hinv / (y^T*Hinv s)
    gsl_matrix_memcpy(H_new, Hinv);

    // We free the used parameters
    gsl_vector_free(product);
    gsl_vector_free(Hinv_y);
    gsl_vector_free(parant);
    return 0;
}


int quasi_newton_with_hessian(double function(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double eps){
  //  This function calculates the minimum or maximum by a quasi Newton's method with a analytically hessian using Broyden's step
  //   Takes in f, which takes the input vector x, calculates analytically the gradient-vector f(x), the Hessian matrix, secong arguement is vector x: on input contains the starting point, on output becomes the latest approximation to the minima or maxima;
  //  at last it takes double epsilon: the accuracy goal.

// The function is based on the previous function, but modified to use Broyden's step thorugh the functioncall to the hessian-update
// We allocate the needed parameters

    double dx = 1e-4;
    int n = x->size;
    int step = 0;

    gsl_matrix* Hinv = gsl_matrix_alloc(n, n);
    gsl_matrix* extra_mem = gsl_matrix_alloc(n, n);
    gsl_vector* grad_f = gsl_vector_alloc(n);
    gsl_vector* grad_f_dx = gsl_vector_alloc(n);
    gsl_vector* Delta_x = gsl_vector_alloc(n);

    // We start by setting the approximate H⁻1 to be I
    gsl_matrix_set_identity(Hinv);
    // We update the gradient
    double fx = function(x, grad_f);

    do {
        fx = function(x,grad_f);// We update the funciton values.
        gsl_blas_dgemv(CblasNoTrans,1.0,Hinv, grad_f,0.0,Delta_x); // We set the step by the approximatied Hessen matrix.
        double lambda = 1.;

        gsl_vector_scale(Delta_x, -lambda);
        gsl_vector_add(x, Delta_x);
        gsl_vector_scale(Delta_x, -1.00/lambda);

        // We keep halfing lambda until we satisfy the Armijo condition of eq. 9
        double dotprod;
        gsl_blas_ddot(Delta_x, grad_f, &(dotprod));
        while (function(x,grad_f_dx) > (fx + 1e-4*lambda*dotprod) && lambda > dx/10) {
            lambda /= 2.00;
            // If the lambda is to large, we half it and calculate new values
            gsl_vector_scale(Delta_x, lambda);
            gsl_vector_add(x, Delta_x);
            gsl_vector_scale(Delta_x, 1.00/lambda);
        }
        // When the step  is rescaled we find the new value of the function at the new place
        function(x, grad_f_dx);
        double norm_grad_f_dx = gsl_blas_dnrm2(grad_f_dx);
        step++;

        // If we satisfy the minimiumcondition, if not, we progress the hessen-matrix by the Broydens step
        if (norm_grad_f_dx < eps) {break;}
        else {
            gsl_vector_scale(grad_f_dx, -1);
            gsl_vector_add(grad_f, grad_f_dx); // We find the change in the gradients
            gsl_vector_scale(grad_f, -1.); // We change the sign of the gradient crange to satisfy conditions for y of eq 12
            hessian_broyden_update(Hinv, extra_mem, grad_f, Delta_x, lambda); // The extern function solves the Broydens step, and will update the Hessian matrix
            gsl_vector_memcpy(grad_f_dx, grad_f);
        }
    }while(step < 1e5);

    gsl_matrix_free(Hinv);
    gsl_matrix_free(extra_mem);
    gsl_vector_free(Delta_x);
    gsl_vector_free(grad_f);
    gsl_vector_free(grad_f_dx);

    return step;
}
