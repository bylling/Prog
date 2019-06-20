
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include"odefunctions.h"
#include"integrationfunctions.h"


int main(int argc, char** argv){

// Exercise A
fprintf(stdout, "Exercise A has started\n" );
fprintf(stdout, "An Embedded Runge-Kutta ODE integrator have been implemented. \n" );
fprintf(stdout, "For a quick demonstration we reproduce the system of the orbit-exercise 2-III from practical programming. \n" );

// We initialize our values
double t,h, acc = 0., eps = 1e-3;
gsl_vector* y = gsl_vector_alloc(2);
FILE* file = fopen("plotdata.txt", "w");
// We start the stepping procedure in a for-loop from the initial to final values, where the ode is solved between each step
fprintf(stdout, "Runge-Kutta 12 has initialized. \n" );
for (double x = 0; x < 50; x+=0.01) {
    // We consider initial values:
        gsl_vector_set(y, 0, 1);
        gsl_vector_set(y, 1, -0.5);
    // With initial stepsize:
        h = 0.01*x/fabs(x);
    // With the current value reset.
        t = 0;

        driver(&t, x, &h, y, acc, eps, &rkstep12, &ode_orbit);

// At each step we print the result.
        fprintf(file, "%g %g %g\n", t, gsl_vector_get(y, 0), gsl_vector_get(y, 1));
    }
    fprintf(stdout, "Runge-Kutta 12 is done. \n" );

fprintf(stdout, "The Runge-Kutta ODE-solver has returned the data for a Relativistic Newtonian Motion with parameters from Orbit 2-III, is plotted in the figure. One can verify that the results are correct by comparing to the Orbit exercise. \n" );
fprintf(stdout, "With a successfull demonstration of the Runge-Kutta 12 ODE integrator and corresponding driver, the first exercise is done.\n \n \n" );

// Exercise B
fprintf(stdout, "Exercise B has started\n" );
fprintf(stdout, "A Modified driver has been made, in which the path is stored in a matrix. This way, we can leave out a for-loop in the mainfile, and call the driver directly. Similar to exercise A, the Relativistic Newtonian Motion will now be calculated using this modified driver.\n" );

// We initialize our values
double b = 50;
int pathlength = 0;
t=0, acc = 0., eps = 1e-3;
gsl_vector* y2 = gsl_vector_alloc(2);
gsl_matrix* ypath = gsl_matrix_alloc(10000,3);
FILE* file2 = fopen("plotdata2.txt", "w");
// We start the stepping procedure in a for-loop from the initial to final values, where the ode is solved between each step
fprintf(stdout, "Modified Runge-Kutta 12 has initialized. \n" );
    // We consider initial values:
        gsl_vector_set(y2, 0, 1);
        gsl_vector_set(y2, 1, -0.5);
    // With initial stepsize:
        h = 0.01*b/fabs(b);
    // With the current value reset.
        t = 0;
       driver_mod(&t, b, &h, y2, acc, eps, &rkstep12, &ode_orbit,ypath,&pathlength);

    // At each we print the result from the matrix
    for (int i = 0; i < pathlength; ++i) {
    fprintf(file2, "%g %g %g\n", gsl_matrix_get(ypath, i,0),  gsl_matrix_get(ypath, i,1),  gsl_matrix_get(ypath, i,2));
    }


    fprintf(stdout, "Modified Runge-Kutta 12 is done. \n" );

fprintf(stdout, "The Modified Runge-Kutta ODE-solver has returned the data in a convenient matrix. The calculation of a Relativistic Newtonian Motion with parameters from Orbit 2-III, is plotted in the corresponding figure. One can verify that the results are correct by comparing to the Orbit exercise, or directly seing that the function coincides with the plotted function from the previous exercise. \n" );
fprintf(stdout, "With a successfull demonstration of the modified driver, the second exercise is also done.\n \n \n" );


// Exercise C
fprintf(stdout, "Exercise C has started\n" );
fprintf(stdout, "A solver for the calculation of definite integrals using the Runge-Kutta 12 ODE-solver from exercise A have been implemented. A test will be made by calculation of the integral of x^-x from 0 to 1, which is known to give approximately 1.291285. This is not a simple function and will provide a challenge for the integration routine. \n" );

// We initialize our values
t=0, acc = 1e-8, eps = 1e-8;
gsl_vector* y3 = gsl_vector_alloc(2);
double result;
FILE* file3 = fopen("plotdata3.txt", "w");
// We start the stepping procedure in a for-loop from the initial to final values, where the ode is solved between each step
fprintf(stdout, "The definite integral will now be calculated \n" );
for (double x = 0.00; x < 1; x+=0.0005) {
    // We consider initial stepsize:
        h = 0.01;
    // With the current value reset.
        t = 0;
        integration_by_ode(&t, x, &h, acc, eps, &rkstep12, testfunction,&result);
        // We print the step by step solution for the indefinete integrated function.

        fprintf(file3, "%g %g %g \n", t, result,testfunction(t));
      }
    fprintf(stdout, "The calculation is done. \n" );
    fprintf(stdout, "The calculation have given the solution to the definite integral on %g. Compared to the known solution on 1.291285. Furthermore the indefinite integral is plotted with the function in the corresponding figure, since the algorithm is also capable of producing that.\n",result);
  fprintf(stdout,"For comparison the same integral is calculated using the adaptive integrator from the integration-exercise.\n" );

// We compare with the found adaptive integrator from the following exerice:

// We initialze the integration:
int calls1 = 0;
double err1 = 0;
double ai=0, bi=1, acci=1e-8, epsi=1e-8;
// We setup the funciton
double f(double x){calls1++; return pow(x,-x);};
// We calculate the integral
fprintf(stdout, "The calculation is started. \n" );
double Q1 = adapt(f, ai, bi, acci, epsi, &err1);
fprintf(stdout, "The calculation is done. \n" );

// We print the values
fprintf(stdout,"The found values from the adaptive integrator are: integral I=%lg, error=%lg, calls=%i\n\n",Q1,err1, calls1);
fprintf(stdout, "The calculation by the ODE-integrator is hereby supported by the findings of the adaptive recursive integration routine. \n" );

fprintf(stdout, "With a successfull demonstration of the integrating ODE-implementation, and comparison with the adaptive recursive integration routine, the third exercise is also done, which completes the investigation on ODE's. \n \n \n" );




// We free the parameters
gsl_vector_free(y);
gsl_vector_free(y2);
gsl_vector_free(y3);
gsl_matrix_free(ypath);
fclose(file);
fclose(file2);
fclose(file3);
return 0;
}
