#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include"minimizationfunctions.h"
#include"lineqfunctions.h"
#include"rootfunctions.h"
#include"downhillsimplex.h"

int main(int argc, char** argv){

// Exercise A
fprintf(stdout, "Exercise A has started\n" );

fprintf(stdout, "The Newton's method for back-tracking linesearch has been implemented, where the user provides the analytical gradient and the analytical Hessian matrix.\n" );

// We allocate the parameters
gsl_vector* x_rosen = gsl_vector_alloc(2);
gsl_vector* x_himmel = gsl_vector_alloc(2);
gsl_vector* x_rosen2 = gsl_vector_alloc(2);
gsl_vector* x_himmel2 = gsl_vector_alloc(2);
gsl_vector* x_rosen3 = gsl_vector_alloc(2);
gsl_vector* x_himmel3 = gsl_vector_alloc(2);
gsl_vector* functioncall = gsl_vector_alloc(3);
gsl_vector* x_fit = gsl_vector_alloc(3);


// RosenBrock
fprintf(stdout, "We test the implementation by finding the minimum of the Rosenbrock function, \n" );
// We calculate
gsl_vector_set(x_rosen, 0, 0);
gsl_vector_set(x_rosen, 1, 0);
vector_print("We use initial point at x0 = ",x_rosen);
int steprosen1 = newton_with_hessian(&function_rosenbrock_with_H, x_rosen, 1e-6);
fprintf(stdout,"The Newton method used %i steps.\n",steprosen1);
vector_print("The minima is supposed to be at x=[1,1] and found to be x = ",x_rosen);


// Himmelblau
fprintf(stdout, "We also test the implementation by finding the minimum of the Himmelblau's function, \n" );
// We calculate
gsl_vector_set(x_himmel, 0, 2);
gsl_vector_set(x_himmel, 1, 2);
vector_print("We use initial point at x0 = ",x_himmel);
int stephimmel1 = newton_with_hessian(&function_himmel_with_H, x_himmel, 1e-6);
fprintf(stdout,"The Newton method used %i steps.\n",stephimmel1);

vector_print("Solution is supposed to be one of four local minima at f(3.0,2.0)=0.0,f(−2.805118,3.131312)= 0.0, f(−3.779310,−3.283186)=0.0 or f(3.584428,-1.848126)=0.0, and found to be x = ",x_himmel);
fprintf(stdout, "However the function could also have found a local maxima, since the gradient is also zero at maxima. \n" );
fprintf(stdout, "This concludes the first exercise.  \n \n \n" );


// Exercise B
fprintf(stdout, "Exercise B has started\n" );

fprintf(stdout, "A modified Quasi Newton's method for back-tracking linesearch has been implemented using Broyden's update, where the user provides the analytical gradient and the analytical Hessian matrix.\n" );
fprintf(stdout, "Futhermore the root-finding method of the previous exercise have been reimplemented to compare, since root-finding and minimization are the same, when the minima of the function corresponds to a root-value.\n" );
fprintf(stdout, "We test the implementations, by comparing with the results of the previous calculations. All with the same accuracy.\n" );

fprintf(stdout, "Now we find the solution for the minimum of the Rosenbrock funciton.\n");
//For the Rosenbrock we calculate:
gsl_vector_set(x_rosen2, 0, 0);
gsl_vector_set(x_rosen2, 1, 0);
gsl_vector_set(x_rosen3, 0, 0);
gsl_vector_set(x_rosen3, 1, 0);
int steprosen2 = quasi_newton_with_hessian(&function_rosenbrock_without_H, x_rosen2, 1e-6);
gsl_vector_set(functioncall, 1, 0);
int steprosen3 = newton_with_jacobian(&function_rosenbrock_with_j, x_rosen3, 1e-6,functioncall);
fprintf(stdout,"The Quasi Newtons method used %i steps and the root-finding used %i steps, compared to the Newtons method of %i steps. \n",steprosen2, steprosen3,steprosen1);
vector_print("The minima is supposed to be x=[1,1] and found by the Quasi Newton method to be x = ",x_rosen2);
vector_print("Similarly using the root-finding method it is found to be x = ",x_rosen3);

fprintf(stdout, "Now we find the solution for the minimum of the Himmelblau funciton.\n");
//For the Himmelblau function we calculate:
gsl_vector_set(x_himmel2, 0, 2);
gsl_vector_set(x_himmel2, 1, 2);
gsl_vector_set(x_himmel3, 0, 2);
gsl_vector_set(x_himmel3, 1, 2);
int stephimmel2 = quasi_newton_with_hessian(&function_himmel_without_H, x_himmel2, 1e-6);
gsl_vector_set(functioncall,2,0);
int stephimmel3 = newton_with_jacobian(&function_himmel_with_j, x_himmel3, 1e-6,functioncall);
fprintf(stdout,"The Quasi Newtons method used %i steps and the root-finding used %i steps, compared to the Newtons method of %i steps. \n",stephimmel2, stephimmel3, stephimmel1);
vector_print("The minima is supposed to be one of four local minima at f(3.0,2.0)=0.0,f(−2.805118,3.131312)= 0.0, f(−3.779310,−3.283186)=0.0 or f(3.584428,-1.848126)=0.0, and found to be x = ",x_himmel2);
vector_print("Similarly using the root-finding method, we find another solution, which is found to be x = ",x_himmel3);
fprintf(stdout, "By comparing the solutions, we see that the Newtons method is clearly the best, but also uses the fact, that we know the analytical Hessen matrix: The rootfinding method is a bit worse, since it only considers an optimisation in one cross-section at a time, and therefore require more steps, or can find another solution. The Quasi Newton method have more parameters including the stepsize, which can be optimized to the system in case. However since the Hessen matrix is build through iterations at each step, this must converge more slowly. \n" );

fprintf(stdout, "At last we want to use the minimization to solve a non-linear least-squares fitting problem.\n" );
fprintf(stdout, "We read in the data and start the fitting:\n" );
// We read in the data
double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
int N = sizeof(t)/sizeof(t[0]);
// We set up the inital values in a vector as previously
gsl_vector_set(x_fit, 0, 5.0);
gsl_vector_set(x_fit, 1, 2.0);
gsl_vector_set(x_fit, 2, 1.0);
vector_print("We use initial guess at [A,T,B] = ",x_fit);
int stepfit = quasi_newton_with_hessian(&function_fit_without_H, x_fit, 1e-6); // We set the accuracy drasticly down and calculate the fit
fprintf(stdout,"The Quasi Newtons method used %i steps, to converge to a fit. This is the only method in use, since we have only provided the analytical gradient, and not the analytical Hessen matrix. The fact, that the solution is sufficient can later be verified by inspection of the plot. \n",stepfit);
vector_print("The Quasi Newtons Method returned fited [A,T,B] = ",x_fit);

// We know write the exp data into file data1.txt, this is done from the first to the last element in N steps
FILE* file = fopen("data1.txt", "w");
// We print the values:
for(int i = 0; i < N; i += 1) {
    fprintf(file, "%g %g %g\n", t[i], y[i], e[i]);
}
// We know write the fitted data into file fitdata.txt, this is done from the first to the last element in N steps
FILE* file2 = fopen("fitdata.txt", "w");
// We print the values:
double x_min = 0, x_max = 10;
double stepx = 1.00/100.00;
double data;
// Here we find the calculated fit values:
double A = gsl_vector_get(x_fit,0);
double T = gsl_vector_get(x_fit,1);
double B = gsl_vector_get(x_fit,2);
for(double i = x_min; i < x_max; i += stepx) {
    data = A*exp(-i/T)+B;
    fprintf(file2, "%g %g\n", i, data);
}
fprintf(stdout, "The fitted data is printed to a datafile, and the fit with the experimental data is shown in the corresponding figure. \n" );
fprintf(stdout, "This concludes the Second exercise.  \n \n \n" );


// Exercise C
fprintf(stdout, "Exercise C has started\n" );
fprintf(stdout, "A implementation of the Downhill simplex method has been made from the description in the lecture notes. \n" );
fprintf(stdout, "It will now be used to calculate all three previous systems, and its results will be compared to the other methods. For all of the simplex-es, random startvalues are applied with the first quadrant with no numbers larger than 10.\n" );
// We initalize a triangular set of points for our n-dimentional problem.
int dim = 2;
int dim_fit = 3;
// We set the accuracy goal:
double epsilon_simplex = 1e-6;

// We allocate the simplex:
double** simplex=(double**)calloc(dim+1,sizeof(double*));
for (int collumn=0; collumn<dim+1; collumn++) simplex[collumn] = (double*)calloc(dim,sizeof(double));
double** simplex2=(double**)calloc(dim+1,sizeof(double*));
for (int collumn=0; collumn<dim+1; collumn++) simplex2[collumn] = (double*)calloc(dim,sizeof(double));
double** simplex3=(double**)calloc(dim_fit+1,sizeof(double*));
for (int collumn=0; collumn<dim_fit+1; collumn++) simplex3[collumn] = (double*)calloc(dim_fit,sizeof(double));




double randomnumber;
// RosenBrock
fprintf(stdout, "We start by the Rosenbrock function and initialize our simplex.\n" );
// We set the triangular points
for(int n=0; n<dim+1; n++) {
  for (int i=0; i<dim; i++) {
// We pick out a random number as in the eigenvalues exercise
    randomnumber = ((double) rand())/((double)RAND_MAX)*10;
    simplex[n][i] = randomnumber;
  }
}

// we print the starting position
fprintf(stdout, "The random initialized start-simplex will be\n");
for (int n=0; n<dim+1; n++) {
  for (int i=0; i<dim; i++){
     printf(" %g",simplex[n][i]);
   }
  printf("\n");
}
fprintf(stdout, "We start the minimization of the RosenBrock function\n");
// We call the function
downhill_simplex(function_rosenbrock_simple, simplex, dim, epsilon_simplex);

// We find the solutions
fprintf(stdout, "The minima is supposed to be around x=[1,1], and the found minimum simplex is with corners at\n");
for (int n=0; n<dim+1; n++) {
  fprintf(stdout,"[");
  for (int i=0; i<dim; i++){
     fprintf(stdout," %g ",simplex[n][i]);
   }
  fprintf(stdout,"] \n");
}

// Himmelblau
fprintf(stdout, "With that succes, now we find the solution for the minimum of the Himmelblau funciton.\n");
//For the Himmelblau function we calculate:

// We set the triangular points
for(int n=0; n<dim+1; n++) {
  for (int i=0; i<dim; i++) {
// We pick out a random number as in the eigenvalues exercise
    randomnumber = ((double) rand())/((double)RAND_MAX)*10;
    simplex2[n][i] = randomnumber;
  }
}

// we print the starting position
fprintf(stdout, "The random initialized start-simplex will be\n");
for (int n=0; n<dim+1; n++) {
  for (int i=0; i<dim; i++){
     printf(" %g",simplex2[n][i]);
   }
  printf("\n");
}
fprintf(stdout, "We start the minimization of the Himmelblau function\n");

// We call the function
downhill_simplex(function_himmel_simple, simplex2, dim, epsilon_simplex);

// We find the solutions
fprintf(stdout, "The minima is supposed to be to be one of four local minima at f(3.0,2.0)=0.0,f(−2.805118,3.131312)= 0.0, f(−3.779310,−3.283186)=0.0 or f(3.584428,-1.848126)=0.0, and the found minimum simplex is with corners at x = \n");
for (int n=0; n<dim+1; n++) {
  fprintf(stdout,"[");
  for (int i=0; i<dim; i++){
     fprintf(stdout," %g ",simplex2[n][i]);
   }
  fprintf(stdout,"] \n");
}

fprintf(stdout, "Also with that succes, we now find the solution to the minimization for solving a non-linear least-squares fitting problem .\n");
// We set the triangular points
for(int n=0; n<dim_fit+1; n++) {
  for (int i=0; i<dim_fit; i++) {
// We pick out a random number as in the eigenvalues exercise
    randomnumber = ((double) rand())/((double)RAND_MAX)*10;
    simplex3[n][i] = randomnumber;
  }
}

// we print the starting position
fprintf(stdout, "The random initialized start-simplex will be\n");
for (int n=0; n<dim_fit+1; n++) {
  for (int i=0; i<dim_fit; i++){
     printf(" %g",simplex3[n][i]);
   }
  printf("\n");
}
fprintf(stdout, "We start the minimization of the non-linear least-squares fitting problem. \n");
// We call the function
downhill_simplex(function_fit_simple, simplex3, dim_fit, epsilon_simplex);
// We find the solutions
vector_print("The Quasi Newtons Method returned fited [A,T,B] = ",x_fit);

// We find the solutions
fprintf(stdout, "The found minimum simplex is with corners at\n");
for (int n=0; n<dim_fit+1; n++) {
  fprintf(stdout,"[");
  for (int i=0; i<dim_fit; i++){
     fprintf(stdout," %g ",simplex3[n][i]);
   }
  fprintf(stdout,"] \n");
}


fprintf(stdout, "The downhillsimplex-method have hereby been implemented and are able to reproduce every previous result. \n");

fprintf(stdout, "This ends the third and final exercise. \n");
// We free the parameters

gsl_vector_free(x_rosen);
gsl_vector_free(x_himmel);
gsl_vector_free(x_rosen2);
gsl_vector_free(x_himmel2);
gsl_vector_free(x_rosen3);
gsl_vector_free(x_himmel3);
gsl_vector_free(functioncall);
gsl_vector_free(x_fit);
free(simplex);
free(simplex2);
free(simplex3);
 return 0;
}
