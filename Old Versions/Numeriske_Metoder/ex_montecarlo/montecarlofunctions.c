#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"montecarlofunctions.h"


void randomx(gsl_vector* a, gsl_vector* b, gsl_vector* x){
  // We define our random function, as done in table 1 in the lecture notes.
  // This is however modified to take a gsl_vector instead of an array.
  // The function takes two vectors a and b and generates a random vector between them v.
    for (int i = 0; i < a->size; ++i) { // For each dimension of the gsl_vector
      // we load-out the values
        double ai = gsl_vector_get(a, i);
        double bi = gsl_vector_get(b, i);

        double RND = ((double)rand()/RAND_MAX); // Defined random number between 0 and 1

      // we define x as a random superposition of A and B given as X = A - RND*(B-A)
      // X hereby differs from when RND = 0, X = A to RND = 1, X = B with a uniform distribution between theese
        double xi = ai + RND * (bi - ai);
        // We set the x-vector
        gsl_vector_set(x, i, xi);
    }
}

int plainmc(gsl_vector* a, gsl_vector* b, double f(gsl_vector*), int N, double* result, double* error){
  // PLain Monte carlo integrator as defined in the lecture notes. But modified to take gsl_vectors.
  // The function takes the multidimensional points a and b to integrate inbetween. It takes the function in which we want to integrate, the number of points we want to use in the Montecarlo integration procedure. With this, it returns the result and corresponding error of the integration.


  // We define the volume, in which we want to integrate, found as the distance between a and b in all dimensions multiplied
  double V = 1;
  for (int i = 0; i < a->size; ++i) {
    double ai = gsl_vector_get(a, i);
    double bi = gsl_vector_get(b, i);
    V *= bi-ai;
  }
  // We initialize the parameters:
  gsl_vector* x = gsl_vector_alloc(a->size);
  double sum = 0, sumsq = 0;

  // We call the randomx function to generate a random vector in the space between a and b
  // This is done N times, for each point we want to use for integrating.
  for (int i = 0; i < N; ++i) {
    // We fill the x with random values between a and b
    randomx(a, b, x);
    // WE calculate the funciton at each point
    double fx = f(x);
    // We calculate the sum of all functionvalue for the integral solution, and the sum² for the error-estimation.
    sum += fx;
    sumsq += fx*fx;
  }

  // We find the average value and the variance from the sum and sum²
  double avr = sum/N;
  double var = sumsq/N - avr*avr;
  // We return the result as the average of all of the points times the volume we want to integrate
  *result = avr*V;
  // We estibate the error by the central limit theorem,
  *error = sqrt(var/N)*V;

  // We free the parameters
  gsl_vector_free(x);
  return 0;
}


double findvolume(gsl_vector* a, gsl_vector * b) {
// Takes two vectors and returns the volume in between them
// We define the volume, in which we want to integrate, found as the distance between a and b in all dimensions multiplied
double V = 1;
for (int i = 0; i < a->size; ++i) {
    double ai = gsl_vector_get(a, i);
    double bi = gsl_vector_get(b, i);
    V *= bi-ai;
}
return V;
}

int randommontecarlo(gsl_vector* a, gsl_vector* b, double f(gsl_vector*), int N, double * result, double * var){
  // The function takes two vectors a and b, calculate N random points between a and b.
  // Hereafter it returns an average and a variance based on the description in the lecture notes.

// We start by finding the volume between a and b
double V = findvolume(a,b);
double sum, sumsq;
// We define the parameters
gsl_vector* x = gsl_vector_alloc(a->size);

  // We call the randomx function to generate a random vector in the space between a and b
  // This is done N times, for each point we want to use for integrating.
    for (int i = 0; i < N; ++i) {
        // We fill the x with random values between a and b
          randomx(a, b, x);
        // We calculate the funciton at each point
          double fx = f(x);
        // We calculate the sum of all functionvalue for the integral solution, and the sum² for the error-estimation.
          sum += fx;
          sumsq += fx*fx;
    }
  // We find the average value and the variance from the sum and sum²
      double avr = sum/N;
    *var = sumsq/N - avr*avr;
  // We return the result as the average of all of the points times the volume we want to integrate
      *result = avr*V;
  // We estimate the error by the central limit theorem,
//      *error = sqrt(var/N)*V;
// We free the parameters
    gsl_vector_free(x);

  return 0;
}


int recursivemc(gsl_vector* a, gsl_vector* b, double f(gsl_vector*), int N, double* result, double* error, double acc, double eps){
// The recursive Monte carlo integrator as defined in the lecture notes Table 2 & 4.
// The function takes the multidimensional points a and b to integrate inbetween. It takes the function in which we want to integrate, the number of points we want to use in the Montecarlo integration procedure. With this, it returns the result and corresponding error of the integration from the wanted error input.
// If the error do not fullfill the requirement, we subdivide the regime in two regimes, along the axis of most variance and recalculate the interval with the highest error.


// We initialize the parameters:
    gsl_vector* sub_variance_vector1 = gsl_vector_alloc(a->size);
    gsl_vector* sub_variance_vector2 = gsl_vector_alloc(a->size);
    double result0;
    double var0;
    double result1;
    double var1;
    double result2;
    double var2;

// We sample N random  points  with  plain  Monte Carlo  method and estimate  the  average  and the  error
randommontecarlo(a,b,f,N,&result0,&var0);

//We estimate the error by the central limit theorem,
 double V = findvolume(a,b);
 *result = result0;
 *error = sqrt(var0/(double)N)*V;

// If  the  error  is  acceptable
  if (*error<acc+ eps * *result){
// We free the parameters
 gsl_vector_free(sub_variance_vector1);
 gsl_vector_free(sub_variance_vector2);
 // We return the found error and result
return 0;
}
else{
// If the error is not acceptable:

// We allocate the parameters
 gsl_vector* b_avg = gsl_vector_alloc(a->size);
 gsl_vector* a_avg = gsl_vector_alloc(a->size);
 double a_i;
 double test;
 double b_i;
 double ab_average;
// For each dimension
for (int dim = 0; dim < a->size; dim++) {
// We find the midpoint between a and b:
// We read-out the average value between a and b in this dimension
  a_i = gsl_vector_get(a,dim);
  b_i = gsl_vector_get(b,dim);
 ab_average = (b_i + a_i)/2.;
 gsl_vector_memcpy(b_avg,b);
 gsl_vector_memcpy(a_avg,a);
 // We set theese vectors
  gsl_vector_set(a_avg,dim,ab_average);
  gsl_vector_set(b_avg,dim,ab_average);
  // In order to give the gsl_functions some time, to avoid a segmentation fault, we activate a dummy variable test, that doo nothing
  test = 0.;
  a_i += test;
  //We calculate a quick montecarlo for both of the two intervals individually with modified errors
  randommontecarlo(a,b_avg,f,N,&result1,&var1);
  randommontecarlo(a_avg,b,f,N,&result2,&var2);
  // We take the absolute value of the variances
  if (var1 < 0) {
 var1 = - var1;
}
if (var2 < 0) {
 var2 = - var2;
}
 // We store the errors in the vectors
gsl_vector_set(sub_variance_vector1,dim,var1);
gsl_vector_set(sub_variance_vector2,dim,var2);
}
// We add the vectors to find the dimension of maximal variance
gsl_vector_add(sub_variance_vector1,sub_variance_vector2);
int max_variance_index = gsl_vector_max_index(sub_variance_vector1);
// We copy the vectors back again
gsl_vector_memcpy(b_avg,b);
gsl_vector_memcpy(a_avg,a);
// We read-out the average value between a and b in this dimension
 a_i = gsl_vector_get(a,max_variance_index);
 b_i = gsl_vector_get(b,max_variance_index);
 ab_average = (b_i - a_i)/2. + a_i;
// We set theese vectors
gsl_vector_set(a_avg,max_variance_index,ab_average);
gsl_vector_set(b_avg,max_variance_index,ab_average);
// We calculate the two intervals individually with modified errors recursively
recursivemc(a,b_avg,f,N,&result1,&var1,acc/1.414,eps);
recursivemc(a_avg,b,f,N,&result2,&var2,acc/1.414,eps);
//We find the total result and error
*result = result1 + result2;
*error = sqrt(var1 * var1 + var2 * var2) ;

  // We free the parameters
  gsl_vector_free(b_avg);
  gsl_vector_free(a_avg);
  gsl_vector_free(sub_variance_vector1);
  gsl_vector_free(sub_variance_vector2);
  return 0;
  }
}


void vector_print(const char* s, gsl_vector* v){
    printf("%s\n",s);
    for(int i=0;i<v->size;i++){
      printf("%8.3g",gsl_vector_get(v,i));
      printf("\n");
    }
    printf("\n");
  }

void matrix_print(const char* s, gsl_matrix* A){
  printf("%s\n",s);
  for(int i=0;i<A->size1;i++){
    for(int j=0;j<A->size2;j++){
      printf("%8.3g",gsl_matrix_get(A,i,j));
      printf("\t");
    }
  printf("\n");
  }
}
