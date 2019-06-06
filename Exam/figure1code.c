
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"functionfile.h"
#include"function_eigenval.h"

int main(int argc, char const *argv[]) {

// Here we want to compare the number of iterations each function will use to approximate an eigenvalue and eigenvector for different matrix dimenstions n.

// We initialize the datafile
FILE* file = fopen("plot1.data", "w");

// We sum up over different dimensions
for (int n = 2; n < 31; n++) {

// We initialize
gsl_matrix* A = gsl_matrix_alloc(n,n);
gsl_vector* e = gsl_vector_alloc(n);
gsl_vector* V = gsl_vector_alloc(n);
int iter_pow;
int iter_inv_pow;
int iter_shift_inv = 0;
double s = 0.;
int iter_inv = 0;
int update = 10000;


// We fill the matrix with random real numbers
double randomnumber;
for (int i = 0; i < A->size1; ++i) {
    for (int j = i; j < A->size2; ++j) {
      randomnumber = ((double) rand())/((double)RAND_MAX)*10; // We redefine a random number between 1 and 10
        gsl_matrix_set(A, i, j, randomnumber); // We set the random number
        gsl_matrix_set(A, j, i, randomnumber); // We set the random number
    }
}

// We calculate
iter_pow = power_iteration(A, e, V);
iter_inv_pow = inverse_power_iteration(A, e, V);
iter_shift_inv = shifted_inverse_iteration(A, e, V,s);
iter_inv = inverse_iteration(A, e, V, s, update);

// WE print
fprintf(file, "%i %i %i %i %i\n",n,iter_pow,iter_inv_pow,iter_shift_inv,iter_inv);

// We free the parameters
gsl_matrix_free(A);
gsl_vector_free(e);
gsl_vector_free(V);
}
fclose(file);
return 0;
}
