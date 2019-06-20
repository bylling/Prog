#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"functionfile.h"
#include"function_eigenval.h"

int main(int argc, char const *argv[]) {

// Here we want to compare the inverse iteration method for diffrent update parameters at a function of dimension

// We initialize the datafile
FILE* file = fopen("plot2.data", "w");

// We sum up over different dimensions
for (int n = 2; n < 10; n++) {

// We initialize
gsl_matrix* A = gsl_matrix_alloc(n,n);
gsl_vector* e = gsl_vector_alloc(n);
gsl_vector* V = gsl_vector_alloc(n);
int s = 5;
int iter_inv_1 = 0;
int iter_inv_2 = 0;
int iter_inv_3 = 0;
int iter_inv_4 = 0;
int update_1 = 20;
int update_2 = 40;
int update_3 = 60;
int update_4 = 80;


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
iter_inv_1 = inverse_iteration(A, e, V, s, update_1);
iter_inv_2 = inverse_iteration(A, e, V, s, update_2);
iter_inv_3 = inverse_iteration(A, e, V, s, update_3);
iter_inv_4 = inverse_iteration(A, e, V, s, update_4);

// WE print
fprintf(file, "%i %i %i %i %i\n",n,iter_inv_1,iter_inv_2,iter_inv_3,iter_inv_4);

// We free the parameters
gsl_matrix_free(A);
gsl_vector_free(e);
gsl_vector_free(V);
}

fclose(file);
return 0;
}
