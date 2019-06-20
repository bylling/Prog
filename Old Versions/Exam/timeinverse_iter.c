#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"functionfile.h"

int main(int argc, char * argv[]){// This function times the how long the inverse iteration method takes to find the lowest eigenvalue at different dimensions
int systemsize = atof(argv[1]);
gsl_matrix *A=gsl_matrix_alloc(systemsize,systemsize);
gsl_vector *V=gsl_vector_alloc(systemsize);
gsl_vector *e=gsl_vector_alloc(systemsize);
int s = 0;
int update = 100000000;

// We fill the matrix with random real numbers
double randomnumber;
for (int i = 0; i < A->size1; ++i) {
    for (int j = i; j < A->size2; ++j) {
      randomnumber = ((double) rand())/((double)RAND_MAX)*10;
        gsl_matrix_set(A, i, j, randomnumber);
        gsl_matrix_set(A, j, i, randomnumber);
    }
}

//jacobi_mod_ev(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, 1)
inverse_iteration(A, e, V, s, update);

// We free the parameters
gsl_matrix_free(A);
gsl_vector_free(V);
gsl_vector_free(e);
return 0;
}
