#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include"eigenfunctions.h"

int main(int argc, char * argv[]){// This code calculates the classical and indexed jacobi code, for an arbitraty matrix on random indecies. Made to call and time from the Makefile.
int systemsize = atof(argv[1]);
gsl_matrix *A=gsl_matrix_alloc(systemsize,systemsize);
gsl_matrix *V=gsl_matrix_alloc(systemsize,systemsize);
gsl_vector *e=gsl_vector_alloc(systemsize);
// We fill the matrix with random real numbers
double randomnumber;
for (int i = 0; i < A->size1; ++i) {
    for (int j = i; j < A->size2; ++j) {
      randomnumber = ((double) rand())/((double)RAND_MAX)*10;
        gsl_matrix_set(A, i, j, randomnumber);
        gsl_matrix_set(A, j, i, randomnumber);
    }
}

// We do now initialise the classical sweeping jacobi-method from the lecture notes, but modified with the propper indexing:
int sweeps;
sweeps = jacobi_versc(A,e,V);
sweeps += 0;
// We free the parameters
gsl_matrix_free(A);
gsl_matrix_free(V);
gsl_vector_free(e);
return 0;
}
