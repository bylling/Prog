#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"functionfile.h"
#include"function_eigenval.h"
int main (void)
{
  // Initial text:
fprintf(stdout, "\n\n ------------- The Final Exam -----------------\n" );
fprintf(stdout, "Exercise 10: Inverse iteration algorithm for eigenvalues and eigenvectors. \n\n" );
fprintf(stdout, "Implement the variant of the inverse iteration method that calculates the eigenvalue closest to a given number s (and the corresponding eigenvector).\n\n" );

fprintf(stdout, "The Exercise will now be made, and the analysis of the method are described in the report.pdf file.\n\n" );
fprintf(stdout, "\n\n ------------- Quick Test of the Algorithm -----------------\n" );
fprintf(stdout, "As a test, to check whether the implementation works, an example is hereby made using a random real symmetric matrix.\n" );
fprintf(stdout, "The reason for restricting this test to a symmetric matrix is the fact that any real symmetric matrix is Hermitian and therefore all of its eigenvalies are uniquely determined and are real.\n" );

// We initialize our values
int n = 3;
gsl_matrix* A = gsl_matrix_alloc(n,n);
gsl_matrix* A_classic = gsl_matrix_alloc(n,n);
gsl_matrix* D = gsl_matrix_alloc(n,n);
gsl_matrix* Product = gsl_matrix_alloc(n,n);
gsl_matrix* V_classic = gsl_matrix_alloc(n,n);
gsl_vector* V_power = gsl_vector_alloc(n);
gsl_vector* V_inv_power = gsl_vector_alloc(n);
gsl_vector* V_shift_inv = gsl_vector_alloc(n);
gsl_vector* V_inv = gsl_vector_alloc(n);
gsl_vector* e_classic = gsl_vector_alloc(n);
gsl_vector* e_power = gsl_vector_alloc(n);
gsl_vector* e_inv_power = gsl_vector_alloc(n);
gsl_vector* e_shift_inv = gsl_vector_alloc(n);
gsl_vector* e_inv = gsl_vector_alloc(n);
// We fill the matrix with random real numbers
double randomnumber;
for (int i = 0; i < A->size1; ++i) {
    for (int j = i; j < A->size2; ++j) {
      randomnumber = ((double) rand())/((double)RAND_MAX)*10; // We redefine a random number between 1 and 10
        gsl_matrix_set(A, i, j, randomnumber); // We set the random number
        gsl_matrix_set(A, j, i, randomnumber); // We set the random number
    }
}
// We take a copy for the comparison with the classic Jacobi eigenvalue method.
gsl_matrix_memcpy(A_classic,A);


// Now we print the Test matrix:
fprintf(stdout, "The random initial symmetric matrix with indices between 0 and 10 is:\n");
matrix_print("A=",A);

fprintf(stdout, "\n\n ------------- Previously made Classic Jacobi Method for Comparison -----------------\n" );

fprintf(stdout, "The 'Classic' Jacobi eigenvalue method have been implemented as used in exercise C of the eigenvalue exercise, and will be used as a comparison. Thereby the recieved eigenvalues from the inverse iteration method can be easily checked with the calculations of the 'Classic' Jacobi eigenvalue method.\n");

// We call the modified jacobi_diagonalisation
fprintf(stdout, "For test, we call the classic jacobi_diagonalisation with optimised termination of the largest elements in each row\n");
int sweeps_classic;
sweeps_classic = jacobi_classic(A_classic,e_classic,V_classic);
fprintf(stdout, "The Classic Jacobi method ended in %d sweeps.\n",sweeps_classic);
fprintf(stdout, "The Classic Jacobi-method returned\n");
matrix_print("Eigenvectors V=",V_classic);
vector_print("Eigenvalues e=",e_classic);

// We check the method has worked
fprintf(stdout, "The check the method by calculation of V^T * A * V, which should be a matrix D with all eigenvalues on the diagonal.\n");
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,V_classic,A,0.0,Product); // Product = 1 * V^T * A + 0 * Product
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Product,V_classic,0.0,D); // D = 1 * Product * V + 0 * D
matrix_print("V^T * A * V = D = ",D);

fprintf(stdout, "We clearly see that we are able to reproduce the D-matrix, with the eigenvalues and must hereby have revieved the correct eigenvectors and eigenvalues from the classic jacobi algorithm.\n\n");


fprintf(stdout, "\n\n ------------- Test 1. Power iteration method -----------------\n" );

fprintf(stdout, "At first as a steppingstone to the inverse iteration method, we estimate the largest eigenvalue and eigenvector using the related but simpler power iteration method, which has been implemented.\n");
int iter_pow;
iter_pow = power_iteration(A, e_power, V_power);
fprintf(stdout, "The Power iteration method ended in %d iterations.\n",iter_pow);
fprintf(stdout, "The  Power iteration method returned\n");
vector_print("Largest Eigenvector V=",V_power);
vector_print("Largest Eigenvalue e=",e_power);
fprintf(stdout, "We see that the simple power iteration method can reproduce the largest eigenvalue and eigenvector by comparision with the Classic Jacobi method.\n");

fprintf(stdout, "\n\n ------------- Test 2. Inverse power iteration method -----------------\n" );

fprintf(stdout, "As a further steppingstone we consder the inverse power iteration, that works smilarly to the power interation method, but can find the lowest eigenvalue. This has been implemented and will now be tested. We find the inverse using teh Matrix inverse algorithm by Golub-Kahan-Lanczos bidiagonalization as in the linear equation-C exercise.  \n");
int iter_inv_pow;
iter_inv_pow = inverse_power_iteration(A, e_inv_power, V_inv_power);
fprintf(stdout, "The inverse power iteration method ended in %d iterations.\n",iter_inv_pow);
fprintf(stdout, "The inverse power iteration method returned\n");
vector_print("Smallest Eigenvector V=",V_inv_power);
vector_print("Smallets Eigenvalue e=",e_inv_power);
fprintf(stdout, "We see that the inverse power iteration method can reproduce the smallest eigenvalue and eigenvector by comparision with the Classic Jacobi method. If we want to find an arbitrary eigenvalue the simplest way is to implement a shifted inverse iteration, to check that the procedure works, before we reach the more clever inverse iteration method.\n");

fprintf(stdout, "\n\n ------------- Test 3. Shifted Inverse iteration method -----------------\n" );

int iter_shift_inv = 0;
double s = -5.;
iter_shift_inv = shifted_inverse_iteration(A, e_shift_inv, V_shift_inv,s);
fprintf(stdout, "The shifted inverse iteration around s = %g ended in %i iterations.\n",s,iter_shift_inv);
fprintf(stdout, "The shifted inverse iteration method around s = %g returned \n",s);
vector_print("Closest Eigenvector V=",V_shift_inv);
vector_print("Closest Eigenvalue e=",e_shift_inv);
fprintf(stdout, "We see that the shifted inverse iteration method can reproduce the nearest eigenvalue and eigenvector by comparision with the results of the Classic Jacobi method. This method could be made more efficient by updating the shift, to the approximated eigenvalue for each step, but since it finishes in only %i iterations this is only considered when going to the more clever inverse iteration method.\n\n\n",iter_shift_inv);

fprintf(stdout, "\n\n ------------- Final implementation: Inverse Iteration Method -----------------\n" );
fprintf(stdout, "Finally the inverse iteration method has been made. This do not rely on the Golub-Kahan-Lanczos bidiagonalization for finding the inverse, but instead the more simple QR-decomposition to solve a linear system using backsubstitution instead of finding the inverse. This procedure will if we do not update the estimated eigenvalue cost O(n²) operations per iteration, but with a reestimation of the shift for every iteration the cost will be O(n³). Therefore the implementation has been made, so the user will be able to ask for the desired iterations per update of the estimated eigenvalue. \n");

fprintf(stdout, "To directly see the efficiency of the inverse iteration method, we do not update the estimate of the eigenvalue in the following test, for direct comparison.\n");

int iter_inv = 0;
s = -5;
int update = 3;
iter_inv = inverse_iteration(A, e_inv, V_inv,s,update);
fprintf(stdout, "The inverse iteration around s = %g ended in %i iterations, with an update every %i'th iteration'.\n",s,iter_inv,update);
fprintf(stdout, "The inverse iteration method around s = %g returned \n",s);
vector_print("Closest Eigenvector V=",V_inv);
vector_print("Closest Eigenvalue e=",e_inv);

fprintf(stdout, "Two different implementations of the inverse iteration algorithm has hereby been made to find to nearest eigenvalue and eigenvector. They are both demonstrated, and are relying on finding the inverse through olub-Kahan-Lanczos bidiagonalization, and solving a linear system using QR-factorisation with backsubstitution respectively.\n");
fprintf(stdout, "A in-depth examination of the procedure are made in the corresponding report.pdf where additional information can be found.\n");
fprintf(stdout, "The quick demonstration is hereby done.\n");

// We free the parameters
gsl_matrix_free(A);
gsl_matrix_free(A_classic);
gsl_matrix_free(D);
gsl_matrix_free(Product);
gsl_matrix_free(V_classic);
gsl_vector_free(V_power);
gsl_vector_free(V_inv_power);
gsl_vector_free(V_shift_inv);
gsl_vector_free(V_inv);
gsl_vector_free(e_classic);
gsl_vector_free(e_power);
gsl_vector_free(e_inv_power);
gsl_vector_free(e_shift_inv);
gsl_vector_free(e_inv);


return 0;
}
