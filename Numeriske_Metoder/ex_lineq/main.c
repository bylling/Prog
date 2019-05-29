#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#include<math.h>
#include"lineqfunctions.h"



int main(){
fprintf(stdout, "\nExercise A.1 is started\n\n");

fprintf(stdout, "An algorithm for QR-decomposition have been implemented, as well as a fucntion to sovle a linear system by Gram-Schmidt in-place-backsubstitution. Theese algorithms will now be tested.\n\n");

//First we generate a random 4x3 matrix so we can test a tall matrix.
gsl_matrix *A=gsl_matrix_alloc(4,3);
gsl_matrix *R=gsl_matrix_alloc(3,3);
//Now we insert random numbers from 0 to 10 in the matrix by procedure from the MonteCarlo exercise:
for (int i = 0; i < A->size1; ++i) {
    for (int j = 0; j < A->size2; ++j) {
        gsl_matrix_set(A, i, j, ((double) rand())/((double)RAND_MAX)*10);
    }
}
// Now we print the Beginning matrix:
fprintf(stdout, "The random (from 0 to 10) initial matrix is:\n");
matrix_print("A=",A);

// We start the QR decomposition. using Graham-schmidt method by:.
qr_gs_decomp(A,R);
// We print the results
fprintf(stdout, "QR decomposition gives us:\n");
matrix_print("Othonormal basis Q=",A);
matrix_print("Triangular Coefficients R =",R);
// We check Q^T * Q = I
fprintf(stdout, "We check that Q tansposed multiplied with Q equals I, because of orthonormality\n");
gsl_matrix *ident=gsl_matrix_alloc(3,3);
gsl_blas_dgemm(CblasConjTrans,CblasNoTrans,1.0,A,A,0.0,ident);
matrix_print("Q^T * Q=",ident);
// We check Q*R = A
fprintf(stdout, "We check that we can reproduce A with the R-coefficients\n");
gsl_matrix *a_new=gsl_matrix_alloc(4,3);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,R,0.0,a_new);
matrix_print("Q*R=A=",a_new);

// we free the parameters
gsl_matrix_free(A);
gsl_matrix_free(R);
gsl_matrix_free(ident);
gsl_matrix_free(a_new);

fprintf(stdout, "\n\n\nExercise A.2 is started\n\n");

//First we generate a random 3x3 matrix and vector.
gsl_matrix *A2=gsl_matrix_alloc(3,3);
gsl_matrix *A2copy=gsl_matrix_alloc(3,3);
gsl_matrix *R2=gsl_matrix_alloc(3,3);
gsl_vector *b=gsl_vector_alloc(3);
gsl_vector *x=gsl_vector_alloc(3);
gsl_vector *b_test=gsl_vector_alloc(3);

//Now we insert random numbers from 0 to 10 in the matrix and vector by procedure from the MonteCarlo exercise:
for (int i = 0; i < A2->size1; ++i) {
    for (int j = 0; j < A2->size2; ++j) {
        gsl_matrix_set(A2, i, j, ((double) rand())/((double)RAND_MAX)*10);
    }
}
for (int i = 0; i < b->size; ++i) {
        gsl_vector_set(b, i, ((double) rand())/((double)RAND_MAX)*10);
}
gsl_matrix_memcpy(A2copy,A2);

// Now we print the Beginning matrix and Vector:
fprintf(stdout, "The random (from 0 to 10) initial matrix is:\n");
matrix_print("A=",A2);
fprintf(stdout, "The random (from 0 to 10) vector is:\n");
vector_print("b=",b);

// We start the QR decomposition. using Graham-schmidt method by:.
qr_gs_decomp(A2,R2);
// We print the results
fprintf(stdout, "QR decomposition gives us:\n");
matrix_print("Othonormal basis Q=",A2);
matrix_print("Triangular Coefficients R =",R2);

// We find the solution by backsubstitution
qr_gs_solve(A2,R2,b,x);
vector_print("The solution from backsubstitution is x=",x);

// We check that this is a solution by Ax=b:
fprintf(stdout, "We check for correct solution by A*X = B:\n");
gsl_blas_dgemv(CblasNoTrans,1.0,A2copy,x,0.0,b_test); // b = 1*A^T*x + 0*b
vector_print("The check A^T * x = b =",b_test);

// We free the parameters
gsl_matrix_free(A2);
gsl_matrix_free(A2copy);
gsl_matrix_free(R2);
gsl_vector_free(b);
gsl_vector_free(x);
gsl_vector_free(b_test);

fprintf(stdout, "With this succesfull demonstration of the QR-decomposition algorithm and corresponding Gram-Schmidt in-place backsubsitution algorithm for solving linear systems, the first exercise is hereby done. \n" );

// Exercise B
fprintf(stdout, "\n\n\nExercise B is started\n\n");
fprintf(stdout, "An algorithm using QR-decomposition and Gram-Schmidt for calculating a Matric inverse have been implemented, and will now be tested.\n\n");

// We get the wanted matrices
gsl_matrix *A3=gsl_matrix_alloc(3,3);
gsl_matrix *A3copy=gsl_matrix_alloc(3,3);
gsl_matrix *R3=gsl_matrix_alloc(3,3);
gsl_matrix *b2=gsl_matrix_alloc(3,3);

// We insert random numbers as before
for (int i = 0; i < A3->size1; ++i) {
    for (int j = 0; j < A3->size2; ++j) {
        gsl_matrix_set(A3, i, j, ((double) rand())/((double)RAND_MAX)*10);
    }
}
gsl_matrix_memcpy(A3copy,A3);

// Now we print the Beginning matrix:
fprintf(stdout, "The random (from 0 to 10) initial matrix is:\n");
matrix_print("A=",A3);
// We start the QR decomposition. using Graham-schmidt method by:.
qr_gs_decomp(A3,R3);
// We print the results
fprintf(stdout, "QR decomposition gives us:\n");
matrix_print("Othonormal basis Q=",A3);
matrix_print("Triangular Coefficients R =",R3);

// We do the inverse operation:
qr_gs_inverse(A3,R3,b2);

// We print the solution of the inverse of a:
matrix_print("A inverse is found to be B=",b2);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,b2,A3copy,0.0,R3);
matrix_print("We check this by calculation of A*B=I=",R3);

// We free the parameters
gsl_matrix_free(A3);
gsl_matrix_free(A3copy);
gsl_matrix_free(R3);
gsl_matrix_free(b2);

fprintf(stdout, "With this succesfull demonstration of the algorithm for finding matrix inverse by QR-decomposition and Gram-Schmidt, the second exercise is hereby done. \n" );

// Exercise C
fprintf(stdout, "\n\n\nExercise C is started\n\n");
fprintf(stdout, "An algorithm which uses Golub-Kahan-Lanczos bidiagonalization to be able to solve linear systems, find determinants and inverses are implemented, and will now be tested.\n\n");


gsl_matrix *AC=gsl_matrix_alloc(4,4);
gsl_matrix *prod=gsl_matrix_alloc(4,4);
gsl_matrix *ACcopy=gsl_matrix_alloc(4,4);
gsl_matrix *UC=gsl_matrix_alloc(4,4);
gsl_matrix *VC=gsl_matrix_alloc(4,4);
gsl_matrix *BC=gsl_matrix_alloc(4,4);
gsl_vector *s=gsl_vector_alloc(4);
gsl_vector *XC=gsl_vector_alloc(4);
gsl_vector *btestc=gsl_vector_alloc(4);
gsl_matrix *inver=gsl_matrix_alloc(4,4);
// We insert random numbers as before
for (int i = 0; i < AC->size1; ++i) {
    for (int j = 0; j < AC->size2; ++j) {
        gsl_matrix_set(AC, i, j, ((double) rand())/((double)RAND_MAX)*10);
    }
}
for (int i = 0; i < s->size; ++i) {
        gsl_vector_set(s, i, ((double) rand())/((double)RAND_MAX)*10);
}
gsl_matrix_memcpy(ACcopy,AC);
// Now we print the Beginning matrix and vector:
fprintf(stdout, "The random (from 0 to 10) initial matrix is:\n");
matrix_print("A=",AC);
fprintf(stdout, "The random (from 0 to 10) vector s is:\n");
vector_print("s=",s);
// We start the Golub-Kahan-Lanczos Bidiagonalization Procedure.
gkl_biadiag(UC,AC,VC);
// Calculate B
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,UC,AC,0.0,prod); //prod =1 * U^T*A+0*prod
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,prod,VC,0.0,BC); //B =1 * prod*V+0*B

// Printing the results:
fprintf(stdout, "After  Golub-Kahan-Lanczos Bidiagonalization Procedure  we find\n");
matrix_print("U =",UC);
//matrix_print("A =",AC);
matrix_print("V =",VC);
matrix_print("B =",BC);

gkl_solve(UC, AC, VC, s,XC);
fprintf(stdout, "We solve the system of eq. A*X = S using the Golub-Kahan-Lanczos bidiagonalization\n");
vector_print("The solution is x=",XC);
// We check that this is a solution by Ax=b:
fprintf(stdout, "We check for correct solution by A*X = s:\n");
gsl_blas_dgemv(CblasNoTrans,1.0,ACcopy,XC,0.0,btestc); // b = 1*A^T*x + 0*b
vector_print("The check A^T * x = s =",btestc);
// We find the determinant
double determ = 0;
determ = glk_determinant(UC, AC, VC);
fprintf(stdout, "Using Golub-Kahan-Lanczos bidiagonalization we can find the determinant of A to:\n");
fprintf(stdout, "|det(A)| = %g\n",determ);
// We now find the inverse:
glk_inverse(UC,AC,VC, inver);

// We print the solution of the inverse of a:
matrix_print("A inverse is found to be A^-1=",inver);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,inver,ACcopy,0.0,UC);
matrix_print("We check this by calculation of A*B=I=",UC);
fprintf(stdout, "With a succesfull demonstration of finding the determinant, inverse and solution to a linear system all using Golub-Kahan-Lanczos bidiagonalization, this concludes the third exercise, and therby the examination of linear systems.\n");

// We free the parameters
gsl_matrix_free(AC);
gsl_matrix_free(prod);
gsl_matrix_free(ACcopy);
gsl_matrix_free(UC);
gsl_matrix_free(VC);
gsl_matrix_free(BC);
gsl_matrix_free(inver);
gsl_vector_free(XC);
gsl_vector_free(s);
gsl_vector_free(btestc);

return 0;
}
