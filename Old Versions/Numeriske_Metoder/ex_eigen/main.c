#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include<math.h>
#include"eigenfunctions.h"


int main(){
fprintf(stdout, "\nExercise A is started.\n\n");
fprintf(stdout, "An algorithm using Jacobi diagonalization by cyclic sweeps have been implemented and a test with a n=4 matrix is shown:\n");
// We initialise the wanted parameters:
int systemsize = 4;
gsl_matrix *A=gsl_matrix_alloc(systemsize,systemsize);
gsl_matrix *Acopy=gsl_matrix_alloc(systemsize,systemsize);
gsl_matrix *Acopy2=gsl_matrix_alloc(systemsize,systemsize);
gsl_matrix *AcopyC=gsl_matrix_alloc(systemsize,systemsize);
gsl_matrix *V=gsl_matrix_alloc(systemsize,systemsize);
gsl_matrix *V2=gsl_matrix_alloc(systemsize,systemsize);
gsl_matrix *V3=gsl_matrix_alloc(systemsize,systemsize);
gsl_matrix *VC=gsl_matrix_alloc(systemsize,systemsize);
gsl_matrix *prod=gsl_matrix_alloc(systemsize,systemsize);
gsl_matrix *D=gsl_matrix_alloc(systemsize,systemsize);
gsl_vector *e=gsl_vector_alloc(systemsize);
gsl_vector *e2=gsl_vector_alloc(systemsize);
gsl_vector *e3=gsl_vector_alloc(systemsize);
gsl_vector *eC=gsl_vector_alloc(systemsize);
// We fill the matrix with random real numbers
double randomnumber;
for (int i = 0; i < A->size1; ++i) {
    for (int j = i; j < A->size2; ++j) {
      randomnumber = ((double) rand())/((double)RAND_MAX)*10;
        gsl_matrix_set(A, i, j, randomnumber);
        gsl_matrix_set(A, j, i, randomnumber);
    }
}
gsl_matrix_memcpy(Acopy,A);
gsl_matrix_memcpy(Acopy2,A);
gsl_matrix_memcpy(AcopyC,A);
// Now we print the Beginning matrix:
fprintf(stdout, "The random (from 0 to 10) initial symmetric matrix is:\n");
matrix_print("A=",A);
// We do now initialise the cyclic sweeping jacobi-method from the lecture notes:
int sweeps;
sweeps = jacobi(A,e,V);
fprintf(stdout, "The Cyclic Jacobi method ended in %d sweeps\n",sweeps);
// We print the result
fprintf(stdout, "The Jacobi-method returned\n ");
matrix_print("Eigenvectors V=",V);
vector_print("Eigenvalues e=",e);
// We check the method has worked
fprintf(stdout, "The check the method by calculation of V^T * A * V, which should be a matrix D with all eigenvalues on the diagonal.\n");
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,V,Acopy,0.0,prod); // prod = 1 * V^T * A + 0 * prod
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,prod,V,0.0,D); // prod = 1 * V^T * A + 0 * prod
matrix_print("V^T * A * V =",D);
// We free the parameters
fprintf(stdout, "With a succesfull demonstration hereof, we examine the calculation time for different matrix dimentsions. The data from using the cyclic jacobi method are plotted in the corresponding figure.\n");
fprintf(stdout, "The graph clearly depicts the expected cubed dependency. With this examination, the first exercise is done. \n");




// Exercise B
fprintf(stdout, "\n \n Exercise B is started.\n\n");
fprintf(stdout, "An Jacobi diagonalization method for recieving the eigenvalues one by one have been implemented. A test with the same matrix is shown, so the correct values can easily be checked with the calculations of exercise A:\n");
// We call the modified jacobi_diagonalisation
fprintf(stdout, "For test, we call the modified jacobi_diagonalisation and ask it to find first one eigenvalue:\n");
int sweeps2;
sweeps2 = jacobi_mod_ev(Acopy,e2,V2,1);
fprintf(stdout, "The modified Jacobi method ended in %d calculations\n",sweeps2);
// We print the result
fprintf(stdout, "The Jacobi-method returned\n");
matrix_print("Eigenvector V=",V2);
vector_print("Eigenvalue e=",e2);
fprintf(stdout, "The smallest eigenvalue have returned, since this is how the angle in the rotations for the itterations are defined, through the sign of the sinefunction, and the sign of the fraction, from wich the angle is found in the arctan2 funciton.\n");
fprintf(stdout, "We now modify the jacobi diagonalization once again but to return the highest eigenvalue instead, by transposing the rotaion matrix, and finding the corresponding other solution to the angle phi that solves the rotation equations.\n");
fprintf(stdout, "This is done by: .\n");
int sweeps3;
sweeps3 = jacobi_mod_highev(Acopy2,e3,V3,1);
fprintf(stdout, "The modified Jacobi method ended in %d sweeps\n",sweeps3);
// We print the result
fprintf(stdout, "The Jacobi-method returned\n");
matrix_print("Eigenvectors V=",V3);
vector_print("Eigenvalues e=",e3);
fprintf(stdout, "A comparison between the cyclic model and finding the eigenvalues individually is made and plotted in the second figure.\n");

fprintf(stdout, "From the figure, we see that it is slightly less fast to calculate the elements one by one, but if you would be able to start the next processes by getting the eigenvalues one by one, then there is time to safe by initialising the next part of a calculation with the eigenvalues one by one.\n However the first eigenvalue is still the most difficult to calculate, so there is still a significant wait before continuing the process.\n");

fprintf(stdout, "The function works properly and will produce the correct eigenvalues and vectors. With this conclusion, the second exercise is done.\n");



// Exercise C
fprintf(stdout, "\n \n Exercise C is started.\n\n");
fprintf(stdout, " The 'Classic' Jacobi eigenvalue method have been implemented, and a test with the same matrix is shown, so the correct values can easily be checked with the calculations of exercise A:\n");

// We call the modified jacobi_diagonalisation
fprintf(stdout, "For test, we call the classic jacobi_diagonalisation with optimised termination of the largest elements in each row\n");
int sweepsC;
sweepsC = jacobi_versc(AcopyC,eC,VC);
fprintf(stdout, "The modified Jacobi method ended in %d sweeps\n",sweepsC);
fprintf(stdout, "One shall remember that theese found sweeps is thoughout all of the matrix. This means that the cyclic method cannot sweep-wise be compared, since it calculates every element in a sweep, were this indexed method only calculate one element in each collumn.\n");

// We print the result
fprintf(stdout, "The Classical Jacobi-method returned\n");
matrix_print("Eigenvectors V=",VC);
vector_print("Eigenvalues e=",eC);
fprintf(stdout, "With this succesfull demonstration, we can observe that the number of rotations may increase, but since it uses way fewer operations per rotaion it decreases computing time and are more efficient. We can compare the computation time between the cyclic and classical jacobi method for different matrix dimensions, this is shown in the corresponding figure.  The function work as planned, and returns the correct eigenvectors and values. This concludes the last exercise and the examniation of eigenvalues. \n");

// the parameters are freed
gsl_matrix_free(A);
gsl_matrix_free(Acopy);
gsl_matrix_free(Acopy2);
gsl_matrix_free(AcopyC);
gsl_matrix_free(V);
gsl_matrix_free(V2);
gsl_matrix_free(V3);
gsl_matrix_free(VC);
gsl_matrix_free(prod);
gsl_matrix_free(D);
gsl_vector_free(e);
gsl_vector_free(e2);
gsl_vector_free(e3);
gsl_vector_free(eC);
return 0;
}
