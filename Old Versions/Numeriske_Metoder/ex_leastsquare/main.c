#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include"least_square_functions.h"
#include"lineqfunctions.h"
#include"eigenfunctions.h"




int main(int argc, char** argv){

// Exercise A
// data load
fprintf(stdout, "\nExercise A has started\n\n" );
fprintf(stdout, "An algoritm for making an ordinary least-squares-fit with uncertainties using QR decomposition have been implemented, and will now be tested on the given datapoints.\n" );
int n=atoi(argv[1]);
double xdata[n],ydata[n],dydata[n];
 for(int i=0;i<n;i++){
  scanf("%lg %lg %lg",xdata+i,ydata+i,dydata+i);}

// Load data into vectors
gsl_vector* vec_xdata = gsl_vector_alloc(n);
gsl_vector* vec_ydata = gsl_vector_alloc(n);
gsl_vector* vec_dydata = gsl_vector_alloc(n);

for(int i=0;i<n;i++){
  gsl_vector_set(vec_xdata,i,xdata[i]);
  gsl_vector_set(vec_ydata,i,ydata[i]);
  gsl_vector_set(vec_dydata,i,dydata[i]);
  }

fprintf(stdout, "The data has been loaded.\n" );

fprintf(stdout, "We start the fitting procedure using QR-decomposition.\n" );

int freeparams = 3; // The number of parameters for the function to fit is set
gsl_matrix* S = gsl_matrix_alloc(freeparams,freeparams);
gsl_vector* c = gsl_vector_alloc(freeparams);
lsfitA(funs,freeparams,vec_xdata,vec_ydata,vec_dydata,c,S);

fprintf(stdout, "The fitting has ended.\n" );

vector_print("We have using QR-decomposition found solution vector c = ", c);
matrix_print("We have using QR-decomposition found covariance matrix S = ", S);

// We know write the data into file data1.txt, this is done from the first to the last element in 1000 steps
FILE* file = fopen("data1.txt", "w");
double x_min = gsl_vector_get(vec_xdata, 0), x_max = gsl_vector_get(vec_xdata, (vec_xdata->size)-1)+1e-5;
double stepx = (x_max - x_min)/1000;

// Here we find the calculated fit values:
for(double i = x_min; i < x_max; i += stepx) {
    double f = 0.00;
      for(int j = 0; j < c->size; ++j) {
        double c_j = gsl_vector_get(c, j);
        f += c_j*funs(j,i);
    }
    fprintf(file, "%g %g\n", i, f);
}

fprintf(stdout, "The fitted line has been calculated, and the fitted line is depicted with the experimental data in the first figure.\n" );
fprintf(stdout, "As seen in the plot the fit is accurate. This demonstration concludes exercise A. \n" );


// Exercise B

fprintf(stdout, "\n \n Exercise B has started.\n\n" );
fprintf(stdout, "We reuse the function for exercise A, since here the requested covariance matrix is already found.\n" );
// We reuse the data from exercise A, since that function is already made complete.
// However this time we find the correlations from the covariance matrix. By adding or subtracting this contribution we find the errors of our fit.
// This is done now

fprintf(stdout, "We calculate the new fits, based on the covariance matrix.\n" );
matrix_print("The covariance matrix are: ",S);


// We know write the data into file data2.txt, this is done from the first to the last element in 1000 steps
FILE* file2 = fopen("data2.txt", "w");

// Here we find the calculated fit values along with an addition and substraction of the error:
for(double i = x_min; i < x_max; i += stepx) {
    double f = 0.00;
    double df = 0.00;
      for(int j = 0; j < c->size; ++j) {
        double c_j = gsl_vector_get(c, j);
 // Here we insert an error from  the covariance matrix elements as the sum of all covariances:
          for(int k = 0; k < c->size; ++k){
            df += gsl_matrix_get(S,j,k) * funs(j,i)* funs(k,i);
          }
        f += c_j*funs(j,i);
    }
// We now print the values and remember the square root on the errors
    fprintf(file2, "%g %g %g\n", i, f, sqrt(df));
}

fprintf(stdout, "The calculations have ended and the fitted functions is plotted with the experimental data on the second figure.\n" );
fprintf(stdout, "As seen from the plot, the fit and the corresponding error to the fit is accurate to describe the experimental datapoints.\n" );
fprintf(stdout, "This examination concludes exercise B.\n");

// Exercise C.
fprintf(stdout, "\n \n Exercise C has started. \n\n" );
fprintf(stdout, "A function which makes singular value decompostion, to solve the linear and ordinary squares problem have been implemented. \n" );
// First we allocate all of the needed variables for the SVD, and the tests hereoff:

gsl_matrix *A    = gsl_matrix_alloc(n,freeparams);
gsl_matrix *V    = gsl_matrix_alloc(freeparams,freeparams);
gsl_matrix *D    = gsl_matrix_alloc(freeparams,freeparams);
gsl_matrix *U    = gsl_matrix_alloc(n,freeparams);
gsl_matrix *S2    = gsl_matrix_alloc(freeparams,freeparams);
gsl_matrix *prod    = gsl_matrix_alloc(freeparams,freeparams);
gsl_matrix *Atest    = gsl_matrix_alloc(n,freeparams);


// We insert the wanted valies in A
for(int i=0;i<n;i++){
	double xi  = gsl_vector_get(vec_xdata ,i);
	double dyi = gsl_vector_get(vec_dydata,i);
	for(int k=0;k<freeparams;k++){
    gsl_matrix_set(A,i,k,funs(k,xi)/dyi);}
	}

fprintf(stdout, "We start by finding the singular value decomposition of A\n" );
// The singular value decompostion is found using the functions
SVD(A, V, D, U, S2);
// We print the results
fprintf(stdout, "SVD has ended to find the following matrices:.\n" );
matrix_print("Initial Matrix A=",A);
matrix_print("Eigenmatrix for A^T A denoted V=",V);
matrix_print("Eigenvalues for A^T A denoted on the diagonal D=",D);
matrix_print("Othogonal Basis U=",U);
matrix_print("Singuular Value matrix S=",S2);

fprintf(stdout, "We test this by calculating A from decomposition U S V^T.\n" );

// We test the decomposition
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,S2,V,0,prod);  // prod = S2 V^T
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,U,prod,0,Atest);  // Atest = U  prod

matrix_print("We find succesively U * S * V^T = ", Atest);

// Now we want to use this, for solving a linear equation problem
fprintf(stdout, "We test that the decomposition can solve a linear least squares problem, by a random matrix problem as done previously in the exercise on linear equations.\n" );
//First we generate a random 3x3 matrix and vector.
gsl_matrix *MatrixA=gsl_matrix_alloc(3,3);
gsl_matrix *MatrixAcopy=gsl_matrix_alloc(3,3);
gsl_matrix *MatrixV    = gsl_matrix_alloc(3,3);
gsl_matrix *MatrixD    = gsl_matrix_alloc(3,3);
gsl_matrix *MatrixU    = gsl_matrix_alloc(3,3);
gsl_matrix *MatrixS2    = gsl_matrix_alloc(3,3);
gsl_vector *VectorB=gsl_vector_alloc(3);
gsl_vector *VectorX=gsl_vector_alloc(3);
gsl_vector *VectorBtest=gsl_vector_alloc(3);

//Now we insert random numbers from 0 to 10 in the matrix and vector:
for (int i = 0; i < MatrixA->size1; ++i) {
    for (int j = 0; j < MatrixA->size2; ++j) {
        gsl_matrix_set(MatrixA, i, j, ((double) rand())/((double)RAND_MAX)*10);
    }
}
for (int i = 0; i < VectorB->size; ++i) {
        gsl_vector_set(VectorB, i, ((double) rand())/((double)RAND_MAX)*10);
}
gsl_matrix_memcpy(MatrixAcopy,MatrixA);

// Now we print the Beginning matrix and Vector:
fprintf(stdout, "The random (from 0 to 10) initial matrix is:\n");
matrix_print("A=",MatrixA);
fprintf(stdout, "The random (from 0 to 10) vector is:\n");
vector_print("b=",VectorB);
// We call the function to do both the SV-decomp and the solving of the linear system
fprintf(stdout, "We start the singular value decomposition, and solutions of the linear system:\n");
linsysSVD(MatrixA, MatrixV, MatrixD, MatrixU, MatrixS2, VectorB, VectorX);
// We print the solution
fprintf(stdout, "The singular value decomposition has ended.\n");
vector_print("We find the solution for the linear system as X= ", VectorX);
// We test the solution will give the correct result
gsl_blas_dgemv(CblasNoTrans,1,MatrixAcopy,VectorX,0,VectorBtest);  // Btest = A * X
vector_print("We find the solution correct since B = A*X = ", VectorBtest);
fprintf(stdout, "We can hereby find solutions for a typical linear equation problem using the least squares model from the Singular Value decomposition.\n");


fprintf(stdout, "We now for the ultimate test implement a function that can give the solution for an ordinary least-squares fit, by this singular-value decomposition.\n");
// We allocate the needed new vectors

gsl_matrix* SSV = gsl_matrix_alloc(freeparams,freeparams);
gsl_vector* cSV = gsl_vector_alloc(freeparams);


fprintf(stdout, "Ordinary Least square fit for singular-value decomposition has started.\n");

lsfitC(funs,freeparams,vec_xdata,vec_ydata,vec_dydata,cSV,SSV);
fprintf(stdout, "Ordinary Least square fit for singular-value decomposition has ended.\n");

vector_print("We have using SV-decomposition found solution vector c = ", cSV);
matrix_print("We have using SV-decomposition found covariance matrix S = ", SSV);

// We know write the data into file data3.txt, this is done from the first to the last element in 1000 steps
FILE* file3 = fopen("data3.txt", "w");

// Here we find the calculated fit values along with an addition and substraction of the error:
for(double i = x_min; i < x_max; i += stepx) {
    double f = 0.00;
    double df = 0.00;
      for(int j = 0; j < c->size; ++j) {
        double c_j = gsl_vector_get(cSV, j);
 // Here we insert an error from  the covariance matrix elements as the sum of all covariances:
          for(int k = 0; k < cSV->size; ++k){
            df += gsl_matrix_get(SSV,j,k) * funs(j,i)* funs(k,i);
          }
        f += c_j*funs(j,i);
    }
// We now print the values and remember the square root on the errors
    fprintf(file3, "%g %g %g\n", i, f, sqrt(df));
}

fprintf(stdout, "The ordinary least squares fit has been calculated, and the fitted line is depicted with the experimental data in the third figure.\n" );
fprintf(stdout, "As seen in the plot the fit is accurate, just as in previous exercises. \n  This concludes exercise C and the exermination of least square fits. \n" );

// We free all of the parameters

gsl_vector_free(vec_xdata);
gsl_vector_free(vec_ydata);
gsl_vector_free(vec_dydata);
gsl_vector_free(c);
gsl_matrix_free(S);
gsl_vector_free(cSV);
gsl_matrix_free(SSV);
gsl_matrix_free(A);
gsl_matrix_free(V);
gsl_matrix_free(D);
gsl_matrix_free(U);
gsl_matrix_free(S2);
gsl_matrix_free(prod);
gsl_matrix_free(Atest);
gsl_matrix_free(MatrixA);
gsl_matrix_free(MatrixAcopy);
gsl_matrix_free(MatrixV);
gsl_matrix_free(MatrixD);
gsl_matrix_free(MatrixU);
gsl_matrix_free(MatrixS2);
gsl_vector_free(VectorB);
gsl_vector_free(VectorX);
gsl_vector_free(VectorBtest);
 return 0;
}
