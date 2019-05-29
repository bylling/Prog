#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include"least_square_functions.h"
#include"lineqfunctions.h"
#include"eigenfunctions.h"


double funs(int i, double x){ // The functions are stated as in the exercise description
   switch(i){
   case 0: return log(x); break;
   case 1: return 1.0;   break;
   case 2: return x;     break;
   default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
   }
}

void lsfitA(double funs(int i, double x),int freeparams,gsl_vector* vec_xdata, gsl_vector* vec_ydata, gsl_vector* vec_dydata, gsl_vector* c, gsl_matrix* S ){
  // The least_fit function for exercise A. This is done corresponding to the  python script from the lecture notes.
int n = vec_xdata->size; // The number of datapoints are found.
int m = freeparams; // The number of free parameters in the functions are collected and denoted as required

// We allocate the needed matrices and vectors in the correct sizes corresponding to the values.
gsl_matrix *A    = gsl_matrix_alloc(n,m);
gsl_vector *b    = gsl_vector_alloc(n);
gsl_matrix *R    = gsl_matrix_alloc(m,m);
gsl_matrix *I    = gsl_matrix_alloc(m,m);
gsl_matrix* R_inverse = gsl_matrix_alloc(m,m);
// For the fit, we must fill our x-data into function f(x) to insert into matrix A, and solution into vector B

for(int i=0;i<n;i++){
	double xi  = gsl_vector_get(vec_xdata ,i);
	double yi  = gsl_vector_get(vec_ydata ,i);
	double dyi = gsl_vector_get(vec_dydata,i);
  gsl_vector_set(b,i,yi/dyi);
	for(int k=0;k<m;k++){
    gsl_matrix_set(A,i,k,funs(k,xi)/dyi);}
	}

// We want to solve this using QR decomposition, to solve R * c = Q^T b
// Therefore we decompose A into its QR decomposition by the function from the "linear equation exercise":
qr_gs_decomp(A,R);
// Now A = Q
// By backsubstitution we solve R* C = Q^T * B, this is done by the function from the "linear equation exercise"
qr_gs_solve(A,R,b,c);



// With this we can find the covariance matrix, and we first find R^-1 with a function from the "linear quation exercise"
// For this we must solve a linear equation with the identity as input, we therefore set gsl_matrix_set_identity
gsl_matrix_set_identity(I);
// We find the iverse of r
qr_gs_inverse(I,R,R_inverse);
// Now we find the covariance matrix by  R^-1 * (R^-1)^T
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,R_inverse,R_inverse,0,S); // R⁻1 * R-1

// We free the used matrices, that are not returned
gsl_matrix_free(A);
gsl_vector_free(b);
gsl_matrix_free(R);
gsl_matrix_free(I);
gsl_matrix_free(R_inverse);
}


void SVD(gsl_matrix* A, gsl_matrix* V, gsl_matrix* D, gsl_matrix* U,  gsl_matrix* S) {
// This function finds the singular value decompostion and returns all of the found elements to their matrices.

// First we allocate the used parameters:
gsl_matrix *Aprod = gsl_matrix_alloc(A->size2,A->size2);
gsl_vector *e = gsl_vector_alloc(A->size2);
gsl_matrix *Dminhalf = gsl_matrix_alloc(A->size2,A->size2);
gsl_matrix *VDprod = gsl_matrix_alloc(A->size2,A->size2);


// To find A^T A= V D V^T we first calculate the left side
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,A,A,0,Aprod);  // Aprod = A^T A


//matrix_print("Aprod= ",Aprod);

// Hereafter we solve the Aprod = VDV^T by identities of the eigenvalues on the diagonal in D and the eigenvectors of V
// This is done as in the previous eigenvector exercise, be reusing the simple jacobi method
jacobi(Aprod,e,V);

// We set up the D matrix with the found eigenvalues:
gsl_matrix_set_identity(D);
for (int i = 0; i < A->size2; i++) {
  gsl_matrix_set(D,i,i,gsl_vector_get(e,i));
}

// Now we find the rest of the relevant matrices:
// We start by the singular value matrix S=Sqrt(D), since D only consists of diagonal elements, this is easily done by:
gsl_matrix_set_identity(S);
for (int i = 0; i < A->size2; i++) {
  gsl_matrix_set(S,i,i,sqrt(gsl_vector_get(e,i)));
}


// We can futhermore find the covariance as  Σ= (A^TA)^−1= (VS^2V^T)^−1=VS^−2V^T


// Hereafter we find U as U= A * V * D^-1/2, through a simple calculation, and remembering once again that D is diagonal so:
for (int i = 0; i < A->size2; i++) {
  gsl_matrix_set(Dminhalf,i,i,1.0/sqrt(gsl_vector_get(e,i)));
}
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,Dminhalf,0,VDprod);  // VDPROD = V D^-1/2
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,A,VDprod,0,U);  // U = A VDPROD

// We free the intermediate values.
gsl_matrix_free(Aprod);
gsl_matrix_free(Dminhalf);
gsl_matrix_free(VDprod);
gsl_vector_free(e);
}



void linsysSVD(gsl_matrix* MatrixA,  gsl_matrix* MatrixV,  gsl_matrix* MatrixD,  gsl_matrix* MatrixU,  gsl_matrix* MatrixS2, gsl_vector* VectorB, gsl_vector* VectorX){
// The function to do both the SV-decomp and the solving of the linear system
// First we allocate the needed parameters
gsl_vector *y = gsl_vector_alloc(MatrixA->size2);
gsl_vector *prod = gsl_vector_alloc(MatrixA->size2);
gsl_matrix *Sinverse = gsl_matrix_alloc(MatrixA->size2,MatrixA->size2);


// This function finds the solution to the linear system of equations set by A X = B. This is dony by firstly a single value decomposition of A, afterwards a solving procedure is applied for finding X.
SVD(MatrixA, MatrixV, MatrixD, MatrixU, MatrixS2);
// With the singular value decomposition we now find the solution of the diagonal system S* Y = U ^T *B
// Since S only have the diagonal we can multiply with S^-1 to find Y = S^-1 * U ^T *B
// First we find S^-1
for (int i = 0; i < MatrixA->size2; i++) {
  gsl_matrix_set(Sinverse,i,i,1.00/(gsl_matrix_get(MatrixS2,i,i)));
}
// Then we calculate Y
gsl_blas_dgemv(CblasTrans,1,MatrixU,VectorB,0,prod);  // prod = U^T B
gsl_blas_dgemv(CblasNoTrans,1,Sinverse,prod,0,y);  // Y = S^-1 * prod
// Now we act with V on Y to find the solution X
gsl_blas_dgemv(CblasNoTrans,1,MatrixV,y,0,VectorX);  // Y = S^-1 * prod

gsl_vector_free(y);
gsl_vector_free(prod);
gsl_matrix_free(Sinverse);
}


void lsfitC(double funs(int i, double x),int freeparams,gsl_vector* vec_xdata, gsl_vector* vec_ydata, gsl_vector* vec_dydata, gsl_vector* c, gsl_matrix* S ){ // The least_fit function for exercise C.
  // This function uses singular value decompostion to fit function F on the data from the vectors, to return the best fit, and the propper errors.

  int n = vec_xdata->size; // The number of datapoints are found.
  int m = freeparams; // The number of free parameters in the functions are collected and denoted as required

  // We allocate the needed matrices and vectors in the correct sizes corresponding to the values.
  gsl_matrix *MatrixA    = gsl_matrix_alloc(n,m);
  gsl_matrix *MatrixV    = gsl_matrix_alloc(m,m);
  gsl_matrix *MatrixD    = gsl_matrix_alloc(m,m);
  gsl_matrix *MatrixU    = gsl_matrix_alloc(n,m);
  gsl_matrix *MatrixS    = gsl_matrix_alloc(m,m);
  gsl_matrix *Sinverse    = gsl_matrix_alloc(m,m);
  gsl_matrix *Sinverse2    = gsl_matrix_alloc(m,m);
  gsl_matrix *prodS    = gsl_matrix_alloc(m,m);
  gsl_vector *y    = gsl_vector_alloc(m);
  gsl_vector *VectorB = gsl_vector_alloc(n);
  gsl_vector *prod = gsl_vector_alloc(m);

  // For the fit, we must fill our x-data into function f(x) to insert into matrix A, and solution into vector B

  for(int i=0;i<n;i++){
  	double xi  = gsl_vector_get(vec_xdata ,i);
  	double yi  = gsl_vector_get(vec_ydata ,i);
  	double dyi = gsl_vector_get(vec_dydata,i);
    gsl_vector_set(VectorB,i,yi/dyi);
  	for(int k=0;k<m;k++){
      gsl_matrix_set(MatrixA,i,k,funs(k,xi)/dyi);}
  }
// We call the singular value decompostion, to do the least square fit.
SVD(MatrixA, MatrixV, MatrixD, MatrixU, MatrixS);

// With the singular value decomposition we now find the solution of the diagonal system S* Y = U ^T *B
// Since S only have the diagonal we can multiply with S^-1 to find Y = S^-1 * U ^T *B
// First we find S^-1
for (int i = 0; i < MatrixA->size2; i++) {
  gsl_matrix_set(Sinverse,i,i,1.00/(gsl_matrix_get(MatrixS,i,i)));
}

// Then we calculate Y
gsl_blas_dgemv(CblasTrans,1,MatrixU,VectorB,0,prod);  // prod = U^T B
gsl_blas_dgemv(CblasNoTrans,1,Sinverse,prod,0,y);  // Y = S^-1 * prod
// Now we act with V on Y to find the solution X
gsl_blas_dgemv(CblasNoTrans,1,MatrixV,y,0,c);  // Y = S^-1 * prod

// We can futhermore find the covariance as  Σ= (A^TA)^−1= (VS^2V^T)^−1=VS^−2V^T
// First we find S^-2
for (int i = 0; i < MatrixA->size2; i++) {
  gsl_matrix_set(Sinverse2,i,i,1.00/(gsl_matrix_get(MatrixS,i,i) * gsl_matrix_get(MatrixS,i,i)));
}
// Now we calculate: Σ=VS^−2V^T
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,Sinverse2,MatrixV,0,prodS);  // prods = S^−2V^T
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,MatrixV,prodS,0,S);  // Σ=V prods

// This conlcudes the function
  // We free the used matricies and vectors
  gsl_matrix_free(MatrixA);
  gsl_matrix_free(MatrixV);
  gsl_matrix_free(MatrixD);
  gsl_matrix_free(MatrixU);
  gsl_matrix_free(Sinverse);
  gsl_matrix_free(Sinverse2);
  gsl_matrix_free(prodS);
  gsl_matrix_free(MatrixS);
  gsl_vector_free(VectorB);
  gsl_vector_free(prod);
  gsl_vector_free(y);
}
