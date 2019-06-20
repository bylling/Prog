#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"functionfile.h"
#include"function_eigenval.h"
#include"function_lineq.h"

void vector_print(const char* s, gsl_vector* v){ // The function for printing vectors easily
    printf("%s\n",s);
    for(int i=0;i<v->size;i++){
      printf("%8.3g",gsl_vector_get(v,i));
      printf("\n");
    }
    printf("\n");
  }

void matrix_print(const char* s, gsl_matrix* A){ // The function for printing matrices easily
  printf("%s\n",s);
  for(int i=0;i<A->size1;i++){
    for(int j=0;j<A->size2;j++){
      printf("%8.3g",gsl_matrix_get(A,i,j));
      printf("\t");
    }
  printf("\n");
  }
}

int power_iteration(gsl_matrix* A, gsl_vector* e, gsl_vector* V){ // Function for solution of the highest eigenvalue of a matrix using power iteration.
// This function is made to demonstrate the power_iteration algorithm, since it is closely connected to the inverse iteration method.

// We initialize our parameters
int n = A->size1, iter = 0, changed;
gsl_vector* V_new = gsl_vector_alloc(n);
double denominator = 0, numerator = 0;
double Rayleigh=0., NewRayleigh;

// We make a random vector, that will end up containing the eigenvector corresponding to the largest eigenvalue. Therefore we use V.
// We insert random numbers (between 0 and 10) in V
for (int i = 0; i < n; ++i) {
        gsl_vector_set(V, i, ((double) rand())/((double)RAND_MAX)*10);
}

// We want to calculate the new vector x_i+1=Ax_i and solve the  Rayleigh quotient each time, to see when the Rayleigh quotient do not improve any more and hereby contains the eigenvalue.
// This is done with a loop:

do{
  changed=0;
  iter++; // For each iteration, increase iter by one
  gsl_blas_dgemv(CblasNoTrans, 1., A, V,0., V_new);  //x_i+1= A * x_i (This could be changed to gsl_blas_dsymv for symmetric matrices, however we reamain general)
  // Rayleigh quotient is found λ[xi] =x^T_i * A* x_i / x^T_i * x_i = x^T_i+1 * x_i / x^T_i * x_i
  gsl_blas_ddot(V, V, &(denominator)); // We calculate the denominator  x^T_i * x
  gsl_blas_ddot(V_new, V, &(numerator)); // We calculate the numerator x^T_i+1 * x_i
  // We estimate the  Rayleigh quotient
  NewRayleigh = (numerator)/(denominator);
  // If the Rayleigh quotient changes, we break out of the loop
  if (NewRayleigh!=Rayleigh) {
  changed=1;
  }
  // We update the values for the next iteration
  gsl_vector_memcpy(V,V_new); // V -> V_new
  Rayleigh = NewRayleigh;  // Rayleigh -> NewRayleigh
  } while(changed!=0); //We check if the operation had an effect, if not, we stop.

// We normalize the found eigenvector:
gsl_vector_scale(V,1./gsl_blas_dnrm2(V));

// After we have found the new eigenvector V and eigenvalue approximated through the Rayleigh parameter, we insert the eigenvalue to its vector at the largest place and return the value:
for (int i = 0; i < n; ++i) {
        gsl_vector_set(e, i, 0);
        if (i==n-1) {
          gsl_vector_set(e,i,Rayleigh);
        }
}

// We free our parameters
gsl_vector_free(V_new);

// We return the number of iterations
return iter;
}



int inverse_power_iteration(gsl_matrix* A, gsl_vector* e, gsl_vector* V){ // Function for solution of the lowest eigenvalue of a matrix using inverse power iteration.
// This function is made to demonstrate the inverse power_iteration algorithm, since it is closely connected to the inverse iteration method.
// This function is very similar to the power method, and only differs by multiplying with the inverse of A instead of A.

// We initialize our parameters
int n = A->size1, iter = 0, changed;
gsl_vector* V_new = gsl_vector_alloc(n);
gsl_matrix* A_glk = gsl_matrix_alloc(n,n);
gsl_matrix* A_inv = gsl_matrix_alloc(n,n);
gsl_matrix* U_glk = gsl_matrix_alloc(n,n);
gsl_matrix* V_glk = gsl_matrix_alloc(n,n);
double denominator = 0, numerator = 0;
double Rayleigh=0, NewRayleigh;
gsl_matrix_memcpy(A_glk,A);
// We make a random vector, that will end up containing the eigenvector corresponding to the largest eigenvalue. Therefore we use V.
// We insert random numbers (between 0 and 10) in V
for (int i = 0; i < n; ++i) {
        gsl_vector_set(V, i, ((double) rand())/((double)RAND_MAX)*10);
}

// For this method we need to calculate the inverse of A, this is done by the "Golub-Kahan-Lanczos bidiagonalization to find inverses" as in exercise C of "Linear equations".
gkl_biadiag(U_glk, A_glk, V_glk);
glk_inverse(U_glk, A_glk, V_glk, A_inv);
// We want to calculate the new vector x_i+1=Ainv x_i and solve the  Rayleigh quotient each time, to see when the Rayleigh quotient do not improve any more and hereby contains the eigenvalue.
// This is done with a loop:
do{
  changed=0;
  iter++; // For each iteration, increase iter by one
  gsl_blas_dgemv(CblasNoTrans, 1., A_inv, V,0., V_new);  //x_i+1= Ainv * x_i (This could be changed to gsl_blas_dsymv for symmetric matrices, however we reamain general)

  // Rayleigh quotient is found λ[xi] =x^T_i * Ainv * x_i / x^T_i * x_i = x^T_i+1 * x_i / x^T_i * x_i
  gsl_blas_ddot(V, V, &(denominator)); // We calculate the denominator  x^T_i * x
  gsl_blas_ddot(V_new, V, &(numerator)); // We calculate the numerator x^T_i+1 * x_i
  // We estimate the  Rayleigh quotient
  NewRayleigh = (double)(numerator)/(double)(denominator);

  // If the Rayleigh quotient changes, we break out of the loop
  if (NewRayleigh!=Rayleigh) {
  changed=1;
  }
  // We update the values for the next iteration
  gsl_vector_memcpy(V,V_new); // V -> V_new
  Rayleigh = NewRayleigh;  // Rayleigh -> NewRayleigh
  } while(changed!=0); //We check if the operation had an effect, if not, we stop.

// We normalize the found eigenvector:
gsl_vector_scale(V,1./gsl_blas_dnrm2(V));

// After we have found the new eigenvector V and eigenvalue approximated through the inverse of the Rayleigh parameter, we insert the eigenvalue to its vector at the largest place and return the value:
for (int i = 0; i < n; ++i) {
        gsl_vector_set(e, i, 0);
        if (i==0) {
          gsl_vector_set(e,i,1./Rayleigh);
        }
}

// We free our parameters
gsl_vector_free(V_new);
gsl_matrix_free(A_glk);
gsl_matrix_free(A_inv);
gsl_matrix_free(U_glk);
gsl_matrix_free(V_glk);

// We return the number of iterations
return iter;
}



int shifted_inverse_iteration(gsl_matrix* A, gsl_vector* e, gsl_vector* V,double s){ // Function for solution of the wanted eigenvalue near the guess "s" of a matrix using shifted inverse iteration.
// This function is made to demonstrate the shifted inverse iteration algorithm, since it is closely connected to the inverse iteration method.
// This function is very similar to the inverse method, and only differs by a shift in the matrix A by A-s*I.

// We initialize our parameters
int n = A->size1, iter = 0, changed;
gsl_vector* V_new = gsl_vector_alloc(n);
gsl_matrix* A_glk = gsl_matrix_alloc(n,n);
gsl_matrix* Shift = gsl_matrix_alloc(n,n);
gsl_matrix* A_inv = gsl_matrix_alloc(n,n);
gsl_matrix* U_glk = gsl_matrix_alloc(n,n);
gsl_matrix* V_glk = gsl_matrix_alloc(n,n);
double denominator = 0, numerator = 0;
double Rayleigh=0., NewRayleigh;
gsl_matrix_memcpy(A_glk,A);


// We make a random vector, that will end up containing the eigenvector corresponding to the largest eigenvalue. Therefore we use V.
// We insert random numbers (between 0 and 10) in V
for (int i = 0; i < n; ++i) {
        gsl_vector_set(V, i, ((double) rand())/((double)RAND_MAX)*10);
}


// We shift the new Matric A to A - I*s
gsl_matrix_set_identity(Shift);
gsl_matrix_scale(Shift, s); // Shift = s * I
gsl_matrix_sub(A_glk, Shift); //  A = A-Shift
// For this method we need to calculate the inverse of A, this is done by the "Golub-Kahan-Lanczos bidiagonalization to find inverses" as in exercise C of "Linear equations".
gkl_biadiag(U_glk, A_glk, V_glk);
glk_inverse(U_glk, A_glk, V_glk, A_inv);

// We want to calculate the new vector x_i+1=(A-s*I)inv x_i and solve the  Rayleigh quotient each time, to see when the Rayleigh quotient do not improve any more and hereby contains the eigenvalue.
// This is done with a loop:
do{
  changed=0;
  iter++; // For each iteration, increase iter by one
  gsl_blas_dgemv(CblasNoTrans, 1., A_inv, V,0., V_new);  //x_i+1= (Ainv-s*I) * x_i (This could be changed to gsl_blas_dsymv for symmetric matrices, however we reamain general)
  // Rayleigh quotient is found λ[xi] =x^T_i * Ainv * x_i / x^T_i * x_i = x^T_i+1 * x_i / x^T_i * x_i
  gsl_blas_ddot(V, V, &(denominator)); // We calculate the denominator  x^T_i * x
  gsl_blas_ddot(V_new, V, &(numerator)); // We calculate the numerator x^T_i+1 * x_i
  // We estimate the  Rayleigh quotient
  NewRayleigh = (numerator)/(denominator);
  // If the Rayleigh quotient changes, we break out of the loop
  if (NewRayleigh!=Rayleigh) {
  changed=1;
  }
  // We update the values for the next iteration
  gsl_vector_memcpy(V,V_new); // V -> V_new
  Rayleigh = NewRayleigh;  // Rayleigh -> NewRayleigh
  } while(changed!=0); //We check if the operation had an effect, if not, we stop.

// We normalize the found eigenvector:
gsl_vector_scale(V,1./gsl_blas_dnrm2(V));

// After we have found the new eigenvector V and eigenvalue approximated through the inverse Rayleigh parameter plus the shift, we insert the eigenvalue to its vector at the largest place and return the value:
for (int i = 0; i < n; ++i) {
        gsl_vector_set(e, i, 0);
        if (i==0) {
          gsl_vector_set(e,i,1./Rayleigh + s);
        }
}

// We free our parameters
gsl_vector_free(V_new);
gsl_matrix_free(A_glk);
gsl_matrix_free(A_inv);
gsl_matrix_free(Shift);
gsl_matrix_free(U_glk);
gsl_matrix_free(V_glk);

// We return the number of iterations
return iter;
}

int inverse_iteration(gsl_matrix* A, gsl_vector* e, gsl_vector* V,double s,int update){ // Function for solution of the nearest eigenvalue near the guess "s" of a matrix using shifted inverse iteration by QR-decompostion and a user-defined eigenvalue-estimate update pattern.
// This function is very similar to the previous methods, and differs by a shift in the matrix A by A-s*I, as well as using QR-decomposition to solve a linear system by backsubstitution instead of the

// We initialize our parameters
int n = A->size1, iter = 0, changed;
gsl_vector* V_new = gsl_vector_alloc(n);
gsl_matrix* A_QR = gsl_matrix_alloc(n,n);
gsl_matrix* Shift = gsl_matrix_alloc(n,n);
gsl_matrix* R_QR = gsl_matrix_alloc(n,n);
double denominator = 0, numerator = 0;
double Rayleigh = 0, NewRayleigh;

// We make a random vector, that will end up containing the eigenvector corresponding to the largest eigenvalue. Therefore we use V.
// We insert random numbers (between 0 and 10) in V
for (int i = 0; i < n; ++i) {
        gsl_vector_set(V, i, ((double) rand())/((double)RAND_MAX)*10);
}


// We want to calculate the new vector by (A−s*1)*x_i+1=x and solve the  Rayleigh quotient each time, to see when the Rayleigh quotient do not improve any more and hereby contains the eigenvalue.
// This is done with a loop:
do{
  // If we are asked to update every x'th iteration, being if mod(iteration) = x, then before the iteration we do:

  if (iter % update == 0) { // We choose mod(iter,update) to be 0 so we always do QR-decomposition at the beginning
    // We load the initial A-matrix
    gsl_matrix_memcpy(A_QR,A);
    // We shift the new Matrix A to A - I*s
    gsl_matrix_set_identity(Shift);
      if (iter != 0) {// Shift = s * I in the beginning
        s = 1/(Rayleigh) + s; // We redefine the shift with a shift along the found estimate
      }
      gsl_matrix_scale(Shift,s); // For the rest Shift = Rayleigh (estimate of e) * I
    gsl_matrix_sub(A_QR, Shift); //  A = A-Shift
    // We do the initial QR-decomposition
    qr_gs_decomp(A_QR, R_QR);
  }

  changed=0;
  iter++; // For each iteration, increase iter by one
  qr_gs_solve(A_QR, R_QR, V, V_new); // We solve the system (A−s*1) * X_i+1 = x using QR-decomposition and backsubstitution as done in the linear equation exercise

  // Rayleigh quotient is found λ[xi] =x^T_i * Ainv * x_i / x^T_i * x_i = x^T_i+1 * x_i / x^T_i * x_i
  gsl_blas_ddot(V, V, &(denominator)); // We calculate the denominator  x^T_i * x
  gsl_blas_ddot(V_new, V, &(numerator)); // We calculate the numerator x^T_i+1 * x_i
  // We estimate the  Rayleigh quotient
  NewRayleigh = (numerator)/(denominator);
  // If the Rayleigh quotient changes, we break out of the loop
  if (NewRayleigh!=Rayleigh) {
  changed=1;
  }
  // We update the values for the next iteration
  gsl_vector_memcpy(V,V_new); // V -> V_new
  Rayleigh = NewRayleigh;  // Rayleigh -> NewRayleigh
  } while(changed!=0); //We check if the operation had an effect, if not, we stop.

// We normalize the found eigenvector:
gsl_vector_scale(V,1./gsl_blas_dnrm2(V));

// After we have found the new eigenvector V and eigenvalue approximated through the inverse Rayleigh parameter plus the shift, we insert the eigenvalue to its vector at the largest place and return the value:
for (int i = 0; i < n; ++i) {
        gsl_vector_set(e, i, 0);
        if (i==0) {
          gsl_vector_set(e,i,1/Rayleigh+s);
        }
}

// We free our parameters
gsl_vector_free(V_new);
gsl_matrix_free(A_QR);
gsl_matrix_free(R_QR);
gsl_matrix_free(Shift);

// We return the number of iterations
return iter;
}
