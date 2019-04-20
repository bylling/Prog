#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<stdio.h>
#include"lineqfunctions.h"

// Theese functions are taken direcly from the exercise on linear equaitons.

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R){
// QR-decomposition of matrix A, using Graham Schmidt, Here A is replaced by Q and R is created
  int m = A->size2; //We get the number of columns
  for(int i=0;i<m;i++){
    gsl_vector_view e = gsl_matrix_column(A,i);  //We get each collumn of A as a temperary subset of the matrix
  	double r = gsl_blas_dnrm2(&e.vector); // This functions compute the norm of the vector
  	gsl_matrix_set(R,i,i,r); // We instert the scaling factor in R, so we can multiply Q and R to find the initial A
  	gsl_vector_scale(&e.vector,1/r); //normalization using the found norm

    for(int j=i+1;j<m;j++){ //For the rest of the columns:
  		gsl_vector_view q = gsl_matrix_column(A,j);//We get each collumn out
  		double s=0;  //We define projection value
      gsl_blas_ddot(&e.vector,&q.vector,&s);  //We project the new vector to the first
  		gsl_blas_daxpy(-s,&e.vector,&q.vector); // -s*e+q -> q  hereby modifying q to have no contribution in direction of the first vector
  		gsl_matrix_set(R,i,j,s); // We insert the projection of the new vector on the original on its cross term coupling in the R matrix
  		gsl_matrix_set(R,j,i,0); // We insert a zero on the transposed value, since we hereby decouple the new orthogonal basis.
  		}
  	}


}



void qr_gs_solve(gsl_matrix* Q, gsl_matrix* R, const gsl_vector* b, gsl_vector* x){
// As described we first multiply Q^T with b
gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x); // x = 1*Q^T*b + 0*x
// As described we call the inplace backsubstitution on the R and product.
qr_gs_backsub(R,x);

}

void qr_gs_backsub(gsl_matrix* R, gsl_vector* x){

	int m = R->size1; //Number of rows in R
// This is done as described in the lecture notes:
for (int i = m-1; i >= 0; i--){ // We start from the last collumbn
  double p=0;
  for (int k = i+1; k < m; k++){
    p = p + gsl_matrix_get(R, i, k) * gsl_vector_get(x, k); // We calculate the equation
  }
  gsl_vector_set(x, i, (gsl_vector_get(x, i)-p)/gsl_matrix_get(R, i, i)); //WE set the solution as defined for backsubstitution
}
}

void qr_gs_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
   gsl_matrix* I = gsl_matrix_alloc(Q->size1,Q->size2);
   gsl_matrix_set_identity(I);
   for (int i = 0; i < Q->size1; i++) {
       gsl_vector_view B_i = gsl_matrix_column(B, i);
       gsl_vector_const_view e_i = gsl_matrix_const_column(I, i);
       qr_gs_solve(Q, R, &(e_i.vector), &(B_i.vector));
   }
      gsl_matrix_free(I);
  }

void gkl_biadiag(gsl_matrix* U, gsl_matrix* A, gsl_matrix* V){
  // The procedure is implemented excactly like in the link in the exercise
gsl_vector *V_1=gsl_vector_alloc(A->size2);
gsl_vector *U_1=gsl_vector_alloc(A->size2);
for (int i = 0; i < V_1->size; ++i) {
        gsl_vector_set(V_1, i, 1);
}
double norm = 1/sqrt(V_1->size);
gsl_vector_scale(V_1,norm);
double beta = 0.0;
double alpha;
for (int k = 0; k < A->size2; k++) {
  gsl_matrix_set_col(V,k,V_1);
  gsl_blas_dgemv(CblasNoTrans,1.0,A,V_1,-beta,U_1); // u = 1*A*v - beta*u
  alpha = gsl_blas_dnrm2(U_1);
  gsl_vector_scale(U_1,1/alpha);
  gsl_matrix_set_col(U,k,U_1);
  gsl_blas_dgemv(CblasTrans,1.0,A,U_1,-alpha,V_1); // v = 1*A^T*u - alpha*b
  beta = gsl_blas_dnrm2(V_1);
  gsl_vector_scale(V_1,1/beta);
}
gsl_vector_free(V_1);
gsl_vector_free(U_1);
    }

void gkl_solve(gsl_matrix* U, gsl_matrix* A, gsl_matrix* V, const gsl_vector* s, gsl_vector* x){
  // At first we calculate B, and since solutions to Ax=S is the same as UBV^Tx=S
  // and since U and B are othhogonal it is B*x_mod=U^T S with x = V * xmod  we can easily by bachsubsisitution solve the equations
gsl_matrix *prod=gsl_matrix_alloc(4,4);
gsl_matrix *B=gsl_matrix_alloc(4,4);
gsl_vector *s_new=gsl_vector_alloc(4);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,U,A,0.0,prod); //prod =1 * U^T*A+0*prod
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,prod,V,0.0,B); //B =1 * prod*V+0*B
// We finc the U^T * S
gsl_blas_dgemv(CblasTrans,1.0,U,s,0.0,s_new); // s_new = 1*U^T*s - 0*s_enw

glk_backsub(B,s_new);
gsl_blas_dgemv(CblasNoTrans,1.0,V,s_new,0.0,x); //x = V * xmod + 0*x

gsl_matrix_free(prod);
gsl_matrix_free(B);
gsl_vector_free(s_new);
}


void glk_backsub(gsl_matrix* R, gsl_vector* x){

int m = R->size1; //Number of rows in R
gsl_vector_set(x, m-1, (gsl_vector_get(x, m-1)-0)/gsl_matrix_get(R, m-1, m-1)); // We do the last row, since this is special and contains 1 value

for (int i = m-2; i >= 0; i--){ // We start from the second last collumn
  double p=0;
  for (int k = i+1; k < i+2; k++){ // Now we only consider the two bidiagonal elements
    p = p + gsl_matrix_get(R, i, k) * gsl_vector_get(x, k); // We calculate the equation
  }
  gsl_vector_set(x, i, (gsl_vector_get(x, i)-p)/gsl_matrix_get(R, i, i)); //WE set the solution as defined for backsubstitution
}
}

double glk_determinant(gsl_matrix* U, gsl_matrix* A, gsl_matrix* V){
// This is done by calculating the determinant of B, U and V, to find A using determinant identities.
//We calculate B:
gsl_matrix *prod=gsl_matrix_alloc(4,4);
gsl_matrix *B=gsl_matrix_alloc(4,4);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,U,A,0.0,prod); //prod =1 * U^T*A+0*prod
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,prod,V,0.0,B); //B =1 * prod*V+0*B
// Since det(U*B*V^T) = det(U) , since u and v are orthogonal, they have determinant 1 or -1  so we only need to find B
// Det(B) = diagonal
double determ  = 1;
for (int i = 0; i < B->size1; ++i) {
    determ = determ * gsl_matrix_get(B, i, i);
}

gsl_matrix_free(prod);
gsl_matrix_free(B);
return determ;
}

void glk_inverse(gsl_matrix* U, gsl_matrix* A, gsl_matrix* V, gsl_matrix* B){
// This inverse is found just as before
   gsl_matrix* I = gsl_matrix_alloc(A->size1,A->size2);
   gsl_matrix_set_identity(I);
   for (int i = 0; i < A->size1; i++) {
       gsl_vector_view B_i = gsl_matrix_column(B, i);
       gsl_vector_const_view e_i = gsl_matrix_const_column(I, i);
       gkl_solve(U,A,V, &(e_i.vector), &(B_i.vector));
   }
      gsl_matrix_free(I);
  }


void vector_print(const char* s, gsl_vector* v){
    printf("%s\n",s);
    for(int i=0;i<v->size;i++){
      printf("%8.3g",gsl_vector_get(v,i));
      printf("\n");
    }
    printf("\n");
  }

  void matrix_print(const char* s, gsl_matrix* A){
    printf("%s\n",s);
    for(int i=0;i<A->size1;i++){
      for(int j=0;j<A->size2;j++){
        printf("%8.3g",gsl_matrix_get(A,i,j));
        printf("\t");
      }
    printf("\n");
    }
  }
