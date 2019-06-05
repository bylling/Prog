#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"functionfile.h"

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
