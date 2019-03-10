#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"nvector.h"

nvector* nvector_alloc(int n){
  nvector* v = malloc(sizeof(nvector));
  (*v).size = n;
  (*v).data = malloc(n*sizeof(double));
  if( v==NULL ) fprintf(stderr,"error in nvector_alloc\n");
  return v;
}

void nvector_free(nvector* v){
  free(v->data);
  free(v);
} /* v->data is identical to (*v).data */

void nvector_set(nvector* v, int i, double value){
  (*v).data[i]=value;
}

double nvector_get(nvector* v, int i){
  return (*v).data[i];
}

double   nvector_dot_product (nvector* u, nvector* v){
double result = 0;
 for (int i = 0; i  < u->size; i++) {
 result = result + (*v).data[i] * (*u).data[i];
}
return result;
}

void nvector_print(char *s, nvector * v)
{
	printf("%s", s);
	for (int i = 0; i < v->size; i++)
		printf("%9.3g ", v->data[i]);
	printf("\n");
}

void nvector_set_zero (nvector* v){             /* all elements ← 0 */
    for (int i = 0; i < v->size; i++) {
    (*v).data[i]=0;
}
}

int double_equal(double a, double b)
{
	double TAU = 1e-6, EPS = 1e-6;
	if (fabs(a - b) < TAU)
		return 1;
	if (fabs(a - b) / (fabs(a) + fabs(b)) < EPS / 2)
		return 1;
	return 0;
}


  int nvector_equal(nvector * a, nvector * b)
  {
  	if (a->size != b->size) return 0;
  	for (int i = 0; i < a->size; i++)
  		if (!double_equal(a->data[i], b->data[i]))
  			return 0;
  	return 1;
  }

void nvector_add      (nvector* a, nvector* b){ /* a_i ← a_i + b_i */
  for (int i = 0; i < a->size; i++) {
		double s = nvector_get(a, i) + nvector_get(b, i);
		nvector_set(a, i, s);
  }
}


void nvector_sub      (nvector* a, nvector* b){ /* a_i ← a_i - b_i */
  for (int i = 0; i < a->size; i++) {
    double dif = nvector_get(a,i) - nvector_get(b,i);
    nvector_set(a,i,dif);
  }}

void nvector_scale    (nvector* a, double x){   /* a_i ← x*a_i     */
  for (int i = 0; i < a->size; i++) {
    double prod = nvector_get(a,i)*x ;
    nvector_set(a,i,prod);
    }
}
