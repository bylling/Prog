#include "nvector.h"
#include <stdio.h>
#include <stdlib.h>
#define RND (double)rand()/RAND_MAX

int main()
{
	int n = 5;

	printf("\nmain: testing nvector_alloc ...\n");
	nvector *v = nvector_alloc(n);
	if (v == NULL) printf("test failed\n");
	else printf("test passed\n");

	printf("\nmain: testing nvector_set and nvector_get ...\n");
	double value = RND;
	int i = n / 2;
	nvector_set(v, i, value);
	double vi = nvector_get(v, i);
	if (double_equal(vi, value)) printf("test passed\n");
	else printf("test failed\n");

  printf("\nmain: testing nvector_add ...\n");
	nvector *a = nvector_alloc(n);
	nvector *b = nvector_alloc(n);
	nvector *c = nvector_alloc(n);
	for (int i = 0; i < n; i++) {
		double x = RND, y = RND;
		nvector_set(a, i, x);
		nvector_set(b, i, y);
		nvector_set(c, i, x + y);
	}
	nvector_add(a, b);
	nvector_print("a+b should   = ", c);
	nvector_print("a+b actually = ", a);

	if (nvector_equal(c, a))
		printf("test passed\n");
	else
		printf("test failed\n");



  printf("\nmain: testing nvector_dot_product ...\n");
	nvector *a_1 = nvector_alloc(n);
	nvector *b_1 = nvector_alloc(n);
  double result = 0;
  for (int i = 0; i < n; i++) {
		double x = RND, y = RND;
		nvector_set(a_1, i, x);
		nvector_set(b_1, i, y);
		result = result + x * y;
	}
	double c_1 = nvector_dot_product(a_1, b_1);
	printf("a+b should   = %g\n", result);
	printf("a+b actually = %g\n", c_1);
	if (double_equal(result, c_1))
		printf("test passed\n");
	else
		printf("test failed\n");




  printf("\nmain: testing nvector_set_zero ...\n");
	nvector *a_2 = nvector_alloc(n);
	nvector *b_2 = nvector_alloc(n);
  for (int i = 0; i < n; i++) {
		double x = RND, y = 0;
		nvector_set(a_2, i, x);
		nvector_set(b_2, i, y);
	}
  nvector_set_zero(a_2);
	nvector_print("Zero should   = ", b_2);
	nvector_print("Zero actually = ", a_2);
	if (nvector_equal(b_2, a_2))
		printf("test passed\n");
	else
		printf("test failed\n");



  printf("\nmain: testing nvector_sub ...\n");
	nvector *a_3 = nvector_alloc(n);
	nvector *b_3 = nvector_alloc(n);
	nvector *c_3 = nvector_alloc(n);
	for (int i = 0; i < n; i++) {
		double x = RND, y = RND;
		nvector_set(a_3, i, x);
		nvector_set(b_3, i, y);
		nvector_set(c_3, i, x - y);
	}
	nvector_sub(a_3, b_3);
	nvector_print("a-b should   = ", c);
	nvector_print("a-b actually = ", a);

	if (nvector_equal(c_3, a_3))
		printf("test passed\n");
	else
		printf("test failed\n");



  printf("\nmain: testing nvector_scale ...\n");
	nvector *a_4 = nvector_alloc(n);
	nvector *b_4 = nvector_alloc(n);
  double y_1 = RND;
  for (int i = 0; i < n; i++) {
		double x = RND;
		nvector_set(a_4, i, x);
		nvector_set(b_4, i, x * y_1);
	}
	nvector_scale(a_4, y_1);
	nvector_print("a scaled y should   = ", b_4);
	nvector_print("a scaled y actually = ", a_4);

	if (nvector_equal(a_4, b_4))
		printf("test passed\n");
	else
		printf("test failed\n");


	nvector_free(v);
	nvector_free(a);
	nvector_free(b);
	nvector_free(c);
  nvector_free(a_1);
	nvector_free(b_1);
  nvector_free(a_2);
	nvector_free(b_2);
  nvector_free(a_3);
	nvector_free(b_3);
	nvector_free(c_3);
  nvector_free(b_4);
	return 0;
}
