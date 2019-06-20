#include"komplex.h"
#include"stdio.h"
#include"math.h"
#define TINY 1e-6

int main(){
	komplex a = {1,2}, b = {3,4};
	printf("\ntesting komplex_add...\n");
	komplex r = komplex_add(a,b);
	komplex R = {4,6};
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a+b should   = ", R);
	komplex_print("a+b actually = ", r);

	if( komplex_equal(R,r,TINY,TINY) )
		printf("test 'add' passed :) \n");
	else
		printf("test 'add' failed: debug me, please... \n");


  printf("\ntesting komplex_sub...\n");
  komplex r_1 = komplex_sub(b,a);
  komplex R_1 = {2,2};
  komplex_print("a=",a);
  komplex_print("b=",b);

  komplex_print("a+b should   = ", R_1);
  komplex_print("a+b actually = ", r_1);
  if( komplex_equal(R_1,r_1,TINY,TINY) )
    printf("test 'sub' passed :) \n");
  else
    printf("test 'sub' failed: debug me, please... \n");

  printf("\ntesting komplex_set...\n");
  komplex r_2;
  komplex_set(&r_2,a.re,a.im);
  komplex R_2 = {1,2};
  komplex_print("Set should return  = ", R_2);
  komplex_print("Set is actually = ", r_2);
  if( komplex_equal(R_2,r_2,TINY,TINY) )
    printf("test 'set' passed :) \n");
  else
    printf("test 'set' failed: debug me, please... \n");

  printf("\ntesting komplex_new...\n");
  komplex r_3 = komplex_new(a.re,a.im);
  komplex R_3 = {1,2};
  komplex_print("New should return  = ", R_3);
  komplex_print("New is actually = ", r_3);
  if( komplex_equal(R_2,r_2,TINY,TINY) )
    printf("test 'New' passed :) \n");
  else
    printf("test 'New' failed: debug me, please... \n");


	printf("\ntesting komplex_mul...\n");
	komplex r_4 = komplex_mul(a,b);
	komplex R_4 = {3,8};
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a*b should   = ", R_4);
	komplex_print("a*b actually = ", r_4);

	if( komplex_equal(R_4,r_4,TINY,TINY) )
		printf("test 'mul' passed :) \n");
	else
		printf("test 'mul' failed: debug me, please... \n");

	printf("\ntesting komplex_div...\n");
	komplex r_5 = komplex_div(a,b);
	komplex R_5 = {0.33333333,0.5};
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a/b should   = ", R_5);
	komplex_print("a/b actually = ", r_5);

	if( komplex_equal(R_5,r_5,TINY,TINY) )
		printf("test 'mul' passed :) \n");
	else
		printf("test 'mul' failed: debug me, please... \n");



	printf("\ntesting komplex_conjugate...\n");
	komplex r_6 = komplex_conjugate(a);
	komplex R_6 = {1,-2};
	komplex_print("a=",a);
	komplex_print("Congugate of a should   = ", R_6);
	komplex_print("Congugate of a is actually = ", r_6);

	if( komplex_equal(R_6,r_6,TINY,TINY) )
		printf("test 'mul' passed :) \n");
	else
		printf("test 'mul' failed: debug me, please... \n");


	printf("\ntesting komplex_abs...\n");
	double r_7 = komplex_abs(a);
	double R_7 = {sqrt(5)};
	komplex_print("a=",a);
	printf("Abs of a should   = %g\n", R_7);
	printf("Abs of a is actually = %g\n", r_7);


	printf("\ntesting komplex_exp...\n");
	komplex r_8 = komplex_exp(a);
	komplex R_8 = {-1.13120438375,2.47172667200};
	komplex_print("a=",a);
	komplex_print("Exp of a should   = ", R_8);
	komplex_print("Exp of a is actually = ", r_8);

	if( komplex_equal(R_8,r_8,TINY,TINY) )
		printf("test 'exp' passed :) \n");
	else
		printf("test 'exp' failed: debug me, please... \n");

	printf("\ntesting komplex_sin...\n");
	komplex r_9 = komplex_sin(a);
	komplex R_9 = {3.1657785,1.9596010};
	komplex_print("a=",a);
	komplex_print("Sin of a should   = ", R_9);
	komplex_print("Sin of a is actually = ", r_9);

	if( komplex_equal(R_9,r_9,TINY,TINY) )
		printf("test 'sin' passed :) \n");
	else
		printf("test 'sin' failed: debug me, please... \n");

	printf("\ntesting komplex_cos...\n");
	komplex r_10 = komplex_cos(a);
	komplex R_10 = {2.932723007,3.05189779915};
	komplex_print("a=",a);
	komplex_print("Cos of a should   = ", R_10);
	komplex_print("Cos of a is actually = ", r_10);

	if( komplex_equal(R_10,r_10,TINY,TINY) )
		printf("test 'Cos' passed :) \n");
	else
		printf("test 'Cos' failed: debug me, please... \n");

	printf("\ntesting komplex_sqrt...\n");
	komplex r_11 = komplex_sqrt(a);
	komplex R_11 = {1.2720196495,0.78615137775742};
	komplex_print("a=",a);
	komplex_print("Sqrt of a should   = ", R_11);
	komplex_print("Sqrt of a is actually = ", r_11);

	if( komplex_equal(R_11,r_11,TINY,TINY) )
		printf("test 'Sqrt' passed :) \n");
	else
		printf("test 'Sqrt' failed: debug me, please... \n");







}
