#include <stdio.h>
#include <limits.h>
#include <float.h>
int equal(double a, double b, double tau, double epsilon);

int largest_integer=INT_MAX;
int minimal_integer=INT_MIN;
float eps_float = FLT_EPSILON;
double eps_dbl = DBL_EPSILON;
long double eps_ldbl = LDBL_EPSILON;
int i=1;


int main(){

// Part I i
 printf("Part I i.\n");
 printf("Maximal integer defined by the system %d\n", largest_integer);

 while (i+1>i) {i++;}
 printf("My max integer with while = %i\n",i);

for (i = 0; i < i+1; i++) {}
printf("My max integer with for = %i\n",i);

i = 0;
do {i++;} while(i < i+1);
printf("My max integer with do while = %i\n",i);
printf("All of these values are the same as expected\n\n");
// Part I ii
printf("Part I ii\n");
printf("Minimal integer defined by the system %d\n", minimal_integer);

i = 0;
while (i-1<i) {i--;}
printf("My min integer with while = %i\n",i);


for (i = 0; i-1 < i; i--) {}
printf("My min integer with for = %i\n",i);

i = 0;
do {i--;} while(i-1 < i);
printf("My min integer with do while = %i\n",i);
printf("All of these values are the same as expected\n\n");

// Part I iii
printf("Part I iii\n");
float x=1; while(1+x!=1){x/=2;} x*=2;
float xx=1; for (xx = 1; 1+xx!=1; xx/=2){} xx*=2;
float xxx=1; do {xxx/=2;} while(1+xxx!=1);  xxx*=2;
printf("The mashine epsilon for a floating is = %g\n", eps_float);
printf("My found mashine epsilon for a float using while is = %g\n", x);
printf("My found mashine epsilon for a float using for is = %g\n", xx);
printf("My found mashine epsilon for a float using do-while is = %g\n\n", xxx);

double u=1; while(1+u!=1){u/=2;} u*=2;
double uu=1; for (uu = 1; 1+uu!=1; uu/=2){} uu*=2;
double uuu=1; do {uuu/=2;} while(1+uuu!=1);  uuu*=2;
printf("The mashine epsilon for a double is = %g\n", eps_dbl);
printf("My found mashine epsilon for a double is = %g\n", u);
printf("My found mashine epsilon for a double using for is = %g\n", uu);
printf("My found mashine epsilon for a double using do-while is = %g\n\n", uuu);



long double y=1; while(1+y!=1){y/=2;} y*=2;
long double yy=1; for (yy = 1; 1+yy!=1; yy/=2){} yy*=2;
long double yyy=1; do {yyy/=2;} while(1+yyy!=1);  yyy*=2;
printf("The mashine epsilon for a long double is = %Lg\n", eps_ldbl);
printf("My found mashine epsilon for a long double is = %Lg\n", y);
printf("My found mashine epsilon for a long double using for is = %Lg\n", yy);
printf("My found mashine epsilon for a long double using do-while is = %Lg\n\n", yyy);
printf("All of these values are found to be the expected.\n\n");

// Part II i
int max=INT_MAX/3;
float a=1.0;
float sum_up_float = 0;
float sum_down_float = 0;
for (i = 1.0; i < max; i++) {sum_up_float=sum_up_float+a/i;};
for (i = 0.0; i < max; i++) {sum_down_float=sum_down_float+a/(max-i);};
printf("The sum will run to the term numver of %d\n", max);
printf("The sum up float will hereby give = %g\n", sum_up_float);
printf("The sum down float will hereby give = %g\n", sum_down_float);
// Part II iI
printf("The difference is due to the fact that the downwards sum will add small terms first, and after a long time be neglected \n");
printf("In the upwards sum we will start with high sums and more quickly neglect the small terms\n\n");
// Part II iii
printf("The sum of the reciprocal natural numbers to the maximal integer will converge, but not the sum to infinity which diverges. \n\n");
// Part II iv
double aa=1.0;
double sum_up_double = 0;
double sum_down_double = 0;
for (i = 1.0; i < max; i++) {sum_up_double=sum_up_double+aa/i;};
for (i = 0.0; i < max; i++) {sum_down_double=sum_down_double+aa/(max-i);};
printf("The sum up double will hereby give = %g\n", sum_up_double);
printf("The sum down double will hereby give = %g\n", sum_down_double);
printf("The sum now does not go to small enough terms to compare with the accuracy of the double \n\n");

// Part III
printf("We check if the function works by three tests one for tau, one for epsilon and one for negative\n");
double r = 1;
double b = 4;
double tau = 5;
double epsilon = 0.5;
int result = equal(r,b,tau,epsilon);
printf("The first test which should give a 1 is %d\n",result);
r = 1;
b = 3;
tau = 1;
epsilon = 2;
int result2 = equal(r,b,tau,epsilon);
printf("The second test which should give a 1 is %d\n",result2);

r = 1;
b = 600;
tau = 1;
epsilon = 1;
int result3 = equal(r,b,tau,epsilon);
printf("The third test which should give a 0 is %d\n",result3);
return 0;}
