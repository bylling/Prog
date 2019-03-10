#include <stdio.h>
#include <complex.h>
#include <math.h>

//Variables for Part I
complex double gammares;
complex double besselres;
complex double sqrtres;
complex double expresi;
complex double exprespi;
complex double irese;

//Variables for Part II
float xfloat = 0.1111111111111111111111111111;
double xdouble = 0.1111111111111111111111111111;
long double  xlong = 0.1111111111111111111111111111L;




int main()
{ // PArt I
   gammares = tgamma(5);
   besselres = j1( 0.5);
   sqrtres = csqrt(-2);
   expresi = cexp(I);
   exprespi = cexp(I*M_PI);
   irese = cpow(I, M_E);
  printf("Part I\n");
  printf("Gamma of five = %g + %g i\n" ,creal(gammares),cimag(gammares));
  printf("Bessel first of a half = %g + %g i\n",creal(besselres), cimag(gammares));
  printf("Squareroot of minus two = %g + %g I\n",creal(sqrtres), cimag(sqrtres));
  printf("Exponential of i = %g + I * %g\n",creal(expresi), cimag(expresi));
  printf("Exponential of i pi = %g + I * %g\n",creal(exprespi), cimag(exprespi));
  printf("I to the power of e = %g + I * %g\n\n", creal(irese), cimag(irese));

// PArt II
  printf("Part II\n");
  printf("Significant digits for float %.25g\n",xfloat);
  printf("Significant digits for double %.25lg\n",xdouble);
  printf("Significant digits for long double %.25Lg\n",xlong);

return 0;}
