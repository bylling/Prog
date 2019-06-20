#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"komplex.h"

void komplex_print (char *s, komplex z) {
	printf ("%s (%g,%g)\n", s, z.re, z.im);
}

komplex komplex_new (double x, double y) {
	komplex z = { x, y };
	return z;
}

void komplex_set (komplex* z, double x, double y) {
	(*z).re = x;
	(*z).im = y;
}

komplex komplex_add (komplex a, komplex b) {
	komplex result = { a.re + b.re , a.im + b.im };
	return result;
}

komplex komplex_sub (komplex a, komplex b) {
	komplex result = { a.re - b.re , a.im - b.im };
	return result;
}

int komplex_equal (komplex a, komplex b, double acc, double eps)  {
  if (fabs(a.re - b.re) < acc ) {return 1;}
  if (fabs(a.im - b.im) < acc ) {return 1;}
  if (fabs(a.re - b.re)/(fabs(a.re) + fabs(b.re)) < eps ) {return 1;}
  if (fabs(a.im - b.im)/(fabs(a.im) + fabs(b.im)) < eps ) {return 1;}
  else return 0;
}

komplex komplex_mul (komplex a, komplex b)  {
  komplex result = {a.re * b.re , a.im * b.im };
  return result;
}

komplex komplex_div (komplex a, komplex b)  {
  komplex result = {a.re / b.re , a.im / b.im };
  return result;
}

komplex komplex_conjugate (komplex z)  {
  komplex result = {z.re , -z.im };
  return result;
}

double komplex_abs (komplex z)  {
	double result = {sqrt(pow(fabs(z.re),2) + pow(fabs(z.im),2))};
  return result;
}

komplex komplex_exp (komplex z)  {
	komplex result = {exp(z.re) * cos(z.im), exp(z.re) * sin(z.im)};
  return result;
}

komplex komplex_sin (komplex z)  {
	komplex result = {sin(z.re)*cosh(z.im), cos(z.re)*sinh(z.im)};
  return result;
}

komplex komplex_cos (komplex z)  {
	komplex result = {cos(z.re)*cosh(z.im), sin(z.re)*sinh(z.im)};
  return result;
}

komplex komplex_sqrt (komplex z)  {
	komplex result = {sqrt((z.re + sqrt(pow(z.re,2) + pow(z.im,2)))/2) , (z.im/abs(z.im)) * sqrt((-z.re + sqrt(pow(z.re,2) + pow(z.im,2)))/2)};
  return result;
}
