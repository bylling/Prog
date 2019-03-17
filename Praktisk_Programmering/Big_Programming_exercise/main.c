#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>

double myarctan(double);
double myarccot(double);
double myarctan2(double);
double myarccot2(double);
double myarctan3(double);
double myarctan0(double);
int main(int argc, char* argv[])
{
printf("Exercise 21: Arctangent and Arccotangent\n");
	double a = atof(argv[1]);
	double b = atof(argv[2]);
	double dx = atof(argv[3]);

	for (double x = a; x <= b; x += dx)
fprintf(stderr, "%g %g %g %g %g %g %g\n", x, myarctan(x), atan(x), myarccot(x), atan(1.0/x), myarctan2(x),myarccot2(x));
printf("The plot is made and used in report.pdf\n");
return 0;
}
