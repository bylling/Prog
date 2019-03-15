#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>

double myerr(double);
int main(int argc, char* argv[])
{

	double a = atof(argv[1]);
	double b = atof(argv[2]);
	double dx = atof(argv[3]);

	for (double x = a; x <= b+1e-6; x += dx)
		fprintf(stderr, "%g %g %g\n", x, myerr(x), gsl_sf_erf(x));
printf("The plot is made as plot.pdf and used in report.pdf\n");
return 0;
}
