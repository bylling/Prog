#include<stdio.h>
#include<math.h>
double myorbit(double,double,double);

int main()
{
	double epsilon_1 = 0;
	double uprime_1 = 0;
	double epsilon_2 = 0;
	double uprime_2 = -0.5;
	double epsilon_3 = 0.01;
	double uprime_3 = -0.5;
	for (double x = (0); x < 20 * M_PI; x += 0.05)
		fprintf(stderr, "%g %g %g %g\n", x, myorbit(x,epsilon_1,uprime_1),myorbit(x,epsilon_2,uprime_2),myorbit(x,epsilon_3,uprime_3));
printf("The Plot is made as plot.svg \n");	
return 0;
}
