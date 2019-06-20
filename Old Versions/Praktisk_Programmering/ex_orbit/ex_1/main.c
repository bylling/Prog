#include<stdio.h>
#include<math.h>
double mylog(double);

int main()
{
	for (double x = (0); x < 3; x += 0.1)
		fprintf(stderr, "%g %g %g\n", x, mylog(x), exp(x)/(exp(x)+1));
printf("The plot is made as plot.svg\n");	
return 0;
}
