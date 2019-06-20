#include<gsl/gsl_sf.h>
#include"make_airy_plot.h"
#include<math.h>
void make_airy_plot(void){
	for(double x=-15;x<5;x+=0.05){
		int mode = 1;
		double g=gsl_sf_airy_Ai(x,mode);
		double t=gsl_sf_airy_Bi(x,mode);
		fprintf(stderr,"%g %g %g\n",x,g,t);
	}
}
