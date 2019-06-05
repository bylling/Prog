#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"functionfile.h"

int main (void)
{
  // Initial text:
fprintf(stdout, "\n\n ------------- The Final Exam -----------------\n" );
fprintf(stdout, "Exercise XX: Title \n\n" );
fprintf(stdout, "The Exercise will now be made, and the analysis as well as results of the method are described in the report.pdf file.\n" );

// We initialize our values
double x_plot, x_step = 1., y_plot, y_theory;
int n = 3;
gsl_matrix* I = gsl_matrix_alloc(n,n);
gsl_matrix_set_identity(I);

// We print the values
FILE* file = fopen("plot1.data", "w");
for (int i = 0; i < 10.; i++) {
  // We make the x-values
  x_plot = i*x_step;
  // We make the y-values.
  y_plot = x_plot * x_plot;
  y_theory = x_plot * x_plot;
  fprintf(file, "%g %g %g\n",x_plot,y_plot,y_theory);
}
matrix_print("Test x=", I);

// We free the parameters
gsl_matrix_free(I);
fclose(file);

return 0;
}
