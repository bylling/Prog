double adapt24(double f(double), double a, double b, double acc, double eps, double * err, double f2, double f3, int nrec);
double adapt(double f(double), double a, double b, double acc, double eps, double * err);
double adapt_clenshaw_curtis(double f(double), double a, double b, double acc, double eps, double* err);
double adapt_from_a_to_infinite(double f(double), double a, double acc, double eps, double* err);
