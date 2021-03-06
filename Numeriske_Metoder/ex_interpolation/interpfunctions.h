typedef struct {int n; double *x,*y,*b,*c;} qspline;
typedef struct {int n; double *x,*y,*b,*c,*d;} cubic_spline;
double myinteglin(int n, double *x, double *y, double z);
double mylinterp(int n, double * x, double * y, double z);
cubic_spline* cubic_spline_alloc(int n,double* x,double* y);
double cubic_spline_eval(cubic_spline *s,double z);
void cubic_spline_free(cubic_spline *s);
double cubic_spline_derivative(cubic_spline *s, double z);
double cubic_spline_integral(cubic_spline *s, double z);
qspline* qspline_alloc(int n,double* x,double* y);
double qspline_eval(qspline *s,double z);
void qspline_free (qspline *s);
double qspline_derivative(qspline *s, double z);
double qspline_integral(qspline *s, double z);
