void randomx(gsl_vector* a, gsl_vector* b, gsl_vector* x);
int plainmc(gsl_vector* a, gsl_vector* b, double f(gsl_vector*), int N, double* result, double* error);
int recursivemc(gsl_vector* a, gsl_vector* b, double f(gsl_vector*), int N, double* result, double* error, double acc, double eps);
double findvolume(gsl_vector* a, gsl_vector * b);
int randommontecarlo(gsl_vector* a, gsl_vector* b, double f(gsl_vector*), int N, double * result, double * var);
void matrix_print(const char* s, gsl_matrix* A);
void vector_print(const char* s, gsl_vector* v);
