void matrix_print(const char* s, gsl_matrix* A);
void vector_print(const char* s, gsl_vector* v);
int power_iteration(gsl_matrix* A, gsl_vector* e, gsl_vector* V);
int inverse_power_iteration(gsl_matrix* A, gsl_vector* e, gsl_vector* V);
int shifted_inverse_iteration(gsl_matrix* A, gsl_vector* e, gsl_vector* V,double s);
