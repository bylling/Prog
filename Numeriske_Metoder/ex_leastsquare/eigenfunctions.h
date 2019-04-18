
int jacobi(gsl_matrix* A, gsl_vector* e, gsl_matrix* V);
int jacobi_mod_ev(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int eigennr);
int jacobi_mod_highev(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int eigennr);
int jacobi_versc(gsl_matrix* A, gsl_vector* e, gsl_matrix* V);
int find_max_in_row(gsl_matrix* A, int r);
