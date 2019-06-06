void gkl_biadiag(gsl_matrix* U, gsl_matrix* A, gsl_matrix* V);
void gkl_solve(gsl_matrix* U, gsl_matrix* A, gsl_matrix* V, const gsl_vector* s, gsl_vector* x);
void glk_inverse(gsl_matrix* U, gsl_matrix* A, gsl_matrix* V, gsl_matrix* B);
void glk_backsub(gsl_matrix* R, gsl_vector* x);
void qr_gs_decomp(gsl_matrix*, gsl_matrix*);
void qr_gs_solve(gsl_matrix*,gsl_matrix*,const gsl_vector*,gsl_vector*);
void qr_gs_backsub(gsl_matrix*, gsl_vector*);
