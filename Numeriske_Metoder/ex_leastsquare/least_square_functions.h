double funs(int i, double x);
void lsfitA(double funs(int i, double x),int freeparams, gsl_vector* vec_xdata, gsl_vector* vec_ydata, gsl_vector* vec_dydata, gsl_vector* c, gsl_matrix* S);
void SVD(gsl_matrix* A, gsl_matrix* V, gsl_matrix* D, gsl_matrix* U,  gsl_matrix* S);
void linsysSVD(gsl_matrix* MatrixA,  gsl_matrix* MatrixV,  gsl_matrix* MatrixD,  gsl_matrix* MatrixU,  gsl_matrix* MatrixS2, gsl_vector* VectorB, gsl_vector* VectorX);
void lsfitC(double funs(int i, double x),int freeparams,gsl_vector* vec_xdata, gsl_vector* vec_ydata, gsl_vector* vec_dydata, gsl_vector* c, gsl_matrix* S );
