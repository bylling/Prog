typedef struct {int number_of_hidden_neurons; double (*activation_function)(double); gsl_vector* data;} ann;  // We define the struct ann to hold our network, where  n is the number of neurons in the hidden layer, f is the pointer to the activation function, and gsl_vector* data keeps the parameters {ai, bi, wi}i=1..n.
ann* ann_alloc(int number_of_hidden_neurons, double(*activation_function)(double));
void ann_free(ann* network);
double ann_feed_forward(ann* neural_network, double x);
void ann_train(ann* neural_network, gsl_vector* xlist, gsl_vector* ylist);
double activation_function(double x);
double step_function(double x);
void matrix_print(const char* s, gsl_matrix* A);
void vector_print(const char* s, gsl_vector* v);
double ann_feed_forward_deriv(ann* neural_network, double x, double (*activation_function_deriv)(double));
double ann_feed_forward_antideriv(ann* neural_network, double x);
double activation_function_deriv(double x);
void ann_train_ODE(ann* neural_network, gsl_vector* xlist, double c, double yc, double ycmark);
double ODE_well(double x);
