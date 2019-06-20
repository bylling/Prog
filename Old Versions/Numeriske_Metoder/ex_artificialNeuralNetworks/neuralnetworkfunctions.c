#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_multimin.h>
#include"neuralnetworkfunctions.h"

ann* ann_alloc(int number_of_hidden_neurons, double(*activation_function)(double)){
  // This function should allocate memory for the struct, as well as intialize the struct.
  // The struct is allocated as in the optimization-exercise in Practical Programming.

// We allocate memory for the neural network in the struct with pointer ann*
  ann* neural_network = malloc(sizeof(ann));
// We insert our variables from the allocating function into the neural_network:
  neural_network->number_of_hidden_neurons = number_of_hidden_neurons;
  neural_network->activation_function = activation_function;
// We also need a data-vector to hold the parameters of each neuron {ai, bi, wi}i=1..n.
// We have 3 parameters for each neuron, and therefore allocate the datavector
  neural_network->data= gsl_vector_alloc(3*number_of_hidden_neurons);
  return neural_network;
}


void ann_free(ann* network){
// This function is required to free the allocated struct with its data.

// We free the data-vector in the struct
gsl_vector_free(network->data);
// The rest of the values in the function is simply double, and are not needed to be freed. Therefore we just free the struct:
free(network);
}

double ann_feed_forward(ann* neural_network, double x){
// This function deals with all of the network. It takes the data x and through the identity neuron will send x to all of the hidden neurons.
// After acting with the function on x in all of the neurons with parameters from data, by y=f((x-ai)/bi)*wi, we sum all of the data in the summation neuron, and return the answer y = f(x)

// We initialize the parameters
double y = 0;
double a_i, b_i, w_i, y_i;
// We do the same operation on all of the neurons:
for (int i = 0; i < neural_network->number_of_hidden_neurons; i++) {
  // We pull out the wanted data from the data-vector, which is structed as: a1,a2,a3...an,b1,b2,b3 ...bn,w1,w2,w3...wn
  a_i = gsl_vector_get(neural_network->data,i);
  b_i = gsl_vector_get(neural_network->data,neural_network->number_of_hidden_neurons+i);
  w_i = gsl_vector_get(neural_network->data,neural_network->number_of_hidden_neurons*2+i);
  // We act with the function from the neural network as y=f((x-ai)/bi)*wi
  y_i = neural_network->activation_function((x-a_i)/b_i) * w_i/b_i;
  y += y_i;
}
return y;
}

double ann_feed_forward_deriv(ann* neural_network, double x, double (*activation_function_deriv)(double)){
// This function works just as the ann_feed_forward, but this function returns the derivative of the function. It takes the data in x, and through the trained neural_network returns the derivative.
// After acting with the known analytical derivative of the function on x in all of the neurons with parameters from data, by y'=f'((x-ai)/bi)*wi/bi, we sum all of the data in the summation neuron, and return the answer y' = d/dx f(x)
// Since we do not want to change all other functions, we do not implement the derivative in the neural_network-struct, instead we just read it through the function.

// We initialize the parameters
double y = 0;
double a_i, b_i, w_i, y_i;
// We do the same operation on all of the neurons:
for (int i = 0; i < neural_network->number_of_hidden_neurons; i++) {
  // We pull out the wanted data from the data-vector, which is structed as: a1,a2,a3...an,b1,b2,b3 ...bn,w1,w2,w3...wn
  a_i = gsl_vector_get(neural_network->data,i);
  b_i = gsl_vector_get(neural_network->data,neural_network->number_of_hidden_neurons+i);
  w_i = gsl_vector_get(neural_network->data,neural_network->number_of_hidden_neurons*2+i);
  // We acti with the derivative of the activation function from the neural network as y'=f'((x-ai)/bi)*wi/bi
  y_i = activation_function_deriv((x-a_i)/b_i) * w_i/b_i;
  y += y_i;
}
return y;
}

double ann_feed_forward_antideriv(ann* neural_network, double x){
// This function works just as the ann_feed_forward and ann_feed_forward_deriv, but this function returns the anti derivative of the function. It takes the data in x, and through the trained neural_network returns the anti derivative.
// After acting with the known analytical antiderivative form the gaussian basis  on x in all of the neurons with parameters from data, using a gaussian basis, we sum all of the data in the summation neuron, and return the answer  Y = int f(x) dx
// This function is specific for the gaussian basis, since antiderivative and derivatives are not equally simple to calculate.

// We initialize the parameters
double y = 0;
double a_i, b_i, w_i, y_i;
// We do the same operation on all of the neurons:
for (int i = 0; i < neural_network->number_of_hidden_neurons; i++) {
  // We pull out the wanted data from the data-vector, which is structed as: a1,a2,a3...an,b1,b2,b3 ...bn,w1,w2,w3...wn
  a_i = gsl_vector_get(neural_network->data,i);
  b_i = gsl_vector_get(neural_network->data,neural_network->number_of_hidden_neurons+i);
  w_i = gsl_vector_get(neural_network->data,neural_network->number_of_hidden_neurons*2+i);
  // We acti with the antiderivative gaussian function from the neural network as: Y = -1/2 sqrt(pi) b w erf((a-x)/b)
  y_i =-(1./2.)* sqrt(M_PI)*w_i*b_i* (erf((a_i-x)/b_i));
  y += y_i;
}
return y;
}




void ann_train(ann* neural_network, gsl_vector* xlist, gsl_vector* ylist){
// Given the tabulated function, {xk,yk}k=1..N, the training of the network consists of tuning its parameters to minimize the deviation
// δ(p)=∑k=1..N(Fp(xk)-yk)²,
// So we hereby input a vector of y-values with corresponding x-values.
// We want to minimize the difference between them, which we calculate as a seperate function, the mismatch- function
// The mismatch will be minimize through varying the parameters of the data-vector in the neural network, therfore we include the parameter vectro
  double mismatch(const gsl_vector* parameters, void * params){
    // This nested function calculates the mismatch  δ(p)=∑k=1..N(Fp(xk)-yk)² with given parameters
      // We insert the parameters into the neural network
      gsl_vector_memcpy(neural_network->data,parameters);
      // We initialize
      double delta=0;
      // For every traning point x
  		for(int i=0; i<xlist->size; i++){
        // We load out the parameters
  			double x=gsl_vector_get(xlist,i);
  			double y=gsl_vector_get(ylist,i);
        // We calculate the expected y-values based on the x-input and the paramers of the function
  			double f=ann_feed_forward(neural_network,x);
        // We sum up the mismatch ∑k=1..N(Fp(xk)-yk)²
  			delta=delta+pow(f-y,2);
  		}
      // We weight the mismatch by the number of traning points, so a large training set and a small set can be compared equally per point.
  		return delta/xlist->size;
      // And return the wanted missmatch
  	}
// After definding the function in which we want to train our system to minimize, this being the difference betwen a known set of points with correspomnding values, and the expected points based on the same values, we initialize
// We initialize the parameters
gsl_vector* parameters=gsl_vector_alloc(neural_network->data->size);
// We use the previous paramers of the neural network as the start-guess of the minimization
gsl_vector_memcpy(parameters,neural_network->data);

// We call the GSL_minimization algotirhm from the Minimization exercise in Practical Programming, to minimize the mismatch function by variyng the parameters
// this is used instead of the quasi-newton method, since my quasi-newton method was made for analytical known gradiant functions. Therefore we simply use the GSL method instead
// We initialize
gsl_multimin_function F;
F.f=mismatch;
F.n=neural_network->data->size;
F.params=NULL;

gsl_vector *step = gsl_vector_alloc(F.n);
gsl_vector_set_all(step,0.5);


gsl_multimin_fminimizer *state = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,F.n);
gsl_multimin_fminimizer_set (state, &F, parameters, step);
// We iterate
int iter=0, status;
double acc=0.01;
do{
		iter++;
		int iteration_status = gsl_multimin_fminimizer_iterate(state);
		if(iteration_status != 0){
			fprintf(stderr,"unable to improve\n");
			break;
    }
		status = gsl_multimin_test_size(state->size, acc);
	}while( status == GSL_CONTINUE && iter < 10000);


// The returned p-value is saved the neural network.
gsl_vector_memcpy(neural_network->data,state->x);

// The training has finished.
gsl_vector_free(step);
gsl_multimin_fminimizer_free(state);
gsl_vector_free(parameters);
}


void ann_train_ODE(ann* neural_network, gsl_vector* xlist, double c, double yc, double ycmark){
// Given the boundary conditions, the training of the network to solve a ODE consists of tuning its parameters to minimize the deviation to the boundary conditions as well as the energy in the solution between the boundary conditions.
// This is done as: δ(p)=∫ab|Φ[Fp,x]|²dx +|Fp(c)-yc|²(b-a) +|Fp'(c)-y'c|²(b-a) .
// We want to minimize the difference between the boundary and our fixed boundary, which we calculate as a seperate function, the mismatch- function
// The mismatch will be minimize through varying the parameters of the data-vector in the neural network, therfore we include the parameter vector
  double mismatch(const gsl_vector* parameters, void * params){
    // This nested function calculates the mismatch  δ(p)=∫ab|Φ[Fp,x]|²dx +|Fp(c)-yc|²(b-a) +|Fp'(c)-y'c|²(b-a) . with given parameters
      // We insert the parameters into the neural network
      gsl_vector_memcpy(neural_network->data,parameters);
      // We initialize
      double delta=0;
      // For every traning point x
  		for(int i=0; i<xlist->size; i++){
        // We load out the parameters
  			double x=gsl_vector_get(xlist,i);
        // We calculate the expected y-values based on the x-input and the paramers of the function
  			double f=ann_feed_forward(neural_network,x);

        // If we are at the boundary condition, we handle the problem seperately.
        if (x==c) {
          // We sum up the mismatch δ(p)=|Fp(c)-yc|²(b-a)
    			delta=delta+pow(f-yc,2)*(gsl_vector_max(xlist)-gsl_vector_min(xlist));
          // We sum up the mismatch of the derivative, which is chosen to be the gaussian
          double fderiv=ann_feed_forward_deriv(neural_network,x,activation_function_deriv);
          delta=delta+pow(fderiv-ycmark,2)*(gsl_vector_max(xlist)-gsl_vector_min(xlist));
        }
        // We find the anti-derivative
        double fantideriv = ann_feed_forward_antideriv(neural_network,x);
        // If we are not at the boundary condition, we minmize the antiderivate
  			delta=delta+pow(fantideriv,2);
  		}
      // We weight the mismatch by the number of traning points, so a large training set and a small set can be compared equally per point.
  		return delta/xlist->size;
      // And return the wanted missmatch
  	}
// After definding the function in which we want to train our system to minimize, this being the difference betwen a known set of points with correspomnding values, and the expected points based on the same values, we initialize
// We initialize the parameters
gsl_vector* parameters=gsl_vector_alloc(neural_network->data->size);
// We use the previous paramers of the neural network as the start-guess of the minimization
gsl_vector_memcpy(parameters,neural_network->data);

// We call the GSL_minimization algotrihm from the Minimization exercise in Practical Programming, to minimize the mismatch function by variyng the parameters
// this is used instead of the quasi-newton method, since my quasi-newton method was made for analytical known gradiant functions. Therefore we simply use the GSL method instead
// We initialize
gsl_multimin_function F;
F.f=mismatch;
F.n=neural_network->data->size;
F.params=NULL;
gsl_vector *step = gsl_vector_alloc(F.n);
gsl_vector_set_all(step,0.5);
gsl_multimin_fminimizer *state = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,F.n);
gsl_multimin_fminimizer_set (state, &F, parameters, step);
// We iterate
int iter=0, status;
double acc=0.01;
do{
		iter++;
		int iteration_status = gsl_multimin_fminimizer_iterate(state);
		if(iteration_status != 0){
			fprintf(stderr,"unable to improve\n");
			break;
    }
		status = gsl_multimin_test_size(state->size, acc);
	}while( status == GSL_CONTINUE && iter < 10000);


// The returned p-value is saved the neural network.
gsl_vector_memcpy(neural_network->data,state->x);

// The training has finished.
gsl_vector_free(step);
gsl_multimin_fminimizer_free(state);
gsl_vector_free(parameters);
}


double activation_function(double x){ // The gaussian activation function for the neural network
  return exp(-x*x);
}

double activation_function_deriv(double x){ // The gaussian activation function for the neural network
  return -2. * x * exp(-x*x);
}

double step_function(double x){ // The heavyside step_function which we want to fit, with 1 between 0 and 1 and 0 elsewhere
  if (0<x && x<1) {
    return 1.;
  }
  else{
    if (x==0 || x==1) {
      return 1./2.;
    }
  return 0.;
  }
}


double ODE_well(double x){ // The simple potential well in arbitrary units
  double potential = 0;
  double height = 1;
  if (0>x || x>1) {
    potential =  height;
  }
  else{
    if (x==0 || x==1) {
      potential = 1./2. * height;
    }
  }

  return  potential;

}


void vector_print(const char* s, gsl_vector* v){
    printf("%s\n",s);
    for(int i=0;i<v->size;i++){
      printf("%8.3g",gsl_vector_get(v,i));
      printf("\n");
    }
    printf("\n");
  }

void matrix_print(const char* s, gsl_matrix* A){
  printf("%s\n",s);
  for(int i=0;i<A->size1;i++){
    for(int j=0;j<A->size2;j++){
      printf("%8.3g",gsl_matrix_get(A,i,j));
      printf("\t");
    }
  printf("\n");
  }
}
