#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"neuralnetworkfunctions.h"


int main (void)
{

  // Exercise A
fprintf(stdout, "Exercise A has started\n\n" );
fprintf(stdout, "A simple artificial neural network have been implemented.\n\n" );
fprintf(stdout, "For all of the following calculations we use a gaussian activation-function.\n" );
fprintf(stdout, "As a test, we want to approximate a rectangular function with the value 1 between 0 and 1, and 0 everywhere else. \n" );

// We initialize our neural networks
int number_of_hidden_neurons=2;
int number_of_hidden_neurons_2=7;
int number_of_hidden_neurons_guess=2;
ann* neural_network=ann_alloc(number_of_hidden_neurons,activation_function);
ann* neural_network_2=ann_alloc(number_of_hidden_neurons_2,activation_function);
ann* neural_network_guess=ann_alloc(number_of_hidden_neurons_guess,activation_function);
fprintf(stdout, "We do the training and calculation for a network consisting of %i and %i neurons. For them both, we use the same inital guess and the same number of trainingpoints. \n\n",number_of_hidden_neurons,number_of_hidden_neurons_2 );

// We initialize the training parameters
int number_of_training_points = 500;
gsl_vector* xlist=gsl_vector_alloc(number_of_training_points);
gsl_vector* ylist=gsl_vector_alloc(number_of_training_points);
double int_start = -5.;
double int_stop = 5.;
double fitfunction;
double x;
double deltaint = (int_stop - int_start)/((double)number_of_training_points-1.00);

// We make a bunch of points between -5 and 5 for training the neural_network.
for(int i = 0; i < number_of_training_points; i++) {
  // We make the x-values
  x = int_start + i * deltaint;
  // We insert in the known stepfunction as a test
  fitfunction = step_function(x);
  // We insert the values in vectors
  gsl_vector_set(xlist,i,x);
  gsl_vector_set(ylist,i,fitfunction);
}

// We insert a guess for the paramers in the neural network, for it to start the optimization.
for(int i = 0; i < number_of_hidden_neurons; i++) {
  // We insert guesses for the parameters in the parameter vectors which is structed as: a1,a2,a3...an,b1,b2,b3 ...bn,w1,w2,w3...wn
  gsl_vector_set(neural_network->data,i,-1.);
  gsl_vector_set(neural_network->data,number_of_hidden_neurons+i,1);
  gsl_vector_set(neural_network->data,number_of_hidden_neurons*2+i,2./(double)number_of_hidden_neurons);
}
for(int i = 0; i < number_of_hidden_neurons_2; i++) {
  // We insert guesses for the parameters in the parameter vectors which is structed as: a1,a2,a3...an,b1,b2,b3 ...bn,w1,w2,w3...wn
  gsl_vector_set(neural_network_2->data,i,-1.);
  gsl_vector_set(neural_network_2->data,number_of_hidden_neurons_2+i,1);
  gsl_vector_set(neural_network_2->data,number_of_hidden_neurons_2*2+i,2./(double)number_of_hidden_neurons_2);
}
for(int i = 0; i < number_of_hidden_neurons_guess; i++) {
  // We insert guesses for the parameters in the parameter vectors which is structed as: a1,a2,a3...an,b1,b2,b3 ...bn,w1,w2,w3...wn
  gsl_vector_set(neural_network_guess->data,i,-1.);
  gsl_vector_set(neural_network_guess->data,number_of_hidden_neurons_guess+i,1);
  gsl_vector_set(neural_network_guess->data,number_of_hidden_neurons_guess*2+i,2./(double)number_of_hidden_neurons_guess);
}

fprintf(stdout, "We start the training of the neural networks, using %i trainingpoints.\n\n",number_of_training_points);
// we call the training funciton
ann_train(neural_network,xlist,ylist);
ann_train(neural_network_2,xlist,ylist);
fprintf(stdout, "We training has ended.\n\n");

// After the training and in order to show the trained neural network graphically, we insert various points into the interval and plot them together with the analytical stepfunction
int number_of_plot_points = 1000;
double plot_start = -5.;
double plot_stop = 5.;
double x_plot;
double y_plot, y_plot_2,y_plot_guess;
double y_theory;
double deltaplot = (plot_stop - plot_start)/((double)number_of_plot_points-1.00);
FILE* file = fopen("plotdata1.txt", "w");

for (int i = 0; i < number_of_plot_points; i++) {
  // We make the x-values
  x_plot = plot_start + i * deltaplot;
  // We calculate the function based on the progression through the layer of neurons.
  y_plot = ann_feed_forward(neural_network,x_plot);
  y_plot_2 = ann_feed_forward(neural_network_2,x_plot);
  y_plot_guess = ann_feed_forward(neural_network_guess,x_plot);
  y_theory = step_function(x_plot);
  fprintf(file, "%g %g %g %g %g\n",x_plot,y_plot,y_theory,y_plot_2,y_plot_guess);
}

fprintf(stdout, "A comparison between the two neural networks with their common guess-values are shown in the corresponding figure.\n");
fprintf(stdout, "We see the fact that the neural network with 7 hidden neurons works better at approximating the step-function, than the neural network with 2 hidden neurons. This is due to the fact that we can make an approximate stepfunction on a basis of gaussian functions, which gets better and better the more gaussians we include in the function.\n");
fprintf(stdout, "With the demonstration of the implementation of a simple artificial neural network with corresponding traning-procedure, the first exercise is hereby done.\n\n\n");

// Exercise B
fprintf(stdout, "Exercise B has started\n\n" );
fprintf(stdout, "A derivative and antiderivative feed-forward function has been implemented, and will now be tested on the stepfunction as before. The functions have been implemented such, that the user would have to define the analytical deriavtive and antiderivative of the activationfunction, which in this case is the gaussian derivative and indefinite integration. However this do not change the fact, that based on the proper training of the neural network, we will be able to approximate any function, derivative and antiderivative on a gaussian basis. The error must then be considered with respect to how good the gaussian basis is in the specific case, based on knowledge of the specific function.\n\n");

fprintf(stdout, "We calculate the deriavtive and the antiderivative using the 7 Hidden neurons and 500 Training points.\n\n");

// We allocate the parameters
double y_plot_2_deriv, y_plot_2_antideriv;
FILE* file2 = fopen("plotdata2.txt", "w");

for (int i = 0; i < number_of_plot_points; i++) {
  // We make the x-values
  x_plot = plot_start + i * deltaplot;
  // We calculate the function based on the progression through the layer of neurons.
  y_plot_2 = ann_feed_forward(neural_network_2,x_plot);
  y_plot_2_deriv = ann_feed_forward_deriv(neural_network_2,x_plot,activation_function_deriv);
  y_plot_2_antideriv = ann_feed_forward_antideriv(neural_network_2,x_plot);
  y_theory = step_function(x_plot);
  fprintf(file2, "%g %g %g %g %g\n",x_plot,y_plot_2,y_theory,y_plot_2_deriv,y_plot_2_antideriv);
}
fprintf(stdout, "The calculation as ended.\n\n");
fprintf(stdout, "The analytical stepfunction, stepfunction made with the neural network, derivative and antiderivate of the stepfunction found by the neural network is plotted in the corresponding figure. \n" );

fprintf(stdout, "We find, that the derivate, which should be two delta-function are poorly described in the gaussian basis, however we see, that the sign changes as expected. \n We furthermore see that the anti-derivative is not represented in a good way either. As expected the function increases during the square function, but the wigles at the border of the gausian fit, makes the antiderivate return to a value of 0 instead of increasing to 1, as should be the case for the square function. Theese indescreptancies are hereby not due to the implementation of the artificial neural network, but based on the gausian basis used as the activation function.\n\n" );

fprintf(stdout, "With the demonstration of the implementation of a derivative and antiderivative approximation using the trained artificial neural network, the second exercise is hereby done.\n\n\n");




// Exercise C
fprintf(stdout, "Exercise C has started\n\n" );
fprintf(stdout, "A new traning algorithm have been implemented to train for the solving a differential equation using a given boundary condition.\n\n" );

// We initialize our neural networks
int number_of_hidden_neurons_ode=2;
ann* neural_network_ode=ann_alloc(number_of_hidden_neurons_ode,ODE_well);


// We initialize the training parameters
int number_of_training_points_ode = 500;
gsl_vector* xlist_ODE=gsl_vector_alloc(number_of_training_points);
double int_start_ODE = 0.;
double int_stop_ODE = 1.;
double yc;
double yc_mark;
double x_ODE;
double deltaint_ode = (int_stop_ODE - int_start_ODE)/((double)number_of_training_points-1.00);

// We make a bunch of points between 0 and 1 for training the neural_network.
for(int i = 0; i < number_of_training_points_ode; i++) {
  // We make the x-values
  x_ODE = int_start_ODE + i * deltaint_ode;
  // We insert the values in a vector
  gsl_vector_set(xlist_ODE,i,x_ODE);
}

// We insert a guess for the paramers in the neural network, for it to start the optimization.
for(int i = 0; i < number_of_hidden_neurons_ode; i++) {
  // We insert guesses for the parameters in the parameter vectors which is structed as: a1,a2,a3...an,b1,b2,b3 ...bn,w1,w2,w3...wn
  gsl_vector_set(neural_network_ode->data,i,-1.);
  gsl_vector_set(neural_network_ode->data,number_of_hidden_neurons_ode+i,1);
  gsl_vector_set(neural_network_ode->data,number_of_hidden_neurons_ode*2+i,2./(double)number_of_hidden_neurons_ode);
}

// We set the boundary conditions
double c = 1./2.;
yc = 1.;
yc_mark = 0.;

fprintf(stdout, "We start the training of the neural networks, with the same initial function as in exercise A.\n\n");
// we call the training funciton
ann_train_ODE(neural_network_ode,xlist_ODE,c,yc,yc_mark);
fprintf(stdout, "We training has ended.\n\n");

// After the training and in order to show the trained neural network graphically
int number_of_plot_points_ode = 1000;
double plot_start_ode = -5.;
double plot_stop_ode = 5.;
double x_plot_ode;
double y_plot_ode, ode_ode;
double deltaplot_ode = (plot_stop_ode - plot_start_ode)/((double)number_of_plot_points_ode-1.00);
FILE* file3 = fopen("plotdata3.txt", "w");

for (int i = 0; i < number_of_plot_points_ode; i++) {
  // We make the x-values
  x_plot_ode = plot_start_ode + i * deltaplot_ode;
  // We calculate the function based on the progression through the layer of neurons.
  y_plot_ode = ann_feed_forward(neural_network,x_plot_ode);
  ode_ode = ODE_well(x_plot_ode);
  fprintf(file3, "%g %g %g \n",x_plot_ode,y_plot_ode,ode_ode);
}
fprintf(stdout, "The finite well potential and a solution to such with constant energy is made with the neural network and is plotted in the corresponding figure. \n" );
fprintf(stdout, "We find that the gaussian basis seems to fulfill the boundary conditions at c = 1/2, with y' = 0. The gaussian basis is resonable to express the bound harmonic state in the well potential. The constants are not scaled, making the axis have arbitrary values, since this is only a simple demonstration.\n" );
fprintf(stdout, "With the proper values inserted, one would however be able to solve any differential equation using this method.\n" );
fprintf(stdout, "With the demonstration of the implementation of a differential equation solver using a specially trained artificial neural network, the third exercise is hereby done.\n\n\n");

// We free the parameters
ann_free(neural_network);
ann_free(neural_network_2);
ann_free(neural_network_guess);
ann_free(neural_network_ode);
gsl_vector_free(xlist);
gsl_vector_free(ylist);
fclose(file);
fclose(file2);
fclose(file3);
}
