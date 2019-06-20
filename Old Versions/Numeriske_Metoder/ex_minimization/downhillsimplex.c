#include <math.h>
// All of this function is taken directly from the lecture notes on Downhill simplex

void reflection(double* highest, double* centroid, int dim, double* reflected) {
  // Reflection: the highest point is reflected against the centroid, phi → pre = pce + (pce − phi).
	for(int i=0; i<dim; i++){
    reflected[i] = 2*centroid[i]-highest[i];
  }
}

void expansion(double* highest, double* centroid, int dim, double* expanded) {
  //Expansion: the highest point reflects and then doubles its distance from the centroid, phi → pex = pce + 2(pce − phi).
	for (int i=0; i<dim; i++){
    expanded[i] = 3*centroid[i]-2*highest[i];
  }
}

void contraction(double* highest, double* centroid, int dim, double* contracted) {
//Contraction: the highest point halves its distance from the centroid, phi → pco = pce + 1/2(phi − pce).
  for (int i=0; i<dim; i++){
    contracted[i] = 0.5*centroid[i]+0.5*highest[i];
  }
}

void reduction(double** simplex, int dim, int lo) {
  //Reduction: all points, except for the lowest, move towards the lowest points halving the distance. pk6=lo → 1/2 (pk + plo).
	for (int k=0; k<dim+1; k++){
    if (k!=lo){
		  for (int i=0; i<dim; i++){
			   simplex[k][i] = 0.5*(simplex[k][i]+simplex[lo][i]);
       }
     }
   }
}

double distance(double* a, double* b, int dim) {
// Calculates the distance between a and b in an arbitrary dimension
  double s = 0;
	for(int i=0; i<dim; i++){
   s += pow(a[i]-b[i],2);
  }
  return sqrt(s);
}

double size(double** simplex,int dim) {
  // We check the longest distance between the first and last point.
	double s = 0;
	for (int k=1; k<dim+1; k++) {
		double dist=distance(simplex[0],simplex[k],dim);
		if (dist>s){
      s = dist;
    }
	}
	return s;
}


void simplex_update(double** simplex, double* f_values, int d, int* hi, int* lo, double* centroid){
 // We find the highest, lowest values an centroid points of the simplex
	*hi = 0;
	double highest = f_values[0];
	*lo = 0;
	double lowest = f_values[0];

	for (int k=1; k<d+1; k++) {
		double next = f_values[k];
		if (next>highest) {
      highest = next; *hi = k;
    }
		if (next<lowest) {
      lowest = next; *lo = k;
    }
	}

	for (int i=0; i<d; i++) {
		double sum = 0;
		for (int k=0; k<d+1; k++){
      if (k!=*hi){
         sum += simplex[k][i];
       }
		centroid[i] = sum/d;
    }
	}
}

void simplex_initiate(double (*fun)(double*), double** simplex, double* f_values, int d,int* hi, int* lo, double* centroid){
  // We initiate the simplex to take values for different functioncall-values. And call an update to find the relevant features of the simplex
	for (int k=0; k<d+1; k++){
    f_values[k] = fun(simplex[k]);
  }
	simplex_update(simplex,f_values,d,hi,lo,centroid);
}


int downhill_simplex(double F(double*), double** simplex, int d, double simplex_size_goal){
	int hi, lo, k = 0;
	double centroid[d], F_value[d+1], p1[d], p2[d];

  // We initialize the simplex function
  simplex_initiate(F,simplex,F_value,d,&hi,&lo,centroid);

  // We start the repeating procedure until we reach the simplex-goal
  while (size(simplex,d)>simplex_size_goal){
    // We find the highest, lowest values an centroid points of the simplex
    simplex_update(simplex,F_value,d,&hi,&lo,centroid);
    // We try a reflection
    reflection(simplex[hi],centroid,d,p1);

    double f_re = F(p1);

    if (f_re<F_value[lo]) {// if the reflection looks good: // try expansion
			expansion(simplex[hi],centroid,d,p2);
			double f_ex = F(p2);
			if (f_ex<f_re) { // If expansion worked: accept expansion
				for (int i=0; i<d; i++) {
					simplex[hi][i] = p2[i];
					F_value[hi]=f_ex;
				}
			}
			else { // If expansion did not work: Reject expansion and accept reflection
				for (int i=0; i<d; i++) {
					simplex[hi][i] = p1[i];
					F_value[hi] = f_re;
				}
			}
		}
		else { //  If the reflection happended to not be good
			if (f_re<F_value[hi]) { // It is relative okay, we accept reflection
				for (int i=0; i<d; i++) {
					simplex[hi][i] = p1[i];
					F_value[hi] = f_re;
				}
			}
			else { // If not sufficiant, we try contradiction
				contraction(simplex[hi],centroid,d,p1);
				double f_co=F(p1);
				if (f_co<F_value[hi]) { // If contradiction works we acceot contradiction
					for (int i=0;i<d;i++) {
						simplex[hi][i]=p1[i];
						F_value[hi]=f_co;
					}
				}
				else { // if neither reflection or contraction works we do reductions
					reduction(simplex,d,lo);
					simplex_initiate(F,simplex,F_value,d,&hi,&lo,centroid);
				}
			}
		}
		k++;
	}
	return k;
}
