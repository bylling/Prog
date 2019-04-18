#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<stdio.h>
#include"least_square_functions.h"
#include"eigenfunctions.h"
#include"lineqfunctions.h"

// All of theese functions are taken direcly from the exercise on eigenfunctions, for use in the Single value decomposition
int jacobi(gsl_matrix* A, gsl_vector* e, gsl_matrix* V){ // This function is made just as defined in the lecture notes.
// Cyclic Jacobi diagonalization; upper triangle of A is destroyed, e and V accumulate eigenvalues and eigenvectors
int changed, sweeps=0, n=A->size1; // Initialised values
for(int i=0;i<n;i++)gsl_vector_set(e,i,gsl_matrix_get(A,i,i)); // We copy the diagonal of A to the eigenmatrix
gsl_matrix_set_identity(V); // We define identity matrix
do{
  changed=0;
  sweeps++; // For each iteration, increase sweeps by one
  int p,q;
	for(p=0;p<n;p++)for(q=p+1;q<n;q++){
		double app=gsl_vector_get(e,p); // The selected value  A_PP is picked out
		double aqq=gsl_vector_get(e,q); // The selected valie A_QQ picked out
		double apq=gsl_matrix_get(A,p,q); // The off diagonal element A_PQ = A_QP is picked out
		double phi=0.5*atan2(2*apq,aqq-app); // The angle to minimize the off diagonal element is found
		double c = cos(phi), s = sin(phi); // The rotation transformation is made
		double app1=c*c*app-2*s*c*apq+s*s*aqq; // The transformation acts on the eigenvalues
		double aqq1=s*s*app+2*s*c*apq+c*c*aqq; // The tramformaiton acts on the eigenvalues
		if(app1!=app || aqq1!=aqq){ changed=1; // If the eigenvalue doesnt change, we have reached the point in which we cannot improve
			gsl_vector_set(e,p,app1); // We set the new eigenvalues into the vectors
			gsl_vector_set(e,q,aqq1);
			gsl_matrix_set(A,p,q,0.0); // We eleminate the off diagonal element in the upper diagonal
			for(int i=0;i<p;i++){ // We act with the transformation on the each row
				double aip=gsl_matrix_get(A,i,p);
				double aiq=gsl_matrix_get(A,i,q);
				gsl_matrix_set(A,i,p,c*aip-s*aiq);
				gsl_matrix_set(A,i,q,c*aiq+s*aip); }
			for(int i=p+1;i<q;i++){ // We act with the transformation on rews between p and q
				double api=gsl_matrix_get(A,p,i);
				double aiq=gsl_matrix_get(A,i,q);
				gsl_matrix_set(A,p,i,c*api-s*aiq);
				gsl_matrix_set(A,i,q,c*aiq+s*api); }
			for(int i=q+1;i<n;i++){ // We act with the transformation on rows after q
				double api=gsl_matrix_get(A,p,i);
				double aqi=gsl_matrix_get(A,q,i);
				gsl_matrix_set(A,p,i,c*api-s*aqi);
				gsl_matrix_set(A,q,i,c*aqi+s*api); }
			for(int i=0;i<n;i++){ //We do the similar transformation on the identity to gain the eigenvectors
				double vip=gsl_matrix_get(V,i,p);
				double viq=gsl_matrix_get(V,i,q);
				gsl_matrix_set(V,i,p,c*vip-s*viq);
				gsl_matrix_set(V,i,q,c*viq+s*vip);}
  		}
    }
  } while(changed!=0); //We check if the operation had an effect, if not, we stop.
return sweeps; }

int jacobi_mod_ev(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int eigennr){ // This function is made just as defined previously, but one row at the time.
// Cyclic Jacobi diagonalization, but modified to go from equation to equation, and return only one eigenvalues and eigenvectors that are true
int changed, sweeps=0, n=A->size1; // Initialised values
for(int i=0;i<n;i++)gsl_vector_set(e,i,gsl_matrix_get(A,i,i)); // We copy the diagonal of A to the eigenmatrix
gsl_matrix_set_identity(V); // We define identity matrix
do{
  changed=0;
  sweeps++; // For each iteration, increase sweeps by one
  int p,q;
	for(p=0;p<eigennr;p++)for(q=p+1;q<n;q++){
		double app=gsl_vector_get(e,p); // The selected value  A_PP is picked out
		double aqq=gsl_vector_get(e,q); // The selected valie A_QQ picked out
		double apq=gsl_matrix_get(A,p,q); // The off diagonal element A_PQ = A_QP is picked out
		double phi=0.5*atan2(2*apq,aqq-app); // The angle to minimize the off diagonal element is found
		double c = cos(phi), s = sin(phi); // The rotation transformation is made
		double app1=c*c*app-2*s*c*apq+s*s*aqq; // The transformation acts on the eigenvalues
		double aqq1=s*s*app+2*s*c*apq+c*c*aqq; // The tramformaiton acts on the eigenvalues
		if(app1!=app || aqq1!=aqq){ changed=1; // If the eigenvalue doesnt change, we have reached the point in which we cannot improve
			gsl_vector_set(e,p,app1); // We set the new eigenvalues into the vectors
			gsl_vector_set(e,q,aqq1);
			gsl_matrix_set(A,p,q,0.0); // We eleminate the off diagonal element in the upper diagonal
			for(int i=0;i<p;i++){ // We act with the transformation on the each row
				double aip=gsl_matrix_get(A,i,p);
				double aiq=gsl_matrix_get(A,i,q);
				gsl_matrix_set(A,i,p,c*aip-s*aiq);
				gsl_matrix_set(A,i,q,c*aiq+s*aip); }
			for(int i=p+1;i<q;i++){ // We act with the transformation on the each row between q and p
				double api=gsl_matrix_get(A,p,i);
				double aiq=gsl_matrix_get(A,i,q);
				gsl_matrix_set(A,p,i,c*api-s*aiq);
				gsl_matrix_set(A,i,q,c*aiq+s*api); }
			for(int i=q+1;i<n;i++){ // We act with the transformation on the the last columns
				double api=gsl_matrix_get(A,p,i);
				double aqi=gsl_matrix_get(A,q,i);
				gsl_matrix_set(A,p,i,c*api-s*aqi);
				gsl_matrix_set(A,q,i,c*aqi+s*api); }
			for(int i=0;i<n;i++){ //We do the similar transformation on the identity to gain the eigenvectors
				double vip=gsl_matrix_get(V,i,p);
				double viq=gsl_matrix_get(V,i,q);
				gsl_matrix_set(V,i,p,c*vip-s*viq);
				gsl_matrix_set(V,i,q,c*viq+s*vip); }
  		}
    }
  } while(changed!=0); //We check if the operation had an effect, if not, we stop.
// We set all non-excact values to 0
for (int ii = eigennr; ii < V->size1; ii++) {
  gsl_vector_set(e,ii,0);
  for(int jj = 0; jj < V->size2; jj++) {
    gsl_matrix_set(V,jj,ii,0);
  }
}
return sweeps;
}


int jacobi_mod_highev(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int eigennr){ // This function is made just as defined proviously.
// Cyclic Jacobi diagonalization, but modified to go from equation to equation starting from the highest, and return exact eigenvalues and matrixes
int changed, sweeps=0, n=A->size1; // Initialised values
for(int i=0;i<n;i++)gsl_vector_set(e,i,gsl_matrix_get(A,i,i)); // We copy the diagonal of A to the eigenmatrix
gsl_matrix_set_identity(V); // We define identity matrix
do{
  changed=0;
  sweeps++; // For each iteration, increase sweeps by one
  int p,q;
	for(p=0;p<eigennr;p++)for(q=p+1;q<n;q++){
		double app=gsl_vector_get(e,p); // The selected value  A_PP is picked out
		double aqq=gsl_vector_get(e,q); // The selected valie A_QQ picked out
		double apq=gsl_matrix_get(A,p,q); // The off diagonal element A_PQ = A_QP is picked out
		double phi=0.5*atan2(2*apq,-(aqq-app)); // The angle to minimize the off diagonal element is found, with a solution to maximise the first element, thereby the fraction is approaching zero from the negative side.
		double c = cos(phi), s = -sin(phi); // The rotation transformation is made, and to rotate the other way, we introduce a minus on the sine, to effectively transpose the rotation
		double app1=c*c*app-2*s*c*apq+s*s*aqq; // The transformation acts on the eigenvalues
		double aqq1=s*s*app+2*s*c*apq+c*c*aqq; // The tramformaiton acts on the eigenvalues
		if(app1!=app || aqq1!=aqq){ changed=1; // If the eigenvalue doesnt change, we have reached the point in which we cannot improve
			gsl_vector_set(e,p,app1); // We set the new eigenvalues into the vectors
			gsl_vector_set(e,q,aqq1);
			gsl_matrix_set(A,p,q,0.0); // We eleminate the off diagonal element in the upper diagonal
			for(int i=0;i<p;i++){ // We act with the transformation on the each row
				double aip=gsl_matrix_get(A,i,p);
				double aiq=gsl_matrix_get(A,i,q);
				gsl_matrix_set(A,i,p,c*aip-s*aiq);
				gsl_matrix_set(A,i,q,c*aiq+s*aip); }
			for(int i=p+1;i<q;i++){ // We act with the transformation on the each row between q and p
				double api=gsl_matrix_get(A,p,i);
				double aiq=gsl_matrix_get(A,i,q);
				gsl_matrix_set(A,p,i,c*api-s*aiq);
				gsl_matrix_set(A,i,q,c*aiq+s*api); }
			for(int i=q+1;i<n;i++){ // We act with the transformation on the each last row
				double api=gsl_matrix_get(A,p,i);
				double aqi=gsl_matrix_get(A,q,i);
				gsl_matrix_set(A,p,i,c*api-s*aqi);
				gsl_matrix_set(A,q,i,c*aqi+s*api); }
			for(int i=0;i<n;i++){ //We do the similar transformation on the identity to gain the eigenvectors
				double vip=gsl_matrix_get(V,i,p);
				double viq=gsl_matrix_get(V,i,q);
				gsl_matrix_set(V,i,p,c*vip-s*viq);
				gsl_matrix_set(V,i,q,c*viq+s*vip); }
  		}
    }
  } while(changed!=0); //We check if the operation had an effect, if not, we stop.
// We set all non-correct values to 0
for (int ii = eigennr; ii < V->size1; ii++) {
  gsl_vector_set(e,ii,0);
  for(int jj = 0; jj < V->size2; jj++) {
    gsl_matrix_set(V,jj,ii,0);
  }
}
return sweeps;
}

int jacobi_versc(gsl_matrix* A, gsl_vector* e, gsl_matrix* V){ // This function is made just as previously, but now with the modified indexing
// Classical Jacobi diagonalization, but modified to go from equation to equation starting from the highest, and return upper triangle of A is destroyed, e and V accumulate eigenvalues and eigenvectors
int changed, sweeps=0, n=A->size1; // Initialised values
int maximum[n];
for(int i=0;i<n;i++)gsl_vector_set(e,i,gsl_matrix_get(A,i,i)); // We copy the diagonal of A to the eigenmatrix
gsl_matrix_set_identity(V); // We define identity matrix
for(int i = 0; i < n; i++){
  maximum[i] = find_max_in_row(A,i);
}

do{
  changed=0;
  sweeps++; // For each iteration, increase sweeps by one
  int p,q;
	for(p=0;p<n-1;p++){
    q=maximum[p]; // We make shure to take only the largest q value in each row
		double app=gsl_vector_get(e,p); // The selected value  A_PP is picked out
		double aqq=gsl_vector_get(e,q); // The selected valie A_QQ picked out
		double apq=gsl_matrix_get(A,p,q); // The off diagonal element A_PQ = A_QP is picked out
		double phi=0.5*atan2(2*apq,(aqq-app)); // The angle to minimize the off diagonal element is found, with a solution to maximise the first element, thereby the fraction is approaching zero from the negative side.
		double c = cos(phi), s = sin(phi); // The rotation transformation is made, and to rotate the other way, we introduce a minus on the sine, to effectively transpose the rotation
		double app1=c*c*app-2*s*c*apq+s*s*aqq; // The transformation acts on the eigenvalues
		double aqq1=s*s*app+2*s*c*apq+c*c*aqq; // The tramformaiton acts on the eigenvalues

    if(app1!=app || aqq1!=aqq){ changed=1; // If the eigenvalue doesnt change, we have reached the point in which we cannot improve
			gsl_vector_set(e,p,app1); // We set the new eigenvalues into the vectors
			gsl_vector_set(e,q,aqq1);
			gsl_matrix_set(A,p,q,0.0); // We eleminate the off diagonal element in the upper diagonal
      maximum[p] = find_max_in_row(A,p);
      maximum[q] = find_max_in_row(A,q);
      for(int i=0;i<p;i++){ // We act with the transformation on the each row
				double aip=gsl_matrix_get(A,i,p);
				double aiq=gsl_matrix_get(A,i,q);
				gsl_matrix_set(A,i,p,c*aip-s*aiq);
				gsl_matrix_set(A,i,q,c*aiq+s*aip); }
			for(int i=p+1;i<q;i++){ // We act with the transformation on the rows between q and p
				double api=gsl_matrix_get(A,p,i);
				double aiq=gsl_matrix_get(A,i,q);
				gsl_matrix_set(A,p,i,c*api-s*aiq);
				gsl_matrix_set(A,i,q,c*aiq+s*api); }
			for(int i=q+1;i<n;i++){ // We act with the transformation on the each last row
				double api=gsl_matrix_get(A,p,i);
				double aqi=gsl_matrix_get(A,q,i);
				gsl_matrix_set(A,p,i,c*api-s*aqi);
				gsl_matrix_set(A,q,i,c*aqi+s*api); }
			for(int i=0;i<n;i++){ //We do the similar transformation on the identity to gain the eigenvectors
				double vip=gsl_matrix_get(V,i,p);
				double viq=gsl_matrix_get(V,i,q);
				gsl_matrix_set(V,i,p,c*vip-s*viq);
				gsl_matrix_set(V,i,q,c*viq+s*vip); }
  		}
    }

  } while(changed!=0); //We check if the operation had an effect, if not, we stop.
return sweeps;
}

int find_max_in_row(gsl_matrix* A, int r){ // We find the maximum in each row by comparing the elements one by one
	int n=A->size1;
  int index_max=0;
	double maximum_val=0; // We have initialised
	for(int i = r+1; i<n; i++){ // For each row
		double a_ri = gsl_matrix_get(A,r,i);
		if(a_ri*a_ri >= maximum_val){ // If we have a higher value, we replace the maximum and set the max index
			maximum_val = a_ri*a_ri;
			index_max = i;
		}
	}
	return index_max;
}
