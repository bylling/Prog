#include<assert.h>
#include<math.h>
double mylinterp(int n, double * x, double * y, double z){
  assert(n>1 && z>=x[0] && z<=x[n-1]);
  int i=0, j=n-1; /* We initialise values for use in binary search*/
  while(j-i>1){
    int m=floor((i+j)/2);
    if(z>x[m]){
      i = m;    /* We set the lower boarder to the middle*/
    }
    else{
      j = m;     /* We set the upper boarder to the middle*/
    }
  }
double solution = y[i]+(y[i + 1]- y[i])/(x[i+1] - x[i])*(z-x[i]); /* Using lin. interpolation, we find the value within the i and i+1 to return*/
return solution;


}
