#include <stdlib.h>
int result;
int equal(double a, double b, double tau, double epsilon){
double diff=abs(a-b);
double reldiff=diff/(abs(a)+abs(b));
if (diff < tau) {result = 1;}
else if (reldiff < epsilon/2) {result =1;}
else {result = 0;};
return result;
}
