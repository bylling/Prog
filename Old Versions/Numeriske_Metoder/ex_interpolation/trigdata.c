
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main() {
  double start = 0;
  double stop = 3 * M_PI;
  double step = 0.5;
  for (double x = start; x < stop; x=x+step) {
    printf("%g  %g  %g\n",x,cos(x),sin(x));
  }
return 0;
}
