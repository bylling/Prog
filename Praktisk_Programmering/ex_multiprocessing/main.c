#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <pthread.h>


struct harm {int a,b; double result;};

void* bar(void* param){

struct harm* data = (struct harm*)param;
	int a=data->a;
  int b=data->b;
  double x,y,z;
  int count=0;
  double pi;

  for(int i = 0; i<a; i++) {
    unsigned int seed = time(NULL)+i*b;
    x=rand_r(&(seed))/(double)RAND_MAX;
    y=rand_r(&(seed))/(double)RAND_MAX;
    z=x*x + y*y;
    if (z<=1) {
      count++;
    }
    }
    pi=(double)count/a*4;
	data->result=pi;
  return 0;
}

int main(int argc, char** argv){
	int n=1e6;
	int mid=n/2;
	struct harm data1,data2;
	data1.a=mid;
  data2.a=mid;
  data1.b=1;
  data2.b=2;
	pthread_t traad1, traad2;
	pthread_create(&traad1,NULL,bar,(void*)&data1);
  pthread_create(&traad2,NULL,bar,(void*)&data2);

	pthread_join(traad1,NULL);
  pthread_join(traad2,NULL);

  printf("Monte Carlo calculation for first processor using %d datapoints = %g\n",n/2,data1.result);
  printf("Monte Carlo calculation for second processor using %d datapoints = %g\n",n/2,data2.result);

double resulttotal=(data1.result+data2.result)/2;

	printf("Monte Carlo calculation using %d datapoints = %g\n",n,resulttotal);

/*
	bar((void*)&data);
	printf("sum from %i to %i = %g\n",1,n,data.result);
*/
return EXIT_SUCCESS;
}
