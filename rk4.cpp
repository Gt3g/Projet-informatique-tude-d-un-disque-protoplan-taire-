#include <iostream>
#include <math.h>
#include<stdlib.h>
#include<random>

#include "nuage_b.h"
#include "particule_b.h"
#include "rk4.h"

using namespace std;

void rk4(void (*sd)(double*,double,double*,int),
	 double* q, double t, double dt, int n) {
  int i, k, p;
  /* Allocations et initialisations */
  static const int PM = 4;
  static const double c2 = 1./2., c3 = 1./3., c6 = 1./6.;
  static const double a[PM][PM] = {{c2,0,0,0}, {0,c2,0,0}, {0,0,1,0}, {c6,c3,c3,c6}};
  static const double b[PM] = {0,c2,c2,1};
  double** y = (double**)malloc((PM+1) * sizeof(double*)); 
  for(i = 0; i < PM+1; i++)
    y[i] = (double*)malloc(n * sizeof(double));
  double** z = (double**)malloc(PM * sizeof(double*)); 
  for(i = 0; i < PM; i++)
    z[i] = (double*)malloc(n * sizeof(double));
  /* Calcul */
  for(i = 0; i < n; i++)
    y[0][i] = q[i];
  for(p = 1; p <= PM; p++){
    sd(y[p-1], t+b[p-1]*dt, z[p-1], n);
    for(i = 0; i < n; i++)
      y[p][i] = q[i];
    for(k = 0; k < p; k++)
      for(i = 0; i < n; i++)
	y[p][i] = y[p][i] + dt*a[p-1][k]*z[k][i];
  }   
  for(i = 0; i < n; i++)
    q[i] = y[PM][i];
  /* Desallocations */
  for(i = 0; i < PM; i++){
    free(y[i]); free(z[i]);
  }
  free(y[PM]); free(y); free(z);
}
