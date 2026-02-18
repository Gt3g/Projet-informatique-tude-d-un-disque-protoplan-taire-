#ifndef RK4_H
#define RK4_H


#include <iostream>
#include <math.h>
#include<stdlib.h>
#include<random>

#include "nuage_b.h"
#include "particule_b.h"

using namespace std;

void rk4(void (*sd)(double*,double,double*,int),
	 double* , double , double , int );


#endif
