#include <iostream>
#include <math.h>
#include<stdlib.h>
#include<random>

#include "nuage_b.h"
#include "particule_b.h"
#include "rk4.h"



void systeme(double* q, double t, double* qp, int n) { // Systeme d'ED
qp[0] = q[1];
qp[1] = -q[0];
}


int main() {
int i, n = 2, Nt = 10000;
double t = 0, tfin = 100, dt = (tfin - t) / (Nt - 1);
double* q_euler = (double*)malloc(n * sizeof(double));
double* q_rk4 = (double*)malloc(n * sizeof(double));
 fstream fich("comparaison_euler_rk4.res", ios::out);
q_euler[0] = 1.; q_euler[1] = 0.; // Conditions initiales
q_rk4[0] = 1.; q_rk4[1] = 0.;

for (i = 0; i < Nt; i++) { // Calcul des valeurs de q(dt) a q(tfin)
fich << t  << " " << q_rk4[0] << endl;
rk4(systeme, q_rk4, t, dt, n);
t += dt; }
free(q_rk4);
fich.close();

return 0; }
