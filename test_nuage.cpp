#include <iostream>
#include <math.h>
#include<stdlib.h>
#include<random>

#include "nuage_b.h"
#include "particule_b.h"

using namespace std;

int main(){
  int n = 10;
  int L = 1000;

  Nuage nuage_test;

  nuage_test.init_random(n , 10 , 1000 ,1, 0,0 ,L);

  nuage_test.affiche_nuage();
}
