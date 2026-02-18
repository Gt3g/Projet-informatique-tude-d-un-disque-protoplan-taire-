#ifndef NUAGE_H
#define NUAGE_H

#include <vector>
#include "particule_b.h"

using namespace std;

class Nuage{
 public:

  vector<Particule> vecteur_de_part;

  void init_random(int, double, double , double , double , double ,double );

  void affiche_nuage();
};


#endif
