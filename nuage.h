#ifndef NUAGE_H
#define NUAGE_H

#include <vector>
#include "particule.h"

using namespace std;

class Nuage {
 public:
  vector<Particule> vecteur_de_part;




  //=============================================
  //nuage homo avec étoile
  //==========================================

  void init_nuage_homogene_etoile(int    nbr_de_part,
                                  double M_part,
				  double r_part,
                                  double R_nuage,
                                  double omega_0,
                                  double v_therm,
				  double G,
				  double R_etoile,
				  double M_etoile);
				      

  

  void affiche_nuage();
};

#endif
