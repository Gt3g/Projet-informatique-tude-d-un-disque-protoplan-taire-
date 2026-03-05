#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include "nuage.h"
#include "particule.h"

using namespace std;


// ================================================================
//  

// ================================================================
//  Nuage homogène sphérique en rotation → formation de disque
// ================================================================

// ================================================================
//  Affichage console (limité aux 10 premières particules)
// ================================================================
void Nuage::affiche_nuage() {
    int nmax = min((int)vecteur_de_part.size(), 10);
    for (int i = 0; i < nmax; i++) {
        cout << "part" << i << " ";
        vecteur_de_part[i].afficher();
        cout << endl;
    }
    if ((int)vecteur_de_part.size() > nmax)
        cout << "  ... (" << vecteur_de_part.size()-nmax
             << " particules supplémentaires)\n";
}




// ================================================================
//  Nuage homogène sphérique en rotation avec étoile → formation de disque
// ================================================================
void Nuage::init_nuage_homogene_etoile(int    nbr_de_part,
                                       double m_part,
				       double r_part,
                                       double R_nuage,
                                       double omega_0,
                                       double v_therm,
				       double G,
				       double R_etoile,
				       double M_etoile){

  //////////////////////////////////////////////////
  //etoile
  Particule p;
        p.init(0, 0, 0, 0, 0, 0, M_etoile, 0.0);
        vecteur_de_part.push_back(p);
  //////////////////////////////////////////////////////////////////

    mt19937 gen(42);   // graine fixe pour reproductibilité
    uniform_real_distribution<double> dist_u(-1.0, 1.0);
    normal_distribution<double>       gauss(-1.0, 1.0);


    // Distribution uniforme dans la sphère 
    int    placed = 0;
    while (placed < nbr_de_part) {

        double x = dist_u(gen) * R_nuage;
        double y = dist_u(gen) * R_nuage;
        double z = dist_u(gen) * R_nuage;

        if (x*x + y*y + z*z > R_nuage * R_nuage) continue;
	if (x*x + y*y  < R_etoile * R_etoile) continue;

        // ── Vitesse de rotation solide autour de Z ──
        //    v_rot = omega × r_cyl,  direction tangentielle
        //    vx = -omega * y,   vy = +omega * x,   vz = 0
        double vx = -omega_0 * y + v_therm * gauss(gen);
        double vy =  omega_0 * x + v_therm * gauss(gen);
        double vz =               v_therm * gauss(gen);


	// Calcul du moment cinétique //

	double l = 0 ;



        Particule p;
        p.init(x, y, z, vx, vy, vz, m_part, l  );
        vecteur_de_part.push_back(p);
        ++placed;
    }

   
}
