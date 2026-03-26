#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include "nuage.h"
#include "particule.h"

using namespace std;

// ================================================================
//  Affichage console (limité aux 10 premières particules)
// ================================================================
void Nuage::affiche_nuage() const {
    int nmax = min((int)vecteur_de_part.size(), 10);
    for (int i = 0; i < nmax; i++) {
        cout << "part" << i << " ";
        vecteur_de_part[i].afficher();
    }
    if ((int)vecteur_de_part.size() > nmax)
        cout << "  ... (" << vecteur_de_part.size() - nmax
             << " particules supplémentaires)\n";
}

// ================================================================
//  Initialisation : nuage sphérique homogène en rotation
//                  avec étoile centrale
//
//  Algorithme :
//    1. Place l'étoile à l'origine (indice 0)
//    2. Tire des positions aléatoires uniformes dans la sphère
//       (rejection sampling)
//    3. Exclut la zone R < rayon_etoile (déjà occupée par l'étoile)
//    4. Applique la rotation solide : vx = -omega*y, vy = +omega*x
//    5. Ajoute une agitation thermique gaussienne si v_therm > 0
// ================================================================
void Nuage::init_nuage_homogene_etoile(int    nbr_de_part,
                                       double masse_part,
                                       double rayon_part,
                                       double rayon_nuage,
                                       double omega_0,
                                       double v_therm,
                                       double /*G*/,
                                       double rayon_etoile,
                                       double masse_etoile) {

    // ── Étoile centrale (indice 0) ────────────────────────────────
    vecteur_de_part.push_back(Particule(0, 0, 0, 0, 0, 0, masse_etoile, 0.0));

    // ── Générateur aléatoire (graine fixe pour reproductibilité) ──
    mt19937 gen(42);
    uniform_real_distribution<double> dist_u(-1.0, 1.0);
    normal_distribution<double>       gauss(0.0, 1.0);

    // ── Distribution uniforme dans la sphère (rejection sampling) ─
    int placed = 0;
    while (placed < nbr_de_part) {

        double x = dist_u(gen) * rayon_nuage;
        double y = dist_u(gen) * rayon_nuage;
        double z = dist_u(gen) * rayon_nuage;

        // Rejet si hors de la sphère
        if (x*x + y*y + z*z > rayon_nuage * rayon_nuage) continue;

        // Rejet si trop proche de l'étoile (zone d'exclusion)
        if (x*x + y*y < rayon_etoile * rayon_etoile)     continue;

        // ── Vitesse de rotation solide autour de Z ────────────────
        //    v_rot = omega × r_cyl (direction tangentielle)
        //    vx = -omega * y,   vy = +omega * x,   vz = 0
        double vx = -omega_0 * y + v_therm * gauss(gen);
        double vy =  omega_0 * x + v_therm * gauss(gen);
        double vz =               v_therm * gauss(gen);

        vecteur_de_part.push_back(Particule(x, y, z, vx, vy, vz, masse_part, 0.0));
        ++placed;
    }
}
