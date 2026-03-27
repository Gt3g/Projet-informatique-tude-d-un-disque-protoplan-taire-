#ifndef NUAGE_H
#define NUAGE_H

#include <vector>
#include "particule.h"

// ================================================================
//  Classe Nuage
//  Représente l'ensemble des particules du système :
//    - vecteur_de_part[0]    : étoile centrale (masse fixe, position fixe)
//    - vecteur_de_part[1..N] : super-grains du nuage moléculaire
// ================================================================
class Nuage {
public:
    std::vector<Particule> vecteur_de_part;

    // Constructeur par défaut
    Nuage() = default;

    // -- Initialisation ------------------------------------------
    // Remplit vecteur_de_part avec nbr_de_part grains + 1 étoile centrale.
    // Distribution uniforme dans une sphère creuse de rayon externe rayon_nuage et de rayon interne rayon étoile,
    // en rotation solide autour de Z (vitesse angulaire omega_0),
    // avec agitation thermique gaussienne v_therm.
    void init_nuage_homogene_etoile(int    nbr_de_part,
                                    double masse_part,
                                    double rayon_part,
                                    double rayon_nuage,
                                    double omega_0,
                                    double v_therm,
                                    double G,
                                    double rayon_etoile,
                                    double masse_etoile);

    // -- Affichage console (10 premières particules) -------------
    void affiche_nuage() const;
};

#endif
