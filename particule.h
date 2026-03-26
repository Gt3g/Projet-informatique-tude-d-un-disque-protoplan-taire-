#ifndef PARTICULE_H
#define PARTICULE_H

// ================================================================
//  Classe Particule
//  Représente une super-particule du nuage moléculaire ou l'étoile.
//
//  Convention :
//    - Unités SI : position en m, vitesse en m/s, masse en kg
//    - vecteur_de_part[0] désigne toujours l'étoile centrale
// ================================================================
class Particule {
public:
    double coord_x, coord_y, coord_z;  // position (m)
    double vit_x,   vit_y,   vit_z;   // vitesse (m/s)
    double masse;                       // masse (kg)
    double moment_cin;                  // moment cinétique scalaire (kg·m²/s)

    // ── Constructeur par défaut : particule à l'origine, au repos ─
    Particule();

    // ── Constructeur principal ────────────────────────────────────
    Particule(double x, double y, double z,
              double vx, double vy, double vz,
              double m, double l = 0.0);

    // ── Affichage console (position, vitesse, masse) ──────────────
    void afficher() const;

    // ── Distance à l'origine (m) ──────────────────────────────────
    double dist_origine() const;

    // ── Carré de la distance à l'axe Z (m²) ──────────────────────
    double dist_axe_carre() const;
};

#endif
