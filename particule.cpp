#include <iostream>
#include <cmath>
#include "particule.h"

using namespace std;

// ── Constructeur par défaut ──────────────────────────────────────
Particule::Particule()
    : coord_x(0.0), coord_y(0.0), coord_z(0.0),
      vit_x(0.0),   vit_y(0.0),   vit_z(0.0),
      masse(1.0),   moment_cin(0.0)
{}

// ── Constructeur principal ───────────────────────────────────────
Particule::Particule(double x, double y, double z,
                     double vx, double vy, double vz,
                     double m, double l)
    : coord_x(x), coord_y(y), coord_z(z),
      vit_x(vx),  vit_y(vy),  vit_z(vz),
      masse(m),   moment_cin(l)
{}

// ── Affichage console ────────────────────────────────────────────
void Particule::afficher() const {
    cout << " x="   << coord_x << " y="   << coord_y << " z="   << coord_z   << endl;
    cout << " v_x=" << vit_x   << " v_y=" << vit_y   << " v_z=" << vit_z     << endl;
    cout << " m="   << masse   << " l="   << moment_cin                       << endl;
}

// ── Distance à l'origine ─────────────────────────────────────────
double Particule::dist_origine() const {
    return sqrt(coord_x*coord_x + coord_y*coord_y + coord_z*coord_z);
}

// ── Carré de la distance à l'axe Z ──────────────────────────────
double Particule::dist_axe_carre() const {
    return coord_x*coord_x + coord_y*coord_y;
}
