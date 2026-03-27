#ifndef FORCES_H
#define FORCES_H

#include <vector>
#include "nuage.h"
#include "nuage_solver.h"

// ================================================================
//  Calcul des forces exercées sur les particules
// ================================================================

// ── Gravitation de la particule j sur la particule i ─────────────
// Loi de Newton avec softening de Plummer :
//   a = G * m_j * dr / (|dr|² + eps²)^(3/2)
//
void gravite_softened(int i, int j,
                      const std::vector<double>& s,
                      const Nuage& nuage,
                      const SimParams& params,
                      double& ax, double& ay, double& az);

// ── Amortissement vertical dans le disque ────────────────────────
// Modèle disque mince isotherme irradié :
//   T(R) = T0 * sqrt(R0/R)
//   H(R) = cs(R) / Omega_K(R)   ~ R^(5/4)
//
// Retourne a_z,diss = -F_DISS_DISK * Omega_K * vz
// si la particule est dans le disque (|z| < H), 0 sinon.
double amortissement_vertical(double xi, double yi, double zi,
                               double vzi,
                               double masse_etoile,
                               const SimParams& params);

// ── Circularisation des orbites ───────────────────────────────────
//
//  Modèle physique :
//    La vitesse d'une particule est décomposée en coordonnées cylindriques :
//      v = v_r * r̂_cyl  +  v_φ * φ̂_cyl  +  v_z * ẑ
//
//    • v_r > 0 signifie une orbite excentrique (mouvement radial)
//    • v_φ représente la rotation képlérienne
//
//  Force d'amortissement de l'excentricité :
//      F_e = -(1/t_e) * v_r * r̂_cyl
//      avec   1/t_e = f_circ * Omega_K(R_cyl)
//
//  Cette force :
//    - amortit v_r → 0 : l'orbite devient circulaire
//    - n'affecte pas v_φ : conserve le moment cinétique L_z
//    - n'affecte pas v_z : délégué à amortissement_vertical()
//
//  Modifie uniquement ax et ay (pas az).
void circularisation(int i,
                     const std::vector<double>& s,
                     const SimParams& params,
                     double masse_etoile,
                     double& ax, double& ay);

// ── Pression par méthode Particle-Mesh (PM) ───────────────────────
//
//  Calcule en une seule fois les accélérations de pression pour
//  toutes les particules (sauf l'étoile i=0).
//  Complexité : O(N + Ng³) au lieu de O(N²).
//
//  Paramètres utilisés depuis SimParams :
//    params.n_grid   : résolution Ng de la grille cubique
//    params.cs_pm    : vitesse du son (m/s)
//    params.gamma_pm : indice polytropique
void pression_grille(const std::vector<double>& s,
                     const Nuage& nuage,
                     const SimParams& params,
                     std::vector<double>& ax_press,
                     std::vector<double>& ay_press,
                     std::vector<double>& az_press);

#endif
