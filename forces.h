
#ifndef FORCES_H
#define FORCES_H

#include <vector>
#include "nuage.h"
#include "nuage_solver.h"

// ================================================================
//  Calcul des forces exercées sur les particules
// ================================================================

// ── Gravitation de la particule j_src sur la particule i ─────────
// Loi de Newton avec softening :
//   a = G * m_j * dr / (|dr|² + eps²)^(3/2)
//
// Les accélérations ax, ay, az sont incrémentées (+=).
void gravite_softened(int i, int j_src,
                      const std::vector<double>& s,
                      const Nuage& nuage,
                      const SimParams& params,
                      double& ax, double& ay, double& az);

// ── Amortissement vertical dans le disque ────────────────────────
// Modèle disque mince isotherme irradié :
//   T(R) = T0 * sqrt(R0/R)
//   H(R) = cs(R) / Omega_K(R)   ∝ R^(5/4)
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
//    La vitesse d'une particule est décomposée en trois composantes
//    dans le repère cylindrique local (r̂_cyl, φ̂_cyl, ẑ) :
//
//      v = v_r * r̂_cyl  +  v_φ * φ̂_cyl  +  v_z * ẑ
//
//    • v_r > 0 signifie une orbite excentrique (mouvement radial)
//    • v_φ représente la rotation keplerienne
//    • v_z est la composante verticale (hors plan)
//
//  Force d'amortissement de l'excentricité :
//
//      F_e = -(1/t_e) * v_r * r̂_cyl
//      avec   1/t_e = f_circ * Omega_K(R_cyl)
//
//  Cette force :
//    - amortit v_r → 0 : l'orbite devient circulaire
//    - n'affecte pas v_φ : conserve le moment cinétique L_z
//    - n'affecte pas v_z : la dissipation verticale est déléguée
//      à amortissement_vertical()
//
//  Le rayon cylindrique R_cyl est utilisé (et non le rayon 3D)
//  car la dynamique keplerienne du disque est dans le plan x-y.
//
//  Référence : Papaloizou & Larwood (2000), Tanaka & Ward (2004)
void circularisation(int i,
                     const std::vector<double>& s,
                     const Nuage& nuage,
                     const SimParams& params,
                     double masse_etoile,
                     double& ax, double& ay, double& az);

#endif
