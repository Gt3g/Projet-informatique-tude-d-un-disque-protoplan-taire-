#ifndef NUAGE_SOLVER_H
#define NUAGE_SOLVER_H

#include "nuage.h"
#include <string>

// ============================================================
//  Solveur RK4 pour un nuage de particules
//
//  Équations résolues pour chaque particule i :
//
//    dx_i/dt  = vx_i
//    dy_i/dt  = vy_i
//    dz_i/dt  = vz_i
//
//    dvx_i/dt = (1/m_i) * Σ_{j≠i} F_grav_x(i,j)
//    dvy_i/dt = (1/m_i) * Σ_{j≠i} F_grav_y(i,j)
//    dvz_i/dt = (1/m_i) * Σ_{j≠i} F_grav_z(i,j)  +  a_z
//
//  Gravitation newtonienne :
//    F_grav(i,j) = G * m_i * m_j / (r_ij² + eps²)  *  r̂_ij
//    (eps = softening pour éviter les divergences à r→0)
//
//  Force d'inertie sur z :
//    a_z = -omega² * z_i   (pseudo-force centrifuge selon z)
//    ou   a_z = cste        (accélération uniforme selon z)
//
// ============================================================

// Constante gravitationnelle (unités SI, adapter si besoin)
constexpr double G_CONST   = 6.674e-11;

// Paramètre de softening (évite r=0)
constexpr double EPS_SOFT  = 1e-3;

// ============================================================
//  Paramètres de simulation
// ============================================================
struct SimParams {
    double G       = G_CONST;   // constante gravitationnelle
    double eps     = EPS_SOFT;  // softening gravitationnel
    double omega   = 0.0;       // pulsation pour force d'inertie centrifuge
    double a_z     = 0.0;       // accélération uniforme selon z (en plus)
    // Force d'inertie totale selon z :  acc_z = -omega²*z + a_z
};

// ============================================================
//  Solveur RK4
//
//  nuage    : état initial (modifié en place à chaque pas)
//  params   : paramètres physiques
//  t0       : temps initial
//  tEnd     : temps final
//  h        : pas de temps
//  filename : fichier de sortie (.fich)
//
//  Format du fichier :
//  t   i   x   y   z   vx   vy   vz
// ============================================================
void rk4_nuage(Nuage&             nuage,
               const SimParams&   params,
               double             t0,
               double             tEnd,
               double             h,
               const std::string& filename);

#endif // NUAGE_SOLVER_H
