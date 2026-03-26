#ifndef NUAGE_SOLVER_H
#define NUAGE_SOLVER_H

#include "nuage.h"
#include "params.h"
#include <string>

// ================================================================
//  Paramètres physiques et numériques passés au solveur RK4
// ================================================================
struct SimParams {
    double G            = G_PHYS;   // Constante gravitationnelle (m³ kg⁻¹ s⁻²)
    double eps          = 1e-3;     // Longueur de softening gravitationnel (m)
    double omega        = 0.0;      // Vitesse angulaire du référentiel (rad/s)
    double a_z          = 0.0;      // Accélération externe selon Z (m/s²)
    int    pas_ecriture = 1;        // Fréquence d'écriture des positions
    int    pas_diag     = 1;        // Fréquence d'écriture des diagnostics

    // Taux d'amortissement de l'excentricité (sans dimension).
    // Le taux effectif vaut f_circ * Omega_K(R) → unité s⁻¹.
    // Voir params.h pour l'interprétation physique.
    double f_circ = F_CIRC;

    // ── Auto-gravité grain ↔ grain ────────────────────────────────
    // Si true : boucle O(N²) sur toutes les paires (i,j), i≠j.
    // Si false : seulement grain ← étoile (test-particle limit).
    bool enable_grain_grain = ENABLE_GRAIN_GRAIN;

    // ── Activation différée des dissipations ──────────────────────
    // Les forces de dissipation sont désactivées pour t < t_diss_start.
    // t_diss_start est exprimé en secondes (converti depuis t_ff dans main).
    // Valeur 0.0 → dissipations actives dès le début.
    double t_diss_start = 0.0;   // (s) — initialisé dans main après calcul de t_ff
};

// ================================================================
//  Solveur RK4 principal
//
//  Intègre le système de N particules de t0 à tEnd avec le pas h.
//  Écrit les positions/vitesses dans filename_traj.
//  Écrit les diagnostics physiques dans filename_diag.
//
//  Diagnostics calculés à chaque pas_diag pas :
//    - E_cin  : énergie cinétique totale
//    - E_grav : énergie potentielle gravitationnelle (étoile + grain-grain)
//    - E_tot  : E_cin + E_grav
//    - L_z    : moment cinétique selon Z
//    - dLz    : variation relative de L_z depuis t=0
// ================================================================
void rk4_nuage(Nuage&             nuage,
               const SimParams&   params,
               double             t0,
               double             tEnd,
               double             h,
               const std::string& filename_traj,
               const std::string& filename_diag);

#endif
