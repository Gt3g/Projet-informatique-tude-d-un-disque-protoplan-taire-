#include <iostream>
#include <cmath>
#include "params.h"
#include "nuage.h"
#include "nuage_solver.h"

// ================================================================
//  Formation d'un disque protoplanétaire
//  à partir d'un nuage moléculaire sphérique homogène
//
//  Physique :
//    1. Nuage sphérique uniforme en rotation solide autour de Z
//       avec beta = E_rot / E_grav = 0.04 (valeur observationnelle)
//    2. L'auto-gravité effondre le nuage (temps ~ t_ff)
//       + auto-gravité grain↔grain si ENABLE_GRAIN_GRAIN = true
//    3. Le moment cinétique est conservé → aplatissement en disque
//    4. Les dissipations (amortissement vertical, circularisation)
//       sont activées à t = T_DISS_START * t_ff
//
//  Paramètres modifiables dans params.h.
//  Validation : vérifier conservation de L_z et E_tot
//               dans diagnostics.fich (sans dissipation).
// ================================================================

int main() {

    // ── Calcul de la vitesse angulaire initiale ───────────────────
    // omega = sqrt(beta * 3*G*M / R³)   [rotation solide uniforme]
    const double omega_0 = std::sqrt(BETA_ROT * 3.0 * G_PHYS * MASSE_ETOILE
                                     / (RAYON_NUAGE * RAYON_NUAGE * RAYON_NUAGE));

    // ── Calcul du temps de chute libre ───────────────────────────
    // t_ff = sqrt(3π / (32 G ρ₀))
    const double rho_0 = MASSE_ETOILE
                       / (4.0/3.0 * M_PI * RAYON_NUAGE * RAYON_NUAGE * RAYON_NUAGE);
    const double t_ff  = std::sqrt(3.0 * M_PI / (32.0 * G_PHYS * rho_0));

    // ── Paramètres temporels ─────────────────────────────────────
    const double t0   = 0.0;
    const double tEnd = DUREE_EN_TFF * t_ff;
    const double h    = t_ff / PAS_EN_TFF;

    // ── Softening gravitationnel ──────────────────────────────────
    const double eps = EPS_FACTOR * RAYON_NUAGE / std::cbrt((double)N_PARTICULES);

    // ── Affichage des paramètres ──────────────────────────────────
    std::cout << "=== Simulation formation de disque protoplanétaire ===\n\n";
    std::cout << "Paramètres physiques :\n";
    std::cout << "  Masse étoile      : " << MASSE_ETOILE / M_SOLEIL << " M_sun\n";
    std::cout << "  Masse disque      : " << MASSE_DISK / M_SOLEIL   << " M_sun"
              << "  (" << FRACTION_DISK*100 << "% M_etoile)\n";
    std::cout << "  Rayon nuage       : " << RAYON_NUAGE / UA << " UA\n";
    std::cout << "  beta (rot/grav)   : " << BETA_ROT << "\n";
    std::cout << "  Temps de chute libre t_ff = "
              << t_ff / AN << " ans\n";
    std::cout << "  Auto-gravité g-g  : " << (ENABLE_GRAIN_GRAIN ? "OUI" : "NON") << "\n";
    std::cout << "  Dissipation dès   : t = " << T_DISS_START << " t_ff = "
              << T_DISS_START * t_ff / AN << " ans\n\n";
    std::cout << "Paramètres numériques :\n";
    std::cout << "  N particules      : " << N_PARTICULES << "\n";
    std::cout << "  Softening eps     : " << eps / UA << " UA\n";
    std::cout << "  Simulation de     : " << t0/AN << " à " << tEnd/AN << " ans\n";
    std::cout << "  Pas de temps h    : " << h/AN << " ans\n\n";

    // ── Initialisation du nuage ───────────────────────────────────
    Nuage nuage;
    nuage.init_nuage_homogene_etoile(N_PARTICULES,
                                     MASSE_PART,
                                     RAYON_PART,
                                     RAYON_NUAGE,
                                     omega_0,
                                     V_THERM,
                                     G_PHYS,
                                     RAYON_ETOILE,
                                     MASSE_ETOILE);

    std::cout << "Particules créées : "
              << nuage.vecteur_de_part.size() << "\n\n";

    // ── Paramètres du solveur ─────────────────────────────────────
    SimParams params;
    params.G              = G_PHYS;
    params.eps            = eps;
    params.omega          = 0.0;         // référentiel inertiel
    params.a_z            = 0.0;
    params.pas_ecriture   = PAS_ECRITURE;
    params.pas_diag       = PAS_DIAG;
    params.enable_grain_grain = ENABLE_GRAIN_GRAIN;

    // Activation différée des dissipations en secondes
    params.t_diss_start   = T_DISS_START * t_ff;

    // ── Simulation ───────────────────────────────────────────────
    rk4_nuage(nuage, params, t0, tEnd, h, FICHIER_SORTIE, FICHIER_DIAG);

    return 0;
}
