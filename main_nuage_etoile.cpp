#include <iostream>
#include <cmath>
#include "nuage.h"
#include "nuage_solver.h"

// ================================================================
//  Formation d'un disque protoplanétaire
//  à partir d'un nuage moléculaire sphérique homogène
//
//  Physique :
//    1. Nuage sphérique uniforme en rotation solide autour de Z
//    2. L'auto-gravité effondre le nuage (temps ~ t_ff)
//    3. Le moment cinétique est conservé
//    4. L'effondrement selon Z est libre → aplatissement
//    5. Dans le plan XY, la rotation freine → un disque se forme
//
//  Paramètres réalistes (cœur de nuage moléculaire) :
//    M_total  = 1 M_soleil = 1.989e30 kg
//    R_cloud  = 1000 UA    = 1.496e14 m
//    beta     = E_rot/|E_grav| = 0.05  (rotation modérée)
//    v_therm  ~ 200 m/s   (son à ~10 K dans H2)
//    t_ff     ~ 5600 ans
//    Simulation : 3 t_ff ~ 17 000 ans
//    Pas h    : t_ff / 200 ~ 28 ans
//    Frames   : 600 pas / 2 pas_ecr = 300 frames
// ================================================================

int main() {

    // ── Constantes ───────────────────────────────────────────
    const double G  = 6.674e-11;
    const double UA = 1.496e11;    // 1 Unité Astronomique (m)
    const double AN = 3.156e7;     // 1 an (s)

    // ── Paramètres du nuage ──────────────────────────────────
    const int    N       = 2000;           // super-grains
    const double M_part = 1;      // 1 M_soleil/100
    const double r_part = 1;
    const double R_nuage = 100.0 * UA;  // rayon initial
    const double M_etoile = 1.989e30;    // 1_M_soleil
    const double R_etoile = 40 * UA;
    

    // Paramètre de rotation beta = E_rot / |E_grav| = 1
    // => omega = sqrt( beta *G*M / R³/3 )
    const double beta    = 4.0;
    const double omega_0 = std::sqrt(beta  * G * M_etoile
                                     / (R_nuage * R_nuage * R_nuage)/3.0);

    // Agitation thermique : vitesse du son dans H2 à ~10 K
    const double v_therm = 0.0;  // m/s

    // ── Calcul du temps de chute libre ───────────────────────
    const double rho_0 = M_etoile
                       / (4.0/3.0 * M_PI * R_nuage * R_nuage * R_nuage);
    const double t_ff  = std::sqrt(3.0 * M_PI / (32.0 * G * rho_0));

    std::cout << "Temps de chute libre t_ff = "
              << t_ff << " s = " << t_ff/AN << " ans\n";

    // ── Paramètres temporels ─────────────────────────────────
    const double t0      = 0.0;
    const double tEnd    = 15.0 * t_ff;    // 3 t_ff ~ 17 000 ans
    const double h       = t_ff / 100.0;  // ~ 28 ans par pas
    const int    pas_ecr = 2;             // 1 frame/2 pas → ~300 frames

    // Softening : fraction du rayon initial moyen inter-particules
    // r_moy ~ R_cloud / N^(1/3) ~ 1000 UA / 8 = 125 UA
    const double eps = 0.3 * R_nuage / std::cbrt((double)N);

    // ── Initialisation ───────────────────────────────────────
    Nuage nuage;
    nuage.init_nuage_homogene_etoile(N,
				     M_part,
				     r_part,
				     R_nuage,
				     omega_0,
				     v_therm,
				     G,
				     R_etoile,
				     M_etoile);

    std::cout << "Particules créées : "
              << nuage.vecteur_de_part.size() << "\n\n";

    // ── Paramètres physiques ─────────────────────────────────
    SimParams params;
    params.G            = G;
    params.eps          = eps;
    params.omega        = 0.0;   // référentiel inertiel
    params.a_z          = 0.0;
    params.pas_ecriture = pas_ecr;

    std::cout << "Softening eps = " << eps << " m = "
              << eps/UA << " UA\n";
    std::cout << "Simulation de " << t0/AN << " à "
              << tEnd/AN << " ans, h = " << h/AN << " ans\n\n";

    // ── Simulation ───────────────────────────────────────────
    rk4_nuage(nuage, params, t0, tEnd, h, "evolution_nuage.fich");

    return 0;
}
