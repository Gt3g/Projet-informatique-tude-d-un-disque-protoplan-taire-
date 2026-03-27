#ifndef PARAMS_H
#define PARAMS_H

// ================================================================
//  Paramètres physiques et numériques de la simulation
//  Formation d'un disque protoplanétaire autour d'une protoétoile
//
//  Unités SI partout (m, kg, s).
//  Ce fichier centralise tous les paramètres modifiables.
//  Après toute modification : recompiler avec make
// ================================================================


// -- Constantes physiques fondamentales --------------------------
constexpr double G_PHYS   = 6.674e-11;  // Constante gravitationnelle (m³ /kg /s²)
constexpr double UA       = 1.496e11;   // Unité Astronomique (m)
constexpr double AN       = 3.156e7;    // Année (s)
constexpr double M_SOLEIL = 1.989e30;   // Masse solaire (kg)


// -- Paramètres de l'étoile centrale -----------------------------
constexpr double MASSE_ETOILE  = M_SOLEIL;      // Masse de l'étoile (kg)
constexpr double RAYON_ETOILE  = 60.0  * UA;    // Rayon d'exclusion (m)


// -- Paramètres du nuage moléculaire -----------------------------
constexpr int    N_PARTICULES  = 1000;
// Rayon typique d'un cœur dense protostellaire : 50–200 UA
constexpr double RAYON_NUAGE   = 150.0 * UA;


// Masse totale du disque / nuage :
//   fraction disque typique f_disk = 0.05–0.20 de la masse stellaire
//   On prend f_disk = 0.10 -> M_disk = 0.10 M_etoile
constexpr double FRACTION_DISK = 0.1;
constexpr double MASSE_DISK    = FRACTION_DISK * MASSE_ETOILE;   // ~2e29 kg
constexpr double MASSE_PART    = MASSE_DISK / N_PARTICULES;       // ~4e26 kg/grain

constexpr double RAYON_PART    = 0.001;  // rayon numérique 


// -- Paramètre de rotation ---------------------------------------
// beta = E_rot / |E_grav|
//
//  E_rot = 40% E_grav -> nuage en quasi-équilibre centrifuge
constexpr double BETA_ROT = 0.4;


// -- Agitation thermique -----------------------------------------
// Vitesse du son dans H2 à ~10 K : cs ~ 200 m/s
constexpr double V_THERM = 0.0;  // m/s (0 = agitation désactivée)


// ── Paramètres temporels ─────────────────────────────────────────
constexpr double DUREE_EN_TFF = 15.0;   // Durée totale (en unités de t_ff)
constexpr double PAS_EN_TFF   = 200.0;  // Nombre de pas par t_ff
constexpr int    PAS_ECRITURE = 5;      // Écriture positions toutes les N pas
constexpr int    PAS_DIAG     = 1;      // Écriture diagnostics tous les N pas


// ── Softening gravitationnel ─────────────────────────────────────
// eps = EPS_FACTOR * R_nuage / N^(1/3)
constexpr double EPS_FACTOR = 0.25;


// ── Auto-gravité grain ↔ grain ───────────────────────────────────
// true  : interactions N-corps complètes O(N**2) — physiquement correct
// false : test-particle (grain ← étoile seulement) — rapide

constexpr bool ENABLE_GRAIN_GRAIN = false;


// ── Activation différée des dissipations ─────────────────────────
// Les forces de dissipation (circularisation, amortissement vertical)
// ne sont physiques qu'une fois un disque formé (t > quelques t_ff).
// Les activer pendant l'effondrement sphérique (t < t_ff) biaiserait
// artificiellement la dynamique.
// (l'effondrement vers un disque prend typiquement 1–3 t_ff)
constexpr double T_DISS_START = 1.0;   // en unités de t_ff


// ── Paramètres du disque : amortissement vertical ────────────────
// Modèle disque mince irradié passif :
//   T(R) = T0 * (R0/R)^(1/2)  
//   H(R) = cs(R) / Omega_K(R) ~  R^(5/4)  (disque évasé)
constexpr double MU_MOL_DISK = 2.0e-3;  // Masse molaire H2 (kg/mol)
constexpr double R_GAZ_DISK  = 8.314;   // Constante des gaz parfaits (J/mol/K)
constexpr double T0_DISK     = 150.0;   // Température à R = 1 UA (K)
                                     
constexpr double R0_DISK     = UA;

// Coefficient d'amortissement vertical :
// f_diss ~ 1–3 
constexpr double F_DISS_DISK = 1.5;


// ── Paramètres du disque : circularisation des orbites ──────────
//
//  F_e = -(f_circ * Omega_K) * v_r * r̂_cyl
//  t_e = 1 / (f_circ * Omega_K)  ≈  (1/f_circ) / (2π)  périodes
//
//  f_circ = 0.01 → ~16 périodes pour circulariser une orbite
constexpr double F_CIRC = 0.05;


// ── Pression par méthode Particle-Mesh (PM) ──────────────────────
//
//  Modélise la pression du gaz sans boucle O(N²).
//  Complexité : O(N + Ng³)  →  bien plus rapide que SPH pour N≥200.
//
//  Algorithme (3 étapes par pas de temps) :
//
//  1. DÉPÔT (particles → grille) :
//       Chaque particule distribue sa masse sur les 8 cellules voisines
//       par interpolation trilinéaire (schéma Cloud-In-Cell, CIC).
//       → grille de densité ρ[ix][iy][iz]
//
//  2. GRILLE (calculs sur Ng³ cellules) :
//       P[ix][iy][iz] = K * ρ^γ          (loi d'état polytropique)
//       dP calculé par différences finies centrées :
//         (∂P/∂x)[i] = (P[i+1] - P[i-1]) / (2 dx)
//
//  3. INTERPOLATION (grille → particles) :
//       Mêmes poids CIC pour lire le gradient sur la grille.
//       a_press,i = -(1/ρ_i) * dP(r_i)
//
//  Loi d'état : P = K * ρ^γ
//    γ = 1 (isotherme) → P = cs² * ρ  →  K = cs²
//    γ = 5/3 (adiabatique)
//
//  Résolution de la grille Ng :
//    Ng > 16 -> système instable, amélioration future
//
//  Vitesse du son cs_pm :
//    À 1 UA dans H2 à ~150 K : cs ≈ 700 m/s.
//    Augmenter cs_pm → pression plus forte → disque plus épais.
constexpr bool   ENABLE_PRESSURE = true;
constexpr int    N_GRID          = 16;      // résolution de la grille (Ng³ cellules)
constexpr double CS_PM           = 100.0;   // vitesse du son isotherme (m/s)
constexpr double GAMMA_PM        = 1.0;     // indice polytropique (1 = isotherme)


// ── Domaine cylindrique de la grille PM ──────────────────────────
//
//  La pression est négligeable en dehors du disque : la grille est
//  donc restreinte au cylindre |R_cyl| ≤ R_DISK_PM, |z| ≤ H_DISK_PM.
//  Les particules hors de ce cylindre ne déposent pas de masse et
//  ne reçoivent pas de force de pression.
//
//  Effets de ce choix :
//    • Résolution effective bien meilleure : la cellule couvre
//      (2*R_disk / Ng) × (2*H_disk / Ng) au lieu de (2*R_nuage / Ng)
//    • Calcul plus rapide (dépôt/interpolation sur moins de particules)
//
//  Valeurs recommandées :
//    R_DISK_PM = 0.5–0.8 * RAYON_NUAGE  (rayon du disque en formation)
//    H_DISK_PM = 0.1–0.2 * RAYON_NUAGE  (hauteur caractéristique du disque)
//      → H/R ~ 0.1–0.2 est typique d'un disque protoplanétaire mince
//
//  Si une particule sort du cylindre en cours de simulation (orbite
//  très excentrique, éjection), elle ne reçoit simplement pas de force
//  de pression
constexpr double R_DISK_PM = 0.75 * RAYON_NUAGE;  // rayon max du cylindre (m)
constexpr double H_DISK_PM = 0.15 * RAYON_NUAGE; // demi-hauteur du cylindre (m)

// ── Fichiers de sortie ────────────────────────────────────────────
constexpr const char* FICHIER_SORTIE = "evolution_nuage.fich";
// Fichier de diagnostics : t, E_cin, E_grav, E_tot, L_z, delta_Lz
// Permet de vérifier la conservation du moment cinétique (sans
// dissipation) et le bilan énergétique au cours de la simulation.
constexpr const char* FICHIER_DIAG   = "diagnostics.fich";

#endif
