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

// ── Constantes physiques fondamentales ──────────────────────────
constexpr double G_PHYS   = 6.674e-11;  // Constante gravitationnelle (m³ kg⁻¹ s⁻²)
constexpr double UA       = 1.496e11;   // Unité Astronomique (m)
constexpr double AN       = 3.156e7;    // Année (s)
constexpr double M_SOLEIL = 1.989e30;   // Masse solaire (kg)

// ── Paramètres de l'étoile centrale ─────────────────────────────
constexpr double MASSE_ETOILE  = M_SOLEIL;      // Masse de l'étoile (kg)
constexpr double RAYON_ETOILE  = 50.0  * UA;    // Rayon d'exclusion (m)

// ── Paramètres du nuage moléculaire ─────────────────────────────
constexpr int    N_PARTICULES  = 500;
// Rayon typique d'un cœur dense protostellaire : 50–200 UA
// (Motte & André 2001 ; André et al. 2014)
constexpr double RAYON_NUAGE   = 120.0 * UA;

// Masse totale du disque / nuage :
//   fraction disque typique f_disk = 0.05–0.20 de la masse stellaire
//   (Williams & Cieza 2011, Annu. Rev. Astron. Astrophys.)
//   On prend f_disk = 0.10 → M_disk = 0.10 M_sun
//
// CORRECTION de l'ancienne valeur MASSE_PART = 0.001 kg :
//   → les grains étaient 10^33 fois trop légers (masse négligeable)
//   → pas d'auto-gravité, pas de disque massif → sans valeur physique
constexpr double FRACTION_DISK = 0.10;
constexpr double MASSE_DISK    = FRACTION_DISK * MASSE_ETOILE;   // ~2e29 kg
constexpr double MASSE_PART    = MASSE_DISK / N_PARTICULES;       // ~4e26 kg/grain

constexpr double RAYON_PART    = 0.001;  // rayon numérique (non physique)

// ── Paramètre de rotation ────────────────────────────────────────
// beta = E_rot / |E_grav|
//
// Valeurs observées dans les cœurs denses protostellaires :
//   beta ~ 0.01 – 0.07   (Goodman et al. 1993, ApJ 406 ;
//                          Caselli et al. 2002, ApJ 572)
//   valeur typique retenue : beta ~ 0.04
//
// CORRECTION de l'ancienne valeur BETA_ROT = 0.40 :
//   → E_rot = 40% E_grav → nuage en quasi-équilibre centrifuge
//   → effondrement trop lent, disque trop étendu, non observationnel
constexpr double BETA_ROT = 0.5;

// ── Agitation thermique ──────────────────────────────────────────
// Vitesse du son dans H₂ à ~10 K : cs ~ 200 m/s
constexpr double V_THERM = 0.0;  // m/s (0 = agitation désactivée)

// ── Paramètres temporels ─────────────────────────────────────────
constexpr double DUREE_EN_TFF = 20.0;   // Durée totale (en unités de t_ff)
constexpr double PAS_EN_TFF   = 300.0;  // Nombre de pas par t_ff
constexpr int    PAS_ECRITURE = 5;      // Écriture positions toutes les N pas
constexpr int    PAS_DIAG     = 1;      // Écriture diagnostics tous les N pas

// ── Softening gravitationnel ─────────────────────────────────────
// eps = EPS_FACTOR * R_nuage / N^(1/3) (Dehnen 2001)
constexpr double EPS_FACTOR = 0.15;

// ── Auto-gravité grain ↔ grain ───────────────────────────────────
// true  : interactions N-corps complètes O(N²) — physiquement correct
// false : test-particle (grain ← étoile seulement) — rapide
//
// Avec MASSE_DISK = 0.10 M_sun, l'auto-gravité est non négligeable :
// critère de Toomre Q = cs*Omega_K / (pi*G*Sigma) ~ 1 dans le disque
// → instabilités gravitationnelles possibles si N-corps désactivé.
// Recommandation : true pour N ≤ 500, false pour tests rapides.
constexpr bool ENABLE_GRAIN_GRAIN = false;

// ── Activation différée des dissipations ─────────────────────────
// Les forces de dissipation (circularisation, amortissement vertical)
// ne sont physiques qu'une fois un disque formé (t > quelques t_ff).
// Les activer pendant l'effondrement sphérique (t < t_ff) biaiserait
// artificiellement la dynamique.
//
// Valeur retenue : t_diss = 2.0 t_ff
// (l'effondrement vers un disque prend typiquement 1–3 t_ff)
constexpr double T_DISS_START = 1.0;   // en unités de t_ff

// ── Paramètres du disque : amortissement vertical ────────────────
// Modèle disque mince irradié passif :
//   T(R) = T0 * (R0/R)^(1/2)   (Kenyon & Hartmann 1987)
//   H(R) = cs(R) / Omega_K(R)  ∝ R^(5/4)  (disque évasé)
constexpr double MU_MOL_DISK = 2.0e-3;  // Masse molaire H₂ (kg/mol)
constexpr double R_GAZ_DISK  = 8.314;   // Constante des gaz parfaits (J/mol/K)
constexpr double T0_DISK     = 150.0;   // Température à R₀ = 1 UA (K)
                                         // (Hayashi 1981 : T(1 UA) ≈ 150–280 K)
constexpr double R0_DISK     = UA;
// Coefficient d'amortissement vertical :
// f_diss ~ 1–3 (Fromang & Nelson 2006, A&A 457)
constexpr double F_DISS_DISK = 2.0;

// ── Paramètres du disque : circularisation des orbites ──────────
//
//  F_e = -(f_circ * Omega_K) * v_r * r̂_cyl
//  t_e = 1 / (f_circ * Omega_K)  ≈  (1/f_circ) / (2π)  périodes
//
//  f_circ = 0.01 → ~16 périodes pour circulariser une orbite
//  (Papaloizou & Larwood 2000 ; Tanaka & Ward 2004)
constexpr double F_CIRC = 1.1;

// ── Fichiers de sortie ────────────────────────────────────────────
constexpr const char* FICHIER_SORTIE = "evolution_nuage.fich";
// Fichier de diagnostics : t, E_cin, E_grav, E_tot, L_z, delta_Lz
// Permet de vérifier la conservation du moment cinétique (sans
// dissipation) et le bilan énergétique au cours de la simulation.
constexpr const char* FICHIER_DIAG   = "diagnostics.fich";

#endif
