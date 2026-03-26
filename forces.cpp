#include <cmath>
#include "forces.h"
#include "params.h"

// ── Indices dans le vecteur d'état ────────────────────────────────
static inline int idx_x (int i) { return 6*i + 0; }
static inline int idx_y (int i) { return 6*i + 1; }
static inline int idx_z (int i) { return 6*i + 2; }

// ================================================================
//  Gravitation softened de j_src sur i
// ================================================================
void gravite_softened(int i, int j_src,
                      const std::vector<double>& s,
                      const Nuage& nuage,
                      const SimParams& params,
                      double& ax, double& ay, double& az) {

    double mj = nuage.vecteur_de_part[j_src].masse;

    double dx = s[idx_x(j_src)] - s[idx_x(i)];
    double dy = s[idx_y(j_src)] - s[idx_y(i)];
    double dz = s[idx_z(j_src)] - s[idx_z(i)];

    double r2    = dx*dx + dy*dy + dz*dz + params.eps * params.eps;
    double r3    = r2 * std::sqrt(r2);
    double coeff = params.G * mj / r3;

    ax += coeff * dx;
    ay += coeff * dy;
    az += coeff * dz;
}

// ================================================================
//  Amortissement vertical dans le disque
// ================================================================
double amortissement_vertical(double xi, double yi, double zi,
                               double vzi,
                               double masse_etoile,
                               const SimParams& params) {

    double R_cyl = std::max(std::sqrt(xi*xi + yi*yi), params.eps);

    // Température locale : T(R) ∝ R^{-1/2}
    double T_loc = T0_DISK * std::sqrt(R0_DISK / R_cyl);

    // Vitesse du son locale : cs² = (R_gaz/mu) * T
    double cs2 = (R_GAZ_DISK / MU_MOL_DISK) * T_loc;

    // Fréquence képlérienne : Omega_K = sqrt(G*M / R³)
    double Omega_K2 = params.G * masse_etoile / (R_cyl * R_cyl * R_cyl);
    double Omega_K  = std::sqrt(Omega_K2);

    // Hauteur caractéristique du disque : H = cs / Omega_K
    double H = std::sqrt(cs2 / Omega_K2);

    // Amortissement actif seulement si la particule est dans le disque
    if (std::abs(zi) < H)
        return -(F_DISS_DISK * Omega_K) * vzi;

    return 0.0;
}

// ================================================================
//  Circularisation des orbites
//
//  Principe physique :
//    On décompose la vitesse dans le repère cylindrique local :
//
//      v = v_r * r̂_cyl  +  v_φ * φ̂_cyl  +  v_z * ẑ
//
//    La composante radiale v_r traduit l'excentricité de l'orbite.
//    On l'amortit avec un taux proportionnel à la fréquence
//    képlérienne locale :
//
//      F_e = -(f_circ * Omega_K) * v_r * r̂_cyl
//
//    Ce choix :
//      • Circularise progressivement l'orbite (v_r → 0)
//      • Conserve le moment cinétique L_z (v_φ inchangé)
//      • Echelle naturellement avec la dynamique locale (Omega_K)
//      • Donne un temps de circularisation t_e = 1/(f_circ * Omega_K)
//        ≈ (1/f_circ) / (2π) périodes orbitales
//
//  Note : la dissipation verticale (v_z) est gérée séparément
//  par amortissement_vertical().
// ================================================================
void circularisation(int i,
                     const std::vector<double>& s,
                     const Nuage& nuage,
                     const SimParams& params,
                     double masse_etoile,
                     double& ax, double& ay, double& az)
{
    double x  = s[idx_x(i)];
    double y  = s[idx_y(i)];
    double vx = s[6*i + 3];
    double vy = s[6*i + 4];

    // ── Rayon cylindrique (plan du disque) ────────────────────────
    double R_cyl2 = x*x + y*y;
    if (R_cyl2 < params.eps * params.eps) return;  // trop proche de l'axe
    double R_cyl = std::sqrt(R_cyl2);

    // ── Fréquence képlérienne locale ──────────────────────────────
    // Omega_K(R) = sqrt(G*M / R³)
    double Omega_K = std::sqrt(params.G * masse_etoile / (R_cyl2 * R_cyl));

    // ── Taux d'amortissement de l'excentricité ────────────────────
    // 1/t_e = f_circ * Omega_K  →  unité : s⁻¹  (physiquement cohérent)
    double t_e_inv = params.f_circ * Omega_K;

    // ── Vecteur radial cylindrique unitaire ───────────────────────
    // r̂_cyl = (x, y, 0) / R_cyl
    double er_x = x / R_cyl;
    double er_y = y / R_cyl;

    // ── Composante radiale de la vitesse ──────────────────────────
    // v_r = v · r̂_cyl  (positif = fuite, négatif = chute)
    double v_rad = vx * er_x + vy * er_y;

    // ── Force d'amortissement radial ──────────────────────────────
    // F_e = -(1/t_e) * v_r * r̂_cyl
    // La composante azimutale v_φ n'est pas touchée → L_z conservé

    if (R_cyl>RAYON_ETOILE*0.3){
    
    ax -= t_e_inv * v_rad * er_x;
    ay -= t_e_inv * v_rad * er_y;
    // az : délégué à amortissement_vertical()
    };
}
