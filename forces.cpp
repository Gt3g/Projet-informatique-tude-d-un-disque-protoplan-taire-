#include <cmath>
#include <algorithm>
#include "forces.h"
#include "params.h"

// ── Indices dans le vecteur d'état ────────────────────────────────
static inline int idx_x (int i) { return 6*i + 0; }
static inline int idx_y (int i) { return 6*i + 1; }
static inline int idx_z (int i) { return 6*i + 2; }

// ================================================================
//  Pression par méthode Particle-Mesh (PM) — domaine cylindrique
//
//  La grille couvre uniquement le cylindre du disque :
//    x dans [-R_disk, R_disk],  y dans [-R_disk, R_disk],  z dans [-H_disk, H_disk]
//
//  Les particules hors du cylindre (R_cyl > R_disk ou |z| > H_disk)
//  sont ignorées : pression du gaz négligeable en dehors du disque.
//
//  Algorithme (4 étapes) :
//    1. Dépôt CIC    : masse des grains → grille rho[Ng³]
//    2. Pression     : P = K * rho^gamma sur chaque cellule
//    3. Gradient dP  : différences finies centrées (Neumann aux bords)
//    4. Interpolation: dP relue aux positions des grains, a = -dP/rho
//
//  Indice plat de la cellule (ix, iy, iz) :
//    flat = ix + Ng * iy + Ng² * iz
// ================================================================

// ── Indice plat ──────────────────────────────────────────────────
static inline int cell(int ix, int iy, int iz, int Ng) {
    return ix + Ng * iy + Ng * Ng * iz;
}

// ── Test d'appartenance au cylindre du disque ────────────────────
static inline bool dans_disque(double x, double y, double z,
                                double r_disk, double h_disk) {
    return (x*x + y*y <= r_disk * r_disk) && (std::abs(z) <= h_disk);
}

void pression_grille(const std::vector<double>& s,
                     const Nuage& nuage,
                     const SimParams& params,
                     std::vector<double>& ax_press,
                     std::vector<double>& ay_press,
                     std::vector<double>& az_press) {

    const int    N   = (int)nuage.vecteur_de_part.size();
    const int    Ng  = params.n_grid;
    const int    Ng3 = Ng * Ng * Ng;
    const double R   = params.r_disk_pm;
    const double H   = params.h_disk_pm;

    // ── Boîte de la grille = boîte englobante du cylindre ────────
    const double xmin = -R, xmax = R;
    const double ymin = -R, ymax = R;
    const double zmin = -H, zmax = H;

    const double dx = (xmax - xmin) / Ng;
    const double dy = (ymax - ymin) / Ng;
    const double dz = (zmax - zmin) / Ng;
    const double vol_cell = dx * dy * dz;

    // ── 1. Dépôt CIC ─────────────────────────────────────────────
    // L'étoile (i=0) est exclue : sa masse (M_etoile >> M_disk)
    // écraserait la densité du gaz et rendrait la pression sans sens.
    std::vector<double> rho_grid(Ng3, 0.0);

    for (int i = 1; i < N; ++i) {
        const double xi = s[idx_x(i)];
        const double yi = s[idx_y(i)];
        const double zi = s[idx_z(i)];

        if (!dans_disque(xi, yi, zi, R, H)) continue;

        const double mi = nuage.vecteur_de_part[i].masse;
        const double fx = (xi - xmin) / dx;
        const double fy = (yi - ymin) / dy;
        const double fz = (zi - zmin) / dz;

        const int ix = (int)fx;
        const int iy = (int)fy;
        const int iz = (int)fz;

        const double tx = fx - ix,  sx = 1.0 - tx;
        const double ty = fy - iy,  sy = 1.0 - ty;
        const double tz = fz - iz,  sz = 1.0 - tz;

        for (int diz = 0; diz <= 1; ++diz)
        for (int diy = 0; diy <= 1; ++diy)
        for (int dix = 0; dix <= 1; ++dix) {
            int cx = std::min(ix + dix, Ng - 1);
            int cy = std::min(iy + diy, Ng - 1);
            int cz = std::min(iz + diz, Ng - 1);
            double w = (dix ? tx : sx) * (diy ? ty : sy) * (diz ? tz : sz);
            rho_grid[cell(cx, cy, cz, Ng)] += mi * w / vol_cell;
        }
    }

    // ── 2. Pression P = K * rho^gamma ────────────────────────────
    const double cs2 = params.cs_pm * params.cs_pm;
    const double gam = params.gamma_pm;
    const double K   = (gam == 1.0) ? cs2 : cs2 / gam;

    std::vector<double> P_grid(Ng3);
    for (int k = 0; k < Ng3; ++k)
        P_grid[k] = K * std::pow(rho_grid[k], gam);

    // ── 3. Gradient de pression (différences finies centrées) ─────
    // Condition de Neumann aux bords (gradient nul).
    // En pratique les bords cylindriques sont dans des zones vides.
    std::vector<double> grad_x(Ng3, 0.0);
    std::vector<double> grad_y(Ng3, 0.0);
    std::vector<double> grad_z(Ng3, 0.0);

    for (int iz = 0; iz < Ng; ++iz)
    for (int iy = 0; iy < Ng; ++iy)
    for (int ix = 0; ix < Ng; ++ix) {
        int ixm = (ix > 0)      ? ix - 1 : ix;
        int ixp = (ix < Ng - 1) ? ix + 1 : ix;
        int iym = (iy > 0)      ? iy - 1 : iy;
        int iyp = (iy < Ng - 1) ? iy + 1 : iy;
        int izm = (iz > 0)      ? iz - 1 : iz;
        int izp = (iz < Ng - 1) ? iz + 1 : iz;

        // Dénominateur : 2h en intérieur, h aux bords (schéma unilatéral)
        double inv2dx = 1.0 / ((ixp - ixm) * dx);
        double inv2dy = 1.0 / ((iyp - iym) * dy);
        double inv2dz = 1.0 / ((izp - izm) * dz);

        int c0 = cell(ix, iy, iz, Ng);
        grad_x[c0] = (P_grid[cell(ixp,iy, iz, Ng)] - P_grid[cell(ixm,iy, iz, Ng)]) * inv2dx;
        grad_y[c0] = (P_grid[cell(ix, iyp,iz, Ng)] - P_grid[cell(ix, iym,iz, Ng)]) * inv2dy;
        grad_z[c0] = (P_grid[cell(ix, iy, izp,Ng)] - P_grid[cell(ix, iy, izm,Ng)]) * inv2dz;
    }

    // ── 4. Interpolation CIC inverse → accélération particule ─────
    for (int i = 1; i < N; ++i) {
        const double xi = s[idx_x(i)];
        const double yi = s[idx_y(i)];
        const double zi = s[idx_z(i)];

        if (!dans_disque(xi, yi, zi, R, H)) continue;

        const double fx = (xi - xmin) / dx;
        const double fy = (yi - ymin) / dy;
        const double fz = (zi - zmin) / dz;

        const int ix = (int)fx;
        const int iy = (int)fy;
        const int iz = (int)fz;

        const double tx = fx - ix,  sx = 1.0 - tx;
        const double ty = fy - iy,  sy = 1.0 - ty;
        const double tz = fz - iz,  sz = 1.0 - tz;

        double gx = 0.0, gy = 0.0, gz = 0.0, rho_i = 0.0;

        for (int diz = 0; diz <= 1; ++diz)
        for (int diy = 0; diy <= 1; ++diy)
        for (int dix = 0; dix <= 1; ++dix) {
            int cx = std::min(ix + dix, Ng - 1);
            int cy = std::min(iy + diy, Ng - 1);
            int cz = std::min(iz + diz, Ng - 1);
            double w = (dix ? tx : sx) * (diy ? ty : sy) * (diz ? tz : sz);
            int c0  = cell(cx, cy, cz, Ng);
            gx    += w * grad_x[c0];
            gy    += w * grad_y[c0];
            gz    += w * grad_z[c0];
            rho_i += w * rho_grid[c0];
        }

        // a = -(1/rho) * grad_P  (protégé contre densité nulle)
        if (rho_i > 0.0) {
            double inv_rho = 1.0 / rho_i;
            ax_press[i] -= inv_rho * gx;
            ay_press[i] -= inv_rho * gy;
            az_press[i] -= inv_rho * gz;
        }
    }
}

// ================================================================
//  Gravitation softened de j sur i
// ================================================================
void gravite_softened(int i, int j,
                      const std::vector<double>& s,
                      const Nuage& nuage,
                      const SimParams& params,
                      double& ax, double& ay, double& az) {

    double mj = nuage.vecteur_de_part[j].masse;

    double dx = s[idx_x(j)] - s[idx_x(i)];
    double dy = s[idx_y(j)] - s[idx_y(i)];
    double dz = s[idx_z(j)] - s[idx_z(i)];

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

    double T_loc = T0_DISK * std::sqrt(R0_DISK / R_cyl);
    double cs2   = (R_GAZ_DISK / MU_MOL_DISK) * T_loc;

    double Omega_K2 = params.G * masse_etoile / (R_cyl * R_cyl * R_cyl);
    double Omega_K  = std::sqrt(Omega_K2);
    double H_disk   = std::sqrt(cs2 / Omega_K2);

    if (std::abs(zi) < H_disk)
        return -(F_DISS_DISK * Omega_K) * vzi;

    return 0.0;
}

// ================================================================
//  Circularisation des orbites
//
//  Amortit la composante radiale de la vitesse (v_r → 0) sans
//  affecter la composante azimutale : le moment cinétique L_z
//  est ainsi conservé. La composante verticale v_z est traitée
//  séparément par amortissement_vertical().
// ================================================================
void circularisation(int i,
                     const std::vector<double>& s,
                     const SimParams& params,
                     double masse_etoile,
                     double& ax, double& ay)
{
    double x  = s[idx_x(i)];
    double y  = s[idx_y(i)];
    double vx = s[6*i + 3];
    double vy = s[6*i + 4];

    double R_cyl2 = x*x + y*y;
    if (R_cyl2 < params.eps * params.eps) return;
    double R_cyl = std::sqrt(R_cyl2);

    double Omega_K  = std::sqrt(params.G * masse_etoile / (R_cyl2 * R_cyl));
    double t_e_inv  = params.f_circ * Omega_K;

    double er_x = x / R_cyl;
    double er_y = y / R_cyl;
    double v_rad = vx * er_x + vy * er_y;

    if (R_cyl > RAYON_ETOILE * 0.3) {
        ax -= t_e_inv * v_rad * er_x;
        ay -= t_e_inv * v_rad * er_y;
    }
}
