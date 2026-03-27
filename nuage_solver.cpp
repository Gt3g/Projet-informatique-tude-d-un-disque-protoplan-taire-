#include "nuage_solver.h"
#include "forces.h"
#include "params.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <stdexcept>

// ================================================================
//  Représentation interne de l'état du système
//
//  Pour N particules, l'état est un vecteur de taille 6*N :
//    s = [ x0, y0, z0, vx0, vy0, vz0,
//          x1, y1, z1, vx1, vy1, vz1, ... ]
//
//  Les masses sont stockées dans le Nuage (invariantes).
// ================================================================
using State = std::vector<double>;

// ── Indices dans le vecteur d'état ──────────────────────────────
inline int idx_x (int i) { return 6*i + 0; }
inline int idx_y (int i) { return 6*i + 1; }
inline int idx_z (int i) { return 6*i + 2; }
inline int idx_vx(int i) { return 6*i + 3; }
inline int idx_vy(int i) { return 6*i + 4; }
inline int idx_vz(int i) { return 6*i + 5; }

// ================================================================
//  Conversion Nuage ↔ State
// ================================================================
static State nuage_to_state(const Nuage& nuage) {
    int N = (int)nuage.vecteur_de_part.size();
    State s(6 * N);
    for (int i = 0; i < N; ++i) {
        const Particule& p = nuage.vecteur_de_part[i];
        s[idx_x (i)] = p.coord_x;
        s[idx_y (i)] = p.coord_y;
        s[idx_z (i)] = p.coord_z;
        s[idx_vx(i)] = p.vit_x;
        s[idx_vy(i)] = p.vit_y;
        s[idx_vz(i)] = p.vit_z;
    }
    return s;
}

static void state_to_nuage(const State& s, Nuage& nuage) {
    int N = (int)nuage.vecteur_de_part.size();
    for (int i = 0; i < N; ++i) {
        Particule& p = nuage.vecteur_de_part[i];
        p.coord_x = s[idx_x (i)];
        p.coord_y = s[idx_y (i)];
        p.coord_z = s[idx_z (i)];
        p.vit_x   = s[idx_vx(i)];
        p.vit_y   = s[idx_vy(i)];
        p.vit_z   = s[idx_vz(i)];
    }
}

// ================================================================
//  Dérivée du système : ds/dt = derivee(s, nuage, params, t)
//
//  Pour chaque grain i (i ≠ étoile) :
//    d/dt [x,y,z]    = [vx, vy, vz]
//    d/dt [vx,vy,vz] = Σ accélérations gravitationnelles
//                     + pression PM        (si enable_pressure)
//                     + amortissement vertical (si t ≥ t_diss_start)
//                     + circularisation        (si t ≥ t_diss_start)
//
//  Gravitation :
//    - grain ← étoile : toujours actif
//    - grain ← grains : actif si params.enable_grain_grain = true
// ================================================================
static State derivee(const State& s,
                     const Nuage& nuage,
                     const SimParams& params,
                     double t) {

    const int N           = (int)nuage.vecteur_de_part.size();
    const int i_etoile    = 0;
    const double M_etoile = nuage.vecteur_de_part[i_etoile].masse;
    const bool diss_on    = (t >= params.t_diss_start);

    State dsdt(6 * N, 0.0);

    // ── Pression PM (une seule passe O(N + Ng³) avant la boucle) ──
    // Calcule les accélérations de pression pour tous les grains.
    // Activée en même temps que les dissipations (après t_diss_start).
    std::vector<double> ax_press(N, 0.0);
    std::vector<double> ay_press(N, 0.0);
    std::vector<double> az_press(N, 0.0);
    if (diss_on && params.enable_pressure)
        pression_grille(s, nuage, params, ax_press, ay_press, az_press);

    for (int i = 0; i < N; ++i) {

        // dx/dt = v
        dsdt[idx_x (i)] = s[idx_vx(i)];
        dsdt[idx_y (i)] = s[idx_vy(i)];
        dsdt[idx_z (i)] = s[idx_vz(i)];

        // L'étoile (i=0) est fixe : aucune force calculée sur elle
        if (i == i_etoile) continue;

        double ax = 0.0, ay = 0.0, az = 0.0;

        // ── Gravitation grain ← étoile (softened) ────────────────
        gravite_softened(i, i_etoile, s, nuage, params, ax, ay, az);

        // ── Gravitation grain ← grains (N-corps, O(N²)) ──────────
        // Activé si enable_grain_grain = true.
        // Physiquement nécessaire dès que M_disk / M_etoile > quelques %
        // (critère de Toomre, instabilités gravitationnelles).
        if (params.enable_grain_grain) {
            for (int j = 1; j < N; ++j) {
                if (j == i) continue;
                gravite_softened(i, j, s, nuage, params, ax, ay, az);
            }
        }

        // ── Dissipations + pression (activées après t_diss_start) ─
        // Ces forces modélisent les interactions gaz-disque.
        // Elles n'ont de sens physique qu'une fois un disque plan formé.
        if (diss_on) {
            // Amortissement vertical : dissipe v_z dans le plan du disque
            az += amortissement_vertical(s[idx_x(i)], s[idx_y(i)], s[idx_z(i)],
                                         s[idx_vz(i)], M_etoile, params);

            // Circularisation : amortit la vitesse radiale v_r → 0
            // Conserve L_z, circularise les orbites elliptiques
            circularisation(i, s, params, M_etoile, ax, ay);

            // Pression PM : accélérations pré-calculées sur la grille
            if (params.enable_pressure) {
                ax += ax_press[i];
                ay += ay_press[i];
                az += az_press[i];
            }
        }

        dsdt[idx_vx(i)] = ax;
        dsdt[idx_vy(i)] = ay;
        dsdt[idx_vz(i)] = az;
    }

    return dsdt;
}

// ================================================================
//  Opérations vectorielles sur State
// ================================================================
static State add(const State& a, const State& b) {
    State r(a.size());
    for (size_t k = 0; k < a.size(); ++k) r[k] = a[k] + b[k];
    return r;
}

static State scale(double c, const State& a) {
    State r(a.size());
    for (size_t k = 0; k < a.size(); ++k) r[k] = c * a[k];
    return r;
}

// ================================================================
//  Un pas RK4 d'ordre 4
// ================================================================
static State rk4_step(const State& s,
                       const Nuage& nuage,
                       const SimParams& params,
                       double t,
                       double h) {

    State k1 = derivee(s,                        nuage, params, t);
    State k2 = derivee(add(s, scale(h/2.0, k1)), nuage, params, t + h/2.0);
    State k3 = derivee(add(s, scale(h/2.0, k2)), nuage, params, t + h/2.0);
    State k4 = derivee(add(s, scale(h,     k3)), nuage, params, t + h);

    State s_next(s.size());
    for (size_t k = 0; k < s.size(); ++k)
        s_next[k] = s[k] + (h / 6.0) * (k1[k] + 2.0*k2[k] + 2.0*k3[k] + k4[k]);

    return s_next;
}

// ================================================================
//  Diagnostics physiques : E_cin, E_grav, E_tot, L_z
//
//  Ces quantités permettent de VALIDER le code :
//    - Sans dissipation : E_tot et L_z doivent être conservés
//    - Avec dissipation : L_z conservé (circularisation), E_tot décroît
//
//  E_cin  = Σ_i  (1/2) m_i v_i²
//  E_grav = Σ_{i>j} -G m_i m_j / r_ij_softened     [toutes paires]
//  L_z    = Σ_i  m_i (x_i vy_i - y_i vx_i)
// ================================================================
struct Diagnostics {
    double E_cin;
    double E_grav;
    double E_tot;
    double L_z;
};

static Diagnostics compute_diagnostics(const State& s, const Nuage& nuage,
                                        const SimParams& params) {
    const int N = (int)nuage.vecteur_de_part.size();

    double E_cin = 0.0, E_grav = 0.0, L_z = 0.0;

    // Énergie cinétique et moment cinétique
    for (int i = 0; i < N; ++i) {
        double mi = nuage.vecteur_de_part[i].masse;
        double vx = s[idx_vx(i)];
        double vy = s[idx_vy(i)];
        double vz = s[idx_vz(i)];
        double x  = s[idx_x(i)];
        double y  = s[idx_y(i)];

        E_cin += 0.5 * mi * (vx*vx + vy*vy + vz*vz);
        L_z   += mi * (x * vy - y * vx);
    }

    // Énergie potentielle gravitationnelle (paires i < j)
    for (int i = 0; i < N; ++i) {
        double mi = nuage.vecteur_de_part[i].masse;
        for (int j = i + 1; j < N; ++j) {
            // N'inclure grain-grain que si l'option est activée
            // (sinon l'énergie ne serait pas cohérente avec la dérivée)
            if (j != 0 && i != 0 && !params.enable_grain_grain) continue;

            double mj = nuage.vecteur_de_part[j].masse;
            double dx = s[idx_x(j)] - s[idx_x(i)];
            double dy = s[idx_y(j)] - s[idx_y(i)];
            double dz = s[idx_z(j)] - s[idx_z(i)];
            double r2 = dx*dx + dy*dy + dz*dz + params.eps * params.eps;
            E_grav   -= params.G * mi * mj / std::sqrt(r2);
        }
    }

    return { E_cin, E_grav, E_cin + E_grav, L_z };
}

// ================================================================
//  Écriture d'un pas de positions dans le fichier trajectoire
//  Format : t  i  x  y  z  vx  vy  vz
// ================================================================
static void write_step(std::fstream& fich, double t,
                        const State& s, int N) {
    fich << std::scientific << std::setprecision(8);
    for (int i = 0; i < N; ++i) {
        fich << std::setw(18) << t
             << std::setw(6)  << i
             << std::setw(20) << s[idx_x (i)]
             << std::setw(20) << s[idx_y (i)]
             << std::setw(20) << s[idx_z (i)]
             << std::setw(20) << s[idx_vx(i)]
             << std::setw(20) << s[idx_vy(i)]
             << std::setw(20) << s[idx_vz(i)]
             << "\n";
    }
}

// ================================================================
//  Écriture d'une ligne de diagnostics physiques
//  Format : t  E_cin  E_grav  E_tot  L_z  dLz_rel
// ================================================================
static void write_diag(std::fstream& fich, double t,
                        const Diagnostics& d, double Lz0) {
    double dLz_rel = (Lz0 != 0.0) ? (d.L_z - Lz0) / std::abs(Lz0) : 0.0;
    fich << std::scientific << std::setprecision(8)
         << std::setw(18) << t
         << std::setw(20) << d.E_cin
         << std::setw(20) << d.E_grav
         << std::setw(20) << d.E_tot
         << std::setw(20) << d.L_z
         << std::setw(16) << dLz_rel
         << "\n";
}

// ================================================================
//  Solveur RK4 principal
// ================================================================
void rk4_nuage(Nuage&             nuage,
               const SimParams&   params,
               double             t0,
               double             tEnd,
               double             h,
               const std::string& filename_traj,
               const std::string& filename_diag) {

    if (h <= 0.0)   throw std::invalid_argument("Le pas h doit être > 0");
    if (tEnd <= t0) throw std::invalid_argument("tEnd doit être > t0");

    int N = (int)nuage.vecteur_de_part.size();
    if (N == 0) throw std::runtime_error("Le nuage est vide");

    // ── Ouverture des fichiers ────────────────────────────────────
    std::fstream f_traj(filename_traj, std::ios::out | std::ios::trunc);
    if (!f_traj.is_open())
        throw std::runtime_error("Impossible d'ouvrir : " + filename_traj);

    std::fstream f_diag(filename_diag, std::ios::out | std::ios::trunc);
    if (!f_diag.is_open())
        throw std::runtime_error("Impossible d'ouvrir : " + filename_diag);

    // ── En-têtes ──────────────────────────────────────────────────
    f_traj << "# Evolution temporelle RK4 — nuage de " << N << " particules\n";
    f_traj << "# G=" << params.G << "  eps=" << params.eps
           << "  grain_grain=" << params.enable_grain_grain
           << "  t_diss=" << params.t_diss_start << " s\n";
    f_traj << "# " << std::string(96, '-') << "\n";
    f_traj << std::setw(18) << "t"
           << std::setw(6)  << "i"
           << std::setw(20) << "x"
           << std::setw(20) << "y"
           << std::setw(20) << "z"
           << std::setw(20) << "vx"
           << std::setw(20) << "vy"
           << std::setw(20) << "vz"
           << "\n"
           << "# " << std::string(96, '-') << "\n";

    f_diag << "# Diagnostics physiques — conservation de l'énergie et du moment cinétique\n";
    f_diag << "# dLz_rel = (Lz(t) - Lz(0)) / |Lz(0)|  : doit rester petit si L_z conservé\n";
    f_diag << "# " << std::string(96, '-') << "\n";
    f_diag << std::setw(18) << "t(s)"
           << std::setw(20) << "E_cin(J)"
           << std::setw(20) << "E_grav(J)"
           << std::setw(20) << "E_tot(J)"
           << std::setw(20) << "L_z(kg.m2/s)"
           << std::setw(16) << "dLz_rel"
           << "\n"
           << "# " << std::string(96, '-') << "\n";

    // ── Boucle de simulation ──────────────────────────────────────
    State s = nuage_to_state(nuage);
    double t = t0;
    int step_count = 0;

    // État initial : diagnostics de référence
    Diagnostics d0 = compute_diagnostics(s, nuage, params);
    double Lz0 = d0.L_z;

    write_step(f_traj, t, s, N);
    write_diag(f_diag, t, d0, Lz0);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "t=0 | E_tot=" << d0.E_tot << " J | L_z=" << Lz0 << " kg.m²/s\n";

    while (t < tEnd - 1e-12) {
        double step = std::min(h, tEnd - t);
        s = rk4_step(s, nuage, params, t, step);
        t += step;
        ++step_count;

        // Diagnostics
        if (step_count % params.pas_diag == 0) {
            Diagnostics d = compute_diagnostics(s, nuage, params);
            write_diag(f_diag, t, d, Lz0);

            std::cout << "t=" << t/AN << " ans"
                      << " | E_tot=" << d.E_tot
                      << " | dLz=" << (Lz0 != 0.0 ? (d.L_z - Lz0)/std::abs(Lz0) : 0.0)
                      << "\n";
        }

        // Positions
        if (step_count % params.pas_ecriture == 0)
            write_step(f_traj, t, s, N);
    }

    state_to_nuage(s, nuage);
    f_traj.close();
    f_diag.close();

    std::cout << "\nSimulation terminée : " << step_count << " pas, "
              << N << " particules\n";
    std::cout << "Trajectoires : " << filename_traj << "\n";
    std::cout << "Diagnostics  : " << filename_diag << "\n";
}
