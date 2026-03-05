#include "nuage_solver.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <stdexcept>

// ============================================================
//  Représentation interne de l'état du système
//  Pour N particules : état = [ x0,y0,z0, vx0,vy0,vz0,
//                                x1,y1,z1, vx1,vy1,vz1, ... ]
//  Taille du vecteur : 6*N
// ============================================================
using State = std::vector<double>;

// ============================================================
//  Extraction de l'état depuis le Nuage
// ============================================================
static State nuage_to_state(const Nuage& nuage) {
    int N = nuage.vecteur_de_part.size();
    State s(6 * N);
    for (int i = 0; i < N; ++i) {
        const Particule& p = nuage.vecteur_de_part[i];
        s[6*i + 0] = p.coord_x;
        s[6*i + 1] = p.coord_y;
        s[6*i + 2] = p.coord_z;
        s[6*i + 3] = p.vit_x;
        s[6*i + 4] = p.vit_y;
        s[6*i + 5] = p.vit_z;
    }
    return s;
}

// ============================================================
//  Réinjection de l'état dans le Nuage
// ============================================================
static void state_to_nuage(const State& s, Nuage& nuage) {
    int N = nuage.vecteur_de_part.size();
    for (int i = 0; i < N; ++i) {
        Particule& p = nuage.vecteur_de_part[i];
        p.coord_x = s[6*i + 0];
        p.coord_y = s[6*i + 1];
        p.coord_z = s[6*i + 2];
        p.vit_x   = s[6*i + 3];
        p.vit_y   = s[6*i + 4];
        p.vit_z   = s[6*i + 5];
    }
}

// ============================================================
//  Dérivée du système : ds/dt = derivee(s, nuage, params)
//
//  Pour chaque particule i :
//    ds/dt[6i+0..2] = vitesse (vx, vy, vz)
//    ds/dt[6i+3..5] = accélération (grav + inertie z)
// ============================================================
static State derivee(const State& s,
                     const Nuage& nuage,
                     const SimParams& params) {

    int N = nuage.vecteur_de_part.size();
    State dsdt(6 * N, 0.0);

    for (int i = 0; i < N; ++i) {
        double xi = s[6*i + 0];
        double yi = s[6*i + 1];
        double zi = s[6*i + 2];
	double vz = s[6*i + 5];
        double mi = nuage.vecteur_de_part[i].masse;

        // dx/dt = vx, dy/dt = vy, dz/dt = vz
        dsdt[6*i + 0] = s[6*i + 3];
        dsdt[6*i + 1] = s[6*i + 4];
        dsdt[6*i + 2] = s[6*i + 5];

        // Accumulation des accélérations
        double ax = 0.0, ay = 0.0, az = 0.0;

        // ---- Gravitation newtonienne (softened) ----

	/*
        for (int j = 0; j < N; ++j) {      //pour toutes les interractions
	*/

	int j =0;        // seulement interraction avec l'étoile


	
            if (j == i) continue;
            double mj = nuage.vecteur_de_part[j].masse;

            double dx = s[6*j + 0] - xi;
            double dy = s[6*j + 1] - yi;
            double dz = s[6*j + 2] - zi;

            double r2    = dx*dx + dy*dy + dz*dz + params.eps * params.eps;
            double r3    = r2 * std::sqrt(r2);          // r² * r = r³
            double coeff_grav = params.G * mj / r3;

            ax += coeff_grav * dx;
            ay += coeff_grav * dy;
            az += coeff_grav * dz;

	    /*
        }                                    //pr tt les interractions
	    */
	

        //  Force centrifuge : l*l / dist_axe**3
	    
	double l = nuage.vecteur_de_part[i].moment_cin;
	double d_axe_carre = xi*xi + yi*yi;
	double coeff_cent = l*l / (d_axe_carre * d_axe_carre);

	//ax += coeff_cent * xi;
	//ay += coeff_cent * yi;

	    





        // ---- Force de pression verticale (disque déjà formé) ───────────
        //
        // Hypothèse disque mince : z ≪ R_cyl  → r ≈ R_cyl
        //
        // Gradient de pression isotherme :
        //   ρ(R,z) = ρ₀(R) · exp(-z²/2H²)
        //   a_z,press = -(1/ρ)·∂P/∂z = +c_s²·z/H²  = +Ω_K²(R)·z
        //
        // Profil de température irradiée : T(R) = T₀·(R₀/R)^{1/2}
        //   c_s²(R) = (R_gaz/μ)·T(R)
        //   H(R)    = c_s/Ω_K  ∝  R^{5/4}   → évasement (flaring) ✓
        //
        // Note : valide uniquement pour z ≪ R_cyl (disque formé).
        //        Pendant l'effondrement (z ~ R), ne pas activer ce terme.
        // ─────────────────────────────────────────────────────────────────

        constexpr double mu_mol  = 2.0e-3;   // kg/mol (H₂)
        constexpr double R_gaz   = 8.314;    // J/(mol·K)
        constexpr double T0_disk = 120.0;    // K à R₀ = 1 UA
        constexpr double R0_disk = 1.496e11; // 1 UA (m)

        double M_etoile  = nuage.vecteur_de_part[0].masse;
        double R_cyl_min = params.eps;
        double R_cyl     = std::max(std::sqrt(d_axe_carre), R_cyl_min);

        // T(R) ∝ R^{-1/2}
        double T_loc   = T0_disk * std::sqrt(R0_disk / R_cyl);

        // c_s²(R) = (R_gaz/μ)·T
        double cs2     = (R_gaz / mu_mol) * T_loc;

        // Ω_K²(R) = G·M / R_cyl³
        double Omega_K2 = params.G * M_etoile / (R_cyl * R_cyl * R_cyl);
        double Omega_K  = std::sqrt(Omega_K2);

        // H(R) = c_s / Ω_K
        double H       = std::sqrt(cs2 / Omega_K2);

        // a_z,press = +c_s²·z/H²  =  +Ω_K²·z
        // (s'oppose à la gravité verticale ; équilibre atteint à |z| ~ H(R))
        //az += cs2 * 1/(zi / (H * H))/100000;

        // ---- Dissipation verticale résiduelle ───────────────────────────
        // Amortit les oscillations autour de z=0, assure la convergence
        // vers l'équilibre hydrostatique H(R) ∝ R^{5/4}.
        constexpr double f_diss = 0.1;
        az += -(f_diss * Omega_K) * vz;


	

       //Collision dissipation de vz

        
	//	double gamma_z = 1e-9;
	//az += -gamma_z * vz;





	
	
        dsdt[6*i + 3] = ax;
        dsdt[6*i + 4] = ay;
        dsdt[6*i + 5] = az;
    }














    
    return dsdt;
}

// ============================================================
//  Addition et multiplication scalaire sur State
// ============================================================
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

// ============================================================
//  Un pas RK4
// ============================================================
static State rk4_step(const State& s,
                       const Nuage& nuage,
                       const SimParams& params,
                       double h) {
    // k1 évalué à t
    State k1 = derivee(s, nuage, params);

    // k2 évalué à t + h/2, état s + h/2*k1
    Nuage tmp = nuage;
    state_to_nuage(add(s, scale(h/2.0, k1)), tmp);
    State k2 = derivee(add(s, scale(h/2.0, k1)), tmp, params);

    // k3 évalué à t + h/2, état s + h/2*k2
    state_to_nuage(add(s, scale(h/2.0, k2)), tmp);
    State k3 = derivee(add(s, scale(h/2.0, k2)), tmp, params);

    // k4 évalué à t + h, état s + h*k3
    state_to_nuage(add(s, scale(h, k3)), tmp);
    State k4 = derivee(add(s, scale(h, k3)), tmp, params);

    // Combinaison finale
    State s_next(s.size());
    for (size_t k = 0; k < s.size(); ++k)
        s_next[k] = s[k] + (h / 6.0) * (k1[k] + 2.0*k2[k] + 2.0*k3[k] + k4[k]);

    return s_next;
}

// ============================================================
//  Écriture d'un pas dans le fichier
// ============================================================
static void write_step(std::fstream& fich, double t,
                        const State& s, int N) {
    fich << std::scientific << std::setprecision(8);
    for (int i = 0; i < N; ++i) {
        fich << std::setw(18) << t
             << std::setw(6)  << i
             << std::setw(20) << s[6*i + 0]   // x
             << std::setw(20) << s[6*i + 1]   // y
             << std::setw(20) << s[6*i + 2]   // z
	  << std::setw(20) << s[6*i + 3]   // vx
	  << std::setw(20) << s[6*i + 4]   // vy
	  << std::setw(20) << s[6*i + 5]   // vz
             << "\n";
    }
}

// ============================================================
//  Solveur RK4 principal
// ============================================================
void rk4_nuage(Nuage&             nuage,
               const SimParams&   params,
               double             t0,
               double             tEnd,
               double             h,
               const std::string& filename) {

    if (h <= 0.0)   throw std::invalid_argument("Le pas h doit être > 0");
    if (tEnd <= t0) throw std::invalid_argument("tEnd doit être > t0");

    int N = nuage.vecteur_de_part.size();
    if (N == 0) throw std::runtime_error("Le nuage est vide");

    // Ouverture du fichier
    std::fstream fich(filename, std::ios::out | std::ios::trunc);
    if (!fich.is_open())
        throw std::runtime_error("Impossible d'ouvrir : " + filename);

    // En-tête
    fich << "# Evolution temporelle RK4 — nuage de " << N << " particules\n";
    fich << "# G=" << params.G << "  eps=" << params.eps
         << "  omega=" << params.omega << "  a_z=" << params.a_z << "\n";
    fich << "# ---------------------------------------------------------------\n";
    fich << std::setw(14) << "t"
         << std::setw(6)  << "i"
         << std::setw(16) << "x"
         << std::setw(16) << "y"
         << std::setw(16) << "z"
      << std::setw(16) << "vx"
      << std::setw(16) << "vy"
      << std::setw(16) << "vz"
         << "\n";
    fich << std::string(100, '-') << "\n";

    // État initial
    State s = nuage_to_state(nuage);

    double t = t0;
    write_step(fich, t, s, N);

    int step_count = 0;

    while (t < tEnd - 1e-12) {
        double step = std::min(h, tEnd - t);

        // On met à jour le nuage avec l'état courant pour accéder aux masses/charges
        state_to_nuage(s, nuage);

        s = rk4_step(s, nuage, params, step);
        t += step;
        ++step_count;

        write_step(fich, t, s, N);
	std::cout<<t<<endl;
    }

    // Mise à jour finale du nuage
    state_to_nuage(s, nuage);

    fich.close();

    std::cout << "Simulation terminée : " << step_count << " pas, "
              << N << " particules\n";
    std::cout << "Fichier écrit       : " << filename << "\n";
}
