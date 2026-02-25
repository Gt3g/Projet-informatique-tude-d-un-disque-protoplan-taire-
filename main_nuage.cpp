#include <iostream>
#include "nuage.h"
#include "nuage_solver.h"

// ============================================================
//  Simulation d'un nuage de particules :
//   - interactions gravitationnelles mutuelles
//   - force d'inertie selon z  (centrifuge + uniforme)
// ============================================================

int main() {

    // --------------------------------------------------------
    //  Création et initialisation du nuage
    //  init_random(n, R_e, m_e, m_p, q_e, q_p, L)
    //    n   : nombre de particules autour de l'électron central
    //    R_e : rayon d'exclusion autour de l'origine
    //    m_e : masse de la particule centrale
    //    m_p : masse des particules du nuage
    //    q_e : charge de la particule centrale
    //    q_p : charge des particules du nuage
    //    L   : demi-côté du cube de tirage
    // --------------------------------------------------------
    Nuage nuage;
    nuage.init_random(
        10,        // n  : 10 particules autour du centre
        0.5,       // R_e: rayon d'exclusion = 0.5 m
        1.0e10,    // m_e: masse centrale (kg)
        1.0e8,     // m_p: masse des particules (kg)
        0.0,       // q_e: charge centrale (non utilisée ici)
        0.0,       // q_p: charge particules (non utilisée ici)
        5.0        // L  : tirage dans [-5, 5] m
    );

    std::cout << "Nuage initial (" << nuage.vecteur_de_part.size()
              << " particules) :\n";
    nuage.affiche_nuage();

    // --------------------------------------------------------
    //  Paramètres physiques
    // --------------------------------------------------------
    SimParams params;
    params.G     = 6.674e-11;  // constante gravitationnelle SI
    params.eps   = 0.1;        // softening (évite r=0)
    params.omega = 0.5;        // pulsation force centrifuge selon z (rad/s)
    params.a_z   = -1.0;       // accélération uniforme selon z (m/s²)

    // --------------------------------------------------------
    //  Lancement de la simulation RK4
    // --------------------------------------------------------
    rk4_nuage(
        nuage,
        params,
        0.0,                    // t0   : temps initial (s)
        50.0,                   // tEnd : temps final   (s)
        0.5,                    // h    : pas de temps  (s)
        "evolution_nuage.fich"  // fichier de sortie
    );

    std::cout << "\nÉtat final du nuage :\n";
    nuage.affiche_nuage();

    return 0;
}
