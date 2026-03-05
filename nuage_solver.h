#ifndef NUAGE_SOLVER_H
#define NUAGE_SOLVER_H

#include "nuage.h"
#include <string>

constexpr double G_CONST  = 6.674e-11;
constexpr double EPS_SOFT = 1e-3;

struct SimParams {
    double G             = G_CONST;
    double eps           = EPS_SOFT;
    double omega         = 0.0;
    double a_z           = 0.0;
    int    pas_ecriture  = 1;
};

void rk4_nuage(Nuage&             nuage,
               const SimParams&   params,
               double             t0,
               double             tEnd,
               double             h,
               const std::string& filename);

#endif
