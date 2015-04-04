#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>

#include "PhysicalConstants.hpp"
#include "SymmetryCoordinates.hpp"
#include "Interpolations.hpp"

const double velocity = 22e3; // m/s
const double temperature = 300;

int main(int argc, char **argv)
{
    std::vector<SymmetryCoordinates::CartesianPoint2> qgrid;
    std::vector<std::array<double, NUMBER_OF_BRANCHES>> taugrid;
    if (argc < 1) {
        std::cout << argv[0] << ": TAUFILE1 [TAUFILE2 ... TAUFILEN]\n";
        std::cout << "\tTAUFILE contains 8 space-separated values; qx qy tau1 ... tau6\n";
        return 0;
    }

    for (size_t i = 1; i < argc; i++) {
        std::ifstream ifs(argv[i]);
        SymmetryCoordinates::CartesianPoint2 cp;
        ifs >> cp.x;
        ifs >> cp.y;
        qgrid.push_back(cp);
        std::array<double, NUMBER_OF_BRANCHES> tau_branch;
        for (size_t j = 0; j < tau_branch.size(); j++) {
            ifs >> tau_branch[j];
        }
        taugrid.push_back(tau_branch);
    }

    double (*w_functions[])(double, double) = {
        w1, w2, w3, w4, w5, w6
    };

    // hack to estimate contribution from each q point (use area fraction)
    const double hexagon_area = (3.0 * sqrt(3.0) / 2.0) * (2.0 / 3.0) * (2.0 / 3.0);
    const double dqA = hexagon_area / qgrid.size();

    double thermal_conductivity = 0;
    for (size_t i = 0; i < qgrid.size(); i++) {
        const SymmetryCoordinates::CartesianPoint2 &q = qgrid[i];
        const SymmetryCoordinates::SymmetryPoint2 sp = SymmetryCoordinates::fromCartesian(q);
        for (size_t j = 0; j < NUMBER_OF_BRANCHES; j++) {
            const double ws = w_functions[j](sp.r, sp.t);
            const double f = PLANCK_CONSTANT * ws / (BOLTZMANN_CONSTANT * temperature);
            const double g = exp(f) / ((exp(f) - 1) * (exp(f) - 1));
            const double magq = std::sqrt(q.x*q.x + q.y*q.y);
            thermal_conductivity += PLANCK_CONSTANT*ws*velocity*taugrid[i][j]*g*magq;
        }
    }

    thermal_conductivity /= (4 * M_PI * BOLTZMANN_CONSTANT * temperature * temperature * GRAPHENE_THICKNESS);
    std::cout << thermal_conductivity;

    return 0;
}
