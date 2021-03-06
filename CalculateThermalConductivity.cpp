#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>

#include "PhysicalConstants.hpp"
#include "SymmetryCoordinates.hpp"
#include "Interpolations.hpp"

const double temperature = 300;

int main(int argc, char **argv)
{
    std::vector<SymmetryCoordinates::CartesianPoint2<double>> qgrid;
    std::vector<std::array<double, NUMBER_OF_BRANCHES>> taugrid;
    std::vector<std::array<double, NUMBER_OF_BRANCHES>> velgrid;
    if (argc < 1) {
        std::cout << argv[0] << ": TAUFILE1 [TAUFILE2 ... TAUFILEN]\n";
        std::cout << "\tTAUFILE contains 8 space-separated values; qx qy tauinv1 ... tauinv6\n";
        return 0;
    }

    for (int i = 1; i < argc; i++) {
        std::ifstream ifs(argv[i]);
        SymmetryCoordinates::CartesianPoint2<double> cp;
        ifs >> cp.x;
        ifs >> cp.y;
        std::array<double, NUMBER_OF_BRANCHES> tau_branch;
        bool add_branch = true;
        for (size_t j = 0; j < tau_branch.size(); j++) {
            double tau_inv;
            ifs >> tau_inv;
            if (tau_inv < 1e-8) {
                std::cerr << "Warning, SMALL tau_inv for point (" << cp.x << ", " << cp.y << ").\n";
            }
            if (tau_inv == 0) {
                std::cerr << "Warning, ZERO tau_inv for point (" << cp.x << ", " << cp.y << "), probably need to refine qprime. Discarding q point.\n";
                add_branch = false;
            }
            if (std::isinf(tau_inv)) {
                tau_branch[j] = 0;
            } else {
                tau_branch[j] = 1.0 / tau_inv;
            }
        }

        std::array<double, NUMBER_OF_BRANCHES> vel_branch;
        for (size_t j = 0; j < vel_branch.size(); j++) {
            double vel;
            ifs >> vel;
            vel_branch[j] = vel;
        }
        if (add_branch) {
            qgrid.push_back(cp);
            taugrid.push_back(tau_branch);
            velgrid.push_back(vel_branch);
        }
    }

    double (*w_functions[])(double, double) = {
        w1, w2, w3, w4, w5, w6
    };

    // hack to estimate contribution from each q point (use area fraction for IBZ)
    const double hexagon_area = (3.0 * sqrt(3.0) / 2.0) * (2.0 / 3.0) * (2.0 / 3.0);
    const double dqA = hexagon_area / qgrid.size();

    double thermal_conductivity = 0;
    std::ofstream local_ofs("local_contributions");
    for (size_t i = 0; i < qgrid.size(); i++) {
        const SymmetryCoordinates::CartesianPoint2<double> &q = qgrid[i];
        const SymmetryCoordinates::SymmetryPoint2<double> sp = SymmetryCoordinates::fromCartesian(q);
        for (size_t j = 0; j < NUMBER_OF_BRANCHES; j++) {
            const double velocity = velgrid[i][j];

            const double ws = w_functions[j](sp.r, sp.t);
            const double f = DIRAC_CONSTANT * ws / (BOLTZMANN_CONSTANT * temperature);
            const double g = exp(f) / ((exp(f) - 1) * (exp(f) - 1));
            const double magq = std::sqrt(q.x*q.x + q.y*q.y);
            const double local_conductivity_value = DIRAC_CONSTANT*ws*velocity*taugrid[i][j]*g*magq;
            local_ofs << q.x << ' ' << q.y << ' ' << local_conductivity_value << '\n';
            thermal_conductivity += local_conductivity_value*dqA;
            // thermal_conductivity += ws*ws*taugrid[i][j]*f*velocity*velocity*dqA;
            if (thermal_conductivity < 0) {
                std::cout << "Thermal conductivity can't be negative.\n";
                std::cout << q.x << ' ' << q.y << '\n';
                std::cout << ws << ' ' << f << ' ' << g << ' ' << magq << ' ' << taugrid[i][j];
                return 0;
            }
        }
    }

    thermal_conductivity /= (4 * M_PI * BOLTZMANN_CONSTANT * temperature * temperature * GRAPHENE_THICKNESS);
    thermal_conductivity *= GRAPHENE_SAMPLE_LENGTH * GRAPHENE_SAMPLE_LENGTH;
    // thermal_conductivity *= (DIRAC_CONSTANT*DIRAC_CONSTANT) / (4 * M_PI * M_PI * BOLTZMANN_CONSTANT * temperature * temperature * GRAPHENE_THICKNESS);
    std::cout << thermal_conductivity;

    return 0;
}
