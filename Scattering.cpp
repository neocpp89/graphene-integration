#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>

#include "fadiff.h"
#include "SymmetryCoordinates.hpp"
#include "PhysicalConstants.hpp"
#include "Interpolations.hpp"


const size_t POINTS_HINT = 4000000;

double delta(double arg, double sigma=10.0)
{
    return std::exp(arg*arg/(2*sigma*sigma));
}

template <typename Real>
Real PlanckDistribution(Real x)
{
    return 1.0 / (exp(x) - 1.0);
}

template <typename Real>
Real NSplitting(Real w, Real w_prime, Real temperature = 300)
{
    const Real factor = PLANCK_CONSTANT * 2 * M_PI / (BOLTZMANN_CONSTANT * temperature);
    Real x = w;
    Real x_prime = w_prime;
    return (1 + PlanckDistribution(factor*x) + PlanckDistribution(factor*(x - x_prime)));
}

template <typename Real>
Real NCombining(Real w, Real w_prime, Real temperature = 300)
{
    const Real factor = PLANCK_CONSTANT * 2 * M_PI / (BOLTZMANN_CONSTANT * temperature);
    Real x = w;
    Real x_prime = w_prime;
    return (PlanckDistribution(factor*x) - PlanckDistribution(factor*(x + x_prime)));
}

template <bool splitting>
std::vector<SymmetryCoordinates::CartesianPoint2> CalculateQPrime(const SymmetryCoordinates::CartesianPoint2 &q, const std::vector<SymmetryCoordinates::CartesianPoint2> &coord)
{
    std::vector<SymmetryCoordinates::CartesianPoint2> curve;
    curve.reserve(coord.size());

    if (splitting) {
        for (auto const &cp : coord) {
            SymmetryCoordinates::CartesianPoint2 diff = q - cp;
            SymmetryCoordinates::SymmetryPoint2 sp = SymmetryCoordinates::fromCartesian(diff);
            if (sp.r > 1) { 
                curve.push_back(cp);
            }
        }
    } else {
        //combining process
        for (auto const &cp : coord) {
            SymmetryCoordinates::CartesianPoint2 diff = q + cp;
            SymmetryCoordinates::SymmetryPoint2 sp = SymmetryCoordinates::fromCartesian(diff);
            if (sp.r > 1) { 
                curve.push_back(cp);
            }
        }
    }

    /* curve contains q' in cartesian coordinates. */
    return curve;
}

template <bool splitting>
std::vector<SymmetryCoordinates::CartesianPoint2> CalculateQDoublePrime(const SymmetryCoordinates::CartesianPoint2 &q, const std::vector<SymmetryCoordinates::CartesianPoint2> &qprime)
{
    std::vector<SymmetryCoordinates::CartesianPoint2> inp(qprime);

    for(auto &cp : inp) {
        auto const diff = q - cp;
        auto const sp = SymmetryCoordinates::fromCartesian(diff);

        /* inp contains q'' in cartesian coordinates. */
        cp.x = sp.r * std::tan(sp.t) / SymmetryCoordinates::L;
        cp.y = std::acos(std::cos(M_PI*sp.r)) / (M_PI * SymmetryCoordinates::L);
    }

    return inp;
}

template <bool splitting = true>
std::vector<double> IntegrateTauInverse(const SymmetryCoordinates::CartesianPoint2 &q, const std::vector<SymmetryCoordinates::CartesianPoint2> &coord)
{
    std::vector<SymmetryCoordinates::CartesianPoint2> curve = CalculateQPrime<splitting>(q, coord);
    std::vector<SymmetryCoordinates::CartesianPoint2> inp = CalculateQDoublePrime<splitting>(q, curve);

    const double hexagon_area = (3.0 * sqrt(3.0) / 2.0) * (2.0 / 3.0) * (2.0 / 3.0);
    const double dqA = hexagon_area / coord.size();
    const double dqx = sqrt(dqA) ;
    const double dqy = dqx;

    auto qsym = SymmetryCoordinates::fromCartesian(q);

    std::vector<SymmetryCoordinates::SymmetryPoint2> curvesym;
    curvesym.reserve(curve.size());

    std::vector<SymmetryCoordinates::SymmetryPoint2> inpsym;
    inpsym.reserve(curve.size());

    for (size_t i = 0; i < curve.size(); i++) {
        /*
            qsym, curvesym, and inpsym are symmetry coordinates of
            q, q', and q'' respectively.
        */
        curvesym.push_back(fromCartesian(curve[i]));
        inpsym.push_back(fromCartesian(inp[i]));
    }

    double (*w_functions[])(double, double) = {
        w1, w2, w3, w4, w5, w6
    };
    const size_t num_w = sizeof(w_functions)/sizeof(w_functions[0]);

    std::vector<double> tau_inv_branch(num_w, 0);

    const double tol = 1e11;
    for (size_t i = 0; i < curve.size(); i++) {
        for (size_t j = 0; j < num_w; j++) {
            const double w = 1e12 * w_functions[j](qsym.r, qsym.t);
            for (size_t k = 0; k < num_w; k++) {
                const double w_prime = 1e12 * w_functions[k](curvesym[i].r, curvesym[i].t);
                for (size_t l = 0; l < num_w; l++) {
                    const double w_doubleprime = 1e12 * w_functions[l](inpsym[i].r, inpsym[i].t);
                    if (splitting) {
                        double delta_w = w - w_prime - w_doubleprime;
                        if (std::fabs(delta_w) < tol) {
                            tau_inv_branch[j] += w*w_prime*w_doubleprime*dqA*NSplitting(w, w_prime);
                        }
                    } else {
                        double delta_w = w + w_prime - w_doubleprime;
                        if (std::fabs(delta_w) < tol) {
                            tau_inv_branch[j] += w*w_prime*w_doubleprime*dqA*NCombining(w, w_prime);
                        }
                    }

                    /*
                    std::cout << "Finished (" << j << ", " << k << ", " << l << ")";
                    if (scurve.size() > 0) {
                        std::cout << ' ' << scurve.size() << " points found.\n";
                    } else {
                        std::cout << ".\n";
                    }
                    */
                }
            }
        }
    }

    for (auto &tau_inv : tau_inv_branch) {
        tau_inv *= PLANCK_CONSTANT * (8 * M_PI * M_PI * M_PI / (GRAPHENE_DENSITY * LATTICE_A * LATTICE_A * LATTICE_A));
    }
    return tau_inv_branch;
}


int main(int argc, char **argv)
{
    if (argc < 4) {
        std::cout << argv[0] << ": GRIDFILE QX QY.\n";
        return 0;
    }
    std::string infile(argv[1]);
    SymmetryCoordinates::CartesianPoint2 q = {std::stod(argv[2]), std::stod(argv[3])};

    std::cout << "Q point is (" << q << ").\n";

    std::vector<SymmetryCoordinates::CartesianPoint2> coord;
    coord.reserve(POINTS_HINT);
    std::ifstream fin(argv[1]);
    while (true) {
        SymmetryCoordinates::CartesianPoint2 cp;
        fin >> cp.x;
        fin >> cp.y;
        if (fin.fail()) {
            break;
        }
        coord.push_back(cp);
    }
    std::cout << "Read " << coord.size() << " grid points.\n";

    auto tau_inv_branch_splitting = IntegrateTauInverse<true>(q, coord);
    auto tau_inv_branch_combining = IntegrateTauInverse<false>(q, coord);

    if (tau_inv_branch_combining.size() != tau_inv_branch_splitting.size()) {
        throw std::runtime_error("Size mismatch in processes.");
    }

    std::vector<double> tau_inv_branch(tau_inv_branch_splitting.begin(), tau_inv_branch_splitting.end());
    for (size_t i = 0; i < tau_inv_branch_combining.size(); i++) {
        tau_inv_branch[i] += tau_inv_branch_combining[i];
    }

    for (size_t i = 0; i < tau_inv_branch.size(); i++) {
        std::cout << "tau_inv_branch[" << i << "] = " << tau_inv_branch[i] << '\n';
    }

    std::string filename("output/q");
    filename = filename+std::to_string(q.x)+'_'+std::to_string(q.y);
    std::ofstream ofs(filename);
    ofs << q.x << ' ' << q.y;
    for (auto tau_inv: tau_inv_branch) {
        if (std::isinf(tau_inv)) {
            ofs << ' ' << 1e300;
        } else {
            ofs << ' ' << tau_inv;
        }
    }

    return 0;
}
